'''

Class to trim and circularise contigs with overlapping edges

Attributes:
-----------
dnaA_sequence : path to file with dnaA, refA, refB sequences (positional)
fasta_file : input fasta file
working_directory : path to working directory (default to current working directory)
contigs : dict of contigs (instead of fasta file)
alignments : pre computed alignments
dnaA_alignments : pre-computed alignments against dnaA (for testing)
overlap_offset: offset from edge that the overlap can start expressed as a % of length (default 49)
overlap_boundary_max : max boundary of overlap expressed as % of length of reference (default 50)
overlap_min_length : minimum length of overlap (default 2KB)
overlap_percent_identity : percent identity of match between ends (default 85)
dnaA_hit_percent_identity : percent identity of match to dnaA (default 80)
dnaA_hit_length_minimum : minimum acceptable hit length to dnaA expressed as % (of dnaA length) (default 95) 
summary_file :  summary file (default circularisation_summary_file.txt)
debug : do not delete temp files if set to true (default false)
			  
Sample usage:
-------------
from bio_assembly_refinement import circularisation

circulariser = circularisation.Circularisation(dnaA_sequence = dnaA_file,
	                                           fasta_file = myfile.fa		
											  )
circulariser.run()

Todo:
-----
1. Consider looking for and removing adaptor sequences for all contigs before circularising
2. Consider running promer to find dnaA as it may be more conserved at the protein level than at the sequence level
3. Extend logic to encompass more edge cases 

'''

import os
import re
from pyfastaq import tasks, sequences
from pyfastaq import utils as fastaqutils
from pymummer import alignment
from bio_assembly_refinement import utils, contig_history

class Circularisation:
	def __init__(self, 
				 dnaA_sequence,
				 fasta_file='file.fa', 
				 working_directory=None, 
				 contigs={},
				 alignments=[],
				 dnaA_alignments=[], # Can be used for testing 
				 overlap_offset=13, 
				 overlap_boundary_max=50, 
				 overlap_min_length=2000,
				 overlap_percent_identity=85,
				 dnaA_hit_percent_identity=80,
				 dnaA_hit_length_minimum=65,
				 summary_file = "circularisation_summary_file.txt",			  
				 debug=False):

		''' Constructor '''
		self.dnaA_sequence = dnaA_sequence
		self.fasta_file = fasta_file
		self.working_directory = working_directory if working_directory else os.getcwd()		
		self.contigs = contigs
		self.alignments = alignments
		self.dnaA_alignments = dnaA_alignments
		self.overlap_offset = overlap_offset * 0.01
		self.overlap_boundary_max = overlap_boundary_max * 0.01
		self.overlap_min_length = overlap_min_length
		self.overlap_percent_identity = overlap_percent_identity
		self.dnaA_hit_percent_identity = dnaA_hit_percent_identity
		self.dnaA_hit_length_minimum = dnaA_hit_length_minimum * 0.01	
		self.summary_file = summary_file
		self.output_file = self._build_final_filename()
		self.debug = debug
		
		# Extract contigs and generate nucmer hits if not provided
		if not self.contigs:
			self.contigs = {}
			tasks.file_to_dict(self.fasta_file, self.contigs) 
			
		# Keep track of what we do to contigs for the summary file	
		self.contig_histories = {}
		for id in self.contigs.keys():
			self.contig_histories[id] = contig_history.ContigHistory(original_id = id, original_length = len(self.contigs[id]))
		
		# Run nucmer
		if not self.alignments:
			self.alignments = utils.run_nucmer(self.fasta_file, self.fasta_file, self._build_alignments_filename(), min_percent_id=self.overlap_percent_identity)
		
		
	def _look_for_overlap_and_trim(self):
		''' Look for the (best) overlap in contigs. If found, trim overlap off the start. Remember contig for circularisation process '''		
# 		TODO: Optimise. Work this out when we parse alignments in clean contigs stage? 
		circularisable_contigs = []
		for contig_id in self.contigs.keys():
			original_sequence = self.contigs[contig_id]
			acceptable_offset = self.overlap_offset * len(original_sequence)
			boundary = self.overlap_boundary_max * len(original_sequence)
			best_overlap = None
			for algn in self.alignments:
				r_coords = [algn.ref_start, algn.ref_end]
				q_coords = [algn.qry_start, algn.qry_end]
				q_coords.sort() # Sometimes overlap can be inverted but we still want to think of them in a linear way
				if algn.qry_name == contig_id and \
				   algn.ref_name == contig_id and \
				   r_coords[0] < acceptable_offset and \
				   r_coords[1] < boundary and \
				   q_coords[0] > boundary and \
				   q_coords[1] > (algn.qry_length - acceptable_offset) and \
				   algn.hit_length_ref > self.overlap_min_length and \
				   algn.percent_identity > self.overlap_percent_identity:
					if not best_overlap or \
					   (r_coords[0] <= best_overlap.ref_start and \
					    q_coords[1] => best_overlap.qry_end ):
					   best_overlap = algn
					   best_overlap.qry_start = q_coords[0]
					   best_overlap.qry_end = q_coords[1]
				   
			if best_overlap:		
				best_q_coords = [best_overlap.qry_start, best_overlap.qry_end]
				best_q_coords.sort() # Sometimes overlap can be inverted
				self.contigs[contig_id] = original_sequence[best_overlap.ref_end+1:best_q_coords[1]+1]
				(self.contig_histories[contig_id]).overlap_length = best_overlap.hit_length_ref
				(self.contig_histories[contig_id]).overlap_location = str(best_overlap.ref_start) + "," + str(best_overlap.ref_end) + "-" + str(best_q_coords[0]) + "," + str(best_q_coords[1])
				circularisable_contigs.append(contig_id)				
		return circularisable_contigs  
		
		
	def _circularise_and_rename(self, contig_ids):
		'''
		Create a temporary multi FASTA file with circularisable contigs 
		Run nucmer with dnaA sequences 
		For each contig, circularise (either to start at dnaA or a random gene in the case of plasmids)
		Create a new name for the contigs
		'''		
		if not self.dnaA_alignments:
			self.dnaA_alignments = utils.run_nucmer(self._build_intermediate_filename(), self.dnaA_sequence, self._build_dnaA_alignments_filename(), min_percent_id=self.dnaA_hit_percent_identity)
		
		plasmid_count = 1
		chromosome_count = 1
		contig_ids.sort()
		 
		for contig_id in contig_ids:
			print(contig_id)
			plasmid = True		   		
			trimmed_sequence = self.contigs[contig_id]
			for algn in self.dnaA_alignments:	
				if algn.ref_name == contig_id and \
				   algn.hit_length_ref > (self.dnaA_hit_length_minimum * algn.qry_length) and \
				   algn.percent_identity > self.dnaA_hit_percent_identity:	     
					plasmid = False
					print(algn)
					if algn.on_same_strand():
						break_point = algn.ref_start						
					else:
						# Reverse complement sequence, circularise using new start of dnaA in the right orientation
						trimmed_sequence = trimmed_sequence.translate(str.maketrans("ATCGatcg","TAGCtagc"))[::-1]
						break_point = (algn.ref_length - algn.ref_start) - 1 #interbase
						(self.contig_histories[contig_id]).dnaA_on_reverse_strand = True

					self.contigs[contig_id] = trimmed_sequence[break_point:] + trimmed_sequence[0:break_point]		
					# Record history
					(self.contig_histories[contig_id]).new_name = 'chromosome' + str(chromosome_count)
					(self.contig_histories[contig_id]).location_of_dnaA = str(algn.ref_start) + "-" + str(algn.ref_end)
					chromosome_count += 1		
					break;
					
			if plasmid:
				# Choose random gene in plasmid, and circularise
				if len(self.contigs[contig_id]) > 20000:
					gene_start = utils.run_prodigal_and_get_start_of_a_gene(self.contigs[contig_id])
					if gene_start:
						self.contigs[contig_id] = trimmed_sequence[gene_start:] + trimmed_sequence[0:gene_start]	
						(self.contig_histories[contig_id]).location_of_gene_on_plasmid = gene_start
				(self.contig_histories[contig_id]).new_name = 'plasmid' + str(plasmid_count)
				plasmid_count += 1
		

	def _write_contigs_to_file(self, contig_ids, out_file):
		output_fw = fastaqutils.open_file_write(out_file)
		for id in contig_ids:
			new_name = (self.contig_histories[id]).new_name
			contig_name = new_name if new_name else id
			print(sequences.Fasta(contig_name, self.contigs[id]), file=output_fw)
		output_fw.close()
	
		
	def _produce_summary(self):
		'''Generate summary text and write to summary file'''
		text = 	'~~Circularisation~~\n'
		for id in self.contig_histories.keys():
			history = self.contig_histories[id]
			if history.overlap_length > 0:
				text += history.pretty_text()
			else:
				text += "Contig: " + id + " (Unable to circularise) \n"		 
		utils.write_text_to_file(text, self.summary_file)
		
			
	def _build_alignments_filename(self):
		return os.path.join(self.working_directory, "nucmer_all_contigs.coords")
		
		
	def _build_dnaA_alignments_filename(self):
		return os.path.join(self.working_directory, "nucmer_matches_to_dnaA.coords")
		
		
	def _build_intermediate_filename(self):
		return os.path.join(self.working_directory, "trimmed.fa")
		
			
	def _build_unsorted_circularised_filename(self):
		input_filename = os.path.basename(self.fasta_file)
		return os.path.join(self.working_directory, "unsorted_circularised_" + input_filename)	
		
		
	def _build_final_filename(self):
		input_filename = os.path.basename(self.fasta_file)
		return os.path.join(self.working_directory, "circularised_" + input_filename)	
			   
			   
	def run(self):	
		original_dir = os.getcwd()
		os.chdir(self.working_directory)	
		circularisable_contigs = self._look_for_overlap_and_trim()		
		self._write_contigs_to_file(self.contigs, self._build_intermediate_filename()) # Write trimmed sequences to file
		self._circularise_and_rename(circularisable_contigs)									
		self._write_contigs_to_file(circularisable_contigs, self._build_unsorted_circularised_filename()) # Write circularisable contigs to new file
		tasks.sort_by_size(self._build_unsorted_circularised_filename(), self.output_file) # Sort contigs in final file according to size
		
		self._produce_summary()
		
		if not self.debug:
			utils.delete(self._build_dnaA_alignments_filename())
			utils.delete(self._build_alignments_filename())
			utils.delete(self._build_intermediate_filename())
			utils.delete(self._build_unsorted_circularised_filename())
		
		os.chdir(original_dir)
