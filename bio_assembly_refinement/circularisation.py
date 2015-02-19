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
3. Extend logic to encompass edge cases (current version only handles very basic, straight forward cases)


'''

import os
import re
from pyfastaq import tasks, sequences
from pyfastaq import utils as fastaqutils
from pymummer import alignment
from bio_assembly_refinement import utils

class Circularisation:
	def __init__(self, 
				 dnaA_sequence,
				 fasta_file='file.fa', 
				 working_directory=None, 
				 contigs={},
				 alignments=[],
				 dnaA_alignments=[], # Can be used for testing 
				 overlap_offset=49, 
				 overlap_boundary_max=50, 
				 overlap_min_length=2000,
				 overlap_percent_identity=85,
				 dnaA_hit_percent_identity=80,
				 dnaA_hit_length_minimum=95,			  
				 debug=False):

		''' Constructor '''
		self.dnaA_sequence = dnaA_sequence
		self.fasta_file = fasta_file
		self.working_directory = working_directory		
		if not self.working_directory:
			self.working_directory = os.getcwd()		
		self.contigs = contigs
		self.alignments = alignments
		self.dnaA_alignments = dnaA_alignments
		self.overlap_offset = overlap_offset * 0.01
		self.overlap_boundary_max = overlap_boundary_max * 0.01
		self.overlap_min_length = overlap_min_length
		self.overlap_percent_identity = overlap_percent_identity
		self.dnaA_hit_percent_identity = dnaA_hit_percent_identity
		self.dnaA_hit_length_minimum = dnaA_hit_length_minimum * 0.01	
		self.debug = debug
		
		# Extract contigs and generate nucmer hits if not provided
		if not self.contigs:
			self.contigs = {}
			tasks.file_to_dict(self.fasta_file, self.contigs) 
		
		if not self.alignments:
			self.alignments = utils.run_nucmer(self.fasta_file, self.fasta_file, self._build_alignments_filename(), min_percent_id=self.overlap_percent_identity)
		
		self.output_file = self._build_final_filename()
		
		
	def _look_for_overlap_and_trim(self):
		''' Look for overlap in contigs. If found, trim overlap/2 off the ends. Remember contig for circularisation process '''		
# 		TODO: Optimise. Work this out when we parse alignments in clean contigs stage? Move check to pymummer?
		circularisable_contigs = []
		for contig_id in self.contigs.keys():
			acceptable_offset = self.overlap_offset * len(self.contigs[contig_id])
			boundary = self.overlap_boundary_max * len(self.contigs[contig_id])
			for algn in self.alignments:			
				if algn.qry_name == contig_id and \
				   algn.ref_name == contig_id and \
				   algn.ref_start < acceptable_offset and \
				   algn.ref_end < boundary and \
				   algn.qry_start > boundary and \
				   algn.qry_end > (algn.qry_length - acceptable_offset) and \
				   algn.hit_length_ref > self.overlap_min_length and \
				   algn.percent_identity > self.overlap_percent_identity:
					trim_value = round(algn.hit_length_ref/2)
					original_sequence = self.contigs[contig_id]
					self.contigs[contig_id] = original_sequence[trim_value:len(original_sequence)-trim_value]
					circularisable_contigs.append(contig_id)		
					break #Just find the biggest overlap from the end and skip any other hits
		return circularisable_contigs  
		
		
	def _circularise(self, contig_ids):
		'''
		Create a temporary multi FASTA file with circularisable contigs (choosing this as opposed to writing one contig to a file each time)
		Run nucmer with dnaA sequences 
		For each contig, circularise if possible
		'''
		
		if not self.dnaA_alignments:
			self.dnaA_alignments = utils.run_nucmer(self._build_intermediate_filename(), self.dnaA_sequence, self._build_dnaA_alignments_filename(), min_percent_id=self.dnaA_hit_percent_identity)
		 
		for contig_id in contig_ids:			   		
			for algn in self.dnaA_alignments:	
				if algn.ref_name == contig_id and \
				   algn.hit_length_ref > (self.dnaA_hit_length_minimum * algn.qry_length) and \
				   algn.percent_identity > self.dnaA_hit_percent_identity:			       
					trimmed_sequence = self.contigs[contig_id]
					self.contigs[contig_id] = trimmed_sequence[algn.ref_start:] + trimmed_sequence[0:algn.ref_start] 
					break;
	  
	def _write_contigs_to_file(self, contig_ids, out_file):
		output_fw = fastaqutils.open_file_write(out_file)
		for id in contig_ids:
			print(sequences.Fasta(id, self.contigs[id]), file=output_fw)
		output_fw.close()
			
			
	def get_contigs(self):
		return self.contigs
			
			
	def get_results_file(self):
		return self.output_file
		
			
	def _build_alignments_filename(self):
		return os.path.join(self.working_directory, "nucmer_all_contigs.coords")
		
		
	def _build_dnaA_alignments_filename(self):
		return os.path.join(self.working_directory, "nucmer_matches_to_dnaA.coords")
		
		
	def _build_intermediate_filename(self):
		return os.path.join(self.working_directory, "trimmed.fa")
		
			
	def _build_final_filename(self):
		input_filename = os.path.basename(self.fasta_file)
		return os.path.join(self.working_directory, "circularised_" + input_filename)	
			   
			   
	def run(self):
	
		original_dir = os.getcwd()
		os.chdir(self.working_directory)
	
		circularisable_contigs = self._look_for_overlap_and_trim()
		
		self._write_contigs_to_file(self.contigs, self._build_intermediate_filename()) # Write trimmed sequences to file
		self._circularise(circularisable_contigs)
								
		# Write all contigs to a file, ordered by size of contig (re-think. should contigs be re-named ti indicate possible chromosomes/plasmids?)
#		self._write_contigs_to_file(sorted(self.contigs, key=lambda id: len(self.contigs[id]), reverse=True), self.output_file)	
		
		self._write_contigs_to_file(circularisable_contigs, self.output_file) # Only write circularisable contigs (some will be chromosomes, some will be plasmids)
		
		if not self.debug:
			utils.delete(self._build_dnaA_alignments_filename())
			utils.delete(self._build_alignments_filename())
			utils.delete(self._build_intermediate_filename())
		
		os.chdir(original_dir)