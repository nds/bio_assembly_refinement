'''
Take (or generate) a set of contigs and nucmer hits
For each contig, decide if it is circularisable
If it is, run nucmer against dnaA/repA/repB sequence
If match found, make that site the start, trim one overlapping end off contig and re-attach
Return file containing all contigs (some now circularised) 

TODO:
1. Consider looking for and removing adaptor sequences for all contigs before circularising
2. Consider running promer to find dnaA as it may be more conserved at the protein level than at the sequence level
3. Extend logic to encompass edge cases (current version only handles very basic, straight forward cases)

'''

import os
import re
from fastaq import tasks, sequences, utils
from pymummer import coords_file, alignment, nucmer

class Circularisation:
	def __init__(self, 
				 dnaA_sequence,
				 fasta_file=None, 
				 output_file="circularised.fa", 
				 contigs={},
				 alignments=[],
				 dnaA_alignments=[], # For testing 
				 offset=12, #Acceptable offset from edge as a percentage
				 max_match_length=0.5, 
				 percent_identity=99,			  
				 debug=False):
				 
		''' Constructor '''
		self.dnaA_sequence = dnaA_sequence
		self.fasta_file = fasta_file
		self.output_file = output_file
		self.contigs = contigs
		self.alignments = alignments
		self.dnaA_alignments = dnaA_alignments
		self.offset = offset
		self.max_match_length = max_match_length
		self.percent_identity = percent_identity
		self.debug = debug
		self.trim_values = {}
		
		# Extract contigs and generate nucmer hits if not provided
		if not self.contigs:
			self.contigs = {}
			tasks.file_to_dict(self.fasta_file, self.contigs) #Generate set of contigs
		
		if not self.alignments:
			results_file = str(self.fasta_file) + ".coords"
			runner = nucmer.Runner(self.fasta_file, self.fasta_file, results_file, coords_header=False, maxmatch=True) # nucmer default breaklength is 200
			runner.run()
			file_reader = coords_file.reader(results_file)
			self.alignments = [coord for coord in file_reader] #Generate set of alignments
		
		
	def _circularisable(self, contig_id):
		''' 
		Returns true if the ends of a contig overlap
		TODO: 
		1. Optimise. We go through each alignment in contigcleanup. Can that information be re-used?
		2. Move this check to pymmumer alignment class?
		'''	
		 
		for algn in self.alignments:
			acceptable_offset = (self.offset * 0.01) * algn.ref_length
			if algn.qry_name == contig_id and \
			   algn.ref_name == contig_id and \
			   algn.ref_start < acceptable_offset and \
			   algn.ref_end < (algn.ref_length * self.max_match_length) and \
			   algn.qry_start > (algn.qry_length * self.max_match_length) and \
			   algn.qry_end > (algn.qry_length - acceptable_offset) and \
			   algn.percent_identity > self.percent_identity:
				self.trim_values[contig_id] = round(algn.hit_length_ref/2)
				return True
		return False		
		  
		
	def _trim_and_circularise(self, contig_ids):
		'''
		Create a temporary FASTA file with circularisable contigs (choosing this as opposed to writing one contig to a file each time)
		Run nucmer with dnaA/refA/refB sequences 
		For each contig, split at site of match and re-join
		'''
		
		if not self.dnaA_alignments:
			temp_fasta_file = "temp.contigs.fa"
			self._write_contigs_to_file(contig_ids, temp_fasta_file)		
			results_file = "assembly_against_dnaA.coords"
			runner = nucmer.Runner(temp_fasta_file, self.dnaA_sequence, results_file, coords_header=False, maxmatch=True) 
			runner.run()
			file_reader = coords_file.reader(results_file)
			self.dnaA_alignments = [coord for coord in file_reader] 
		 
		for contig_id in contig_ids:			   		
			for algn in self.dnaA_alignments:	
				if algn.ref_name == contig_id and \
				   algn.hit_length_ref > 0.95*(algn.qry_length) and\
				   algn.percent_identity > 99:			       
					original_sequence = self.contigs[contig_id]
					# Trim
					cutoff = self.trim_values[contig_id]
					self.contigs[contig_id] = original_sequence[cutoff:len(original_sequence)-cutoff]
					# Circularise
					trimmed_sequence = self.contigs[contig_id]
					self.contigs[contig_id] = trimmed_sequence[(algn.ref_start-1-cutoff):] + trimmed_sequence[0:(algn.ref_start-1-cutoff)] 
					
	  
	def _write_contigs_to_file(self, contig_ids, out_file):
		output_fw = utils.open_file_write(out_file)
		for id in contig_ids:
			print(sequences.Fasta(id, self.contigs[id]), file=output_fw)
			   
			   
	def run(self):
		circularisable_contigs = []
		remaining_contigs = []
		for contig_id in self.contigs.keys():
			if self._circularisable(contig_id):
				circularisable_contigs.append(contig_id)
				self._trim_and_circularise(circularisable_contigs)
			else:
				remaining_contigs.append(contig_id)
								
		# Write all contigs to a file, ordered by size of contig (re-think)
		self._write_contigs_to_file(sorted(self.contigs, key=lambda id: len(self.contigs[id]), reverse=True), self.output_file)
	
#		self._write_contigs_to_file(circularisable_contigs + remaining_contigs, self.output_file)