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
from fastaq import tasks
from pymummer import coords_file, alignment, nucmer
from bio_assembly_refinement.utils import utils

class Circularisation:
	def __init__(self, 
				 fasta_file=None, 
				 output_file="tmp.circularised.fa", 
				 contigs={},
				 alignments=[],
				 offset=12, #Acceptable offset from edge as a percentage
				 max_match_length=0.5, 
				 percent_identity=99,
				 dnaA_sequence=None,			  
				 debug=False):
				 
		''' Constructor '''
		self.fasta_file = fasta_file
		self.output_file = output_file
		self.contigs = contigs
		self.alignments = alignments
		self.offset = offset
		self.max_match_length = max_match_length
		self.percent_identity = percent_identity
		self.dnaA_sequence = "path to file"
		self.debug = debug
		
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
		Sample nucmer hit showing ends overlapping:
		1	2182	4783104	4780922	2182	2183	99.82	4791129	4791129	1	unitig_0|quiver	unitig_0|quiver
		
		TODO: 
		1. Optimise. We go through each alignment in contigcleanup. Can that information be re-used?
		2. Move this check to alignment class?
		'''
		
		for algn in self.alignments:
			acceptable_offset = (self.offset*0.01) * algn.ref_length		
			if algn.qry_name == contig_id and \
			   algn.ref_name == contig_id and \
			   algn.ref_start < acceptable_offset and \
			   algn.ref_end < (algn.ref_length * self.max_match_length) and \
			   algn.qry_start > (algn.qry_length * self.max_match_length) and \
			   algn.qry_end > (algn.qry_length - acceptable_offset) and \
			   algn.percent_identity > self.percent_identity:
				return True  
		return False
		
		
	def _circularise(self, contig_ids):
		'''
		Create a temporary FASTA file with contigs (choosing this as opposed to writing one contig to a file each time)
		Run nucmer with dnaA/refA/refB sequences 
		For each contig, split at site of match (if available), trim one overlapping end and re-join
		Return new contig sequence
		'''
		temp_file = "temp.contigs.fa"
		regex = "|".join(re.escape(id) for id in contig_ids)
		cmd = " ".join(["fastaq filter", self.fasta_file, temp_file, "--regex=\""+regex+"\""])
		print(cmd)
		utils.syscall(cmd)
		
				
				
	def run(self):
		sorted_contig_ids = sorted(self.contigs.keys()) #Sorting so that results are consistent in output file (for testing)
		circularisable_contigs = []
		for contig_id in sorted_contig_ids:
			if self._circularisable(contig_id):
 				print (contig_id)
 				circularisable_contigs.append(contig_id)
 				# Circularise and output contig
# 			else:
				# Just output contig as it is
			   
		self._circularise(circularisable_contigs)