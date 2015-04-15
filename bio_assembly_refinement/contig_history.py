''' 
Class to record the processing done to a contig

Attributes:
-----------

Sample usage:
-------------

'''

import os
import shutil


class ContigHistory:
	def __init__(self, original_id, original_length):
				 
		''' Attributes '''
		self.original_id = original_id
		self.original_length = original_length
		self.new_length = '-'
		self.new_name = ''
		self.coverage = 0
		self.overlap_length = '-'
		self.overlap_location = 'No overlap found'
		self.position_used_for_circularisation = '-' #Start of dnaA in the case of chromosomes, or a predicted gene in the case of plasmids
		self.dnaA_or_gene_on_reverse_strand = False
		self.comment = '-'

					
	def pretty_text(self):
		'''Returns a string with information about contig'''
		# ID New_ID Original_length New_length Overlap Circularisation_point_in_trimmed_sequence dnaA_reversed Comment		
		return "\t".join([self.original_id,
						  self.new_name,
						  str(self.original_length),
						  str(self.new_length),
						  str(self.overlap_length)+"(" + self.overlap_location + ")",
						  str(self.position_used_for_circularisation),
						  str(self.dnaA_or_gene_on_reverse_strand),
						  self.comment,						  
						  ])
		
		
		
# 		return "Contig: " + self.original_id + "\n" + \
# 			   "\t" + "New name: " + self.new_name + "\n" + \
# 		       "\t" + "Original length: " + str(self.original_length) + "\n" + \
# 		       "\t" + "New length: " + str(self.new_length) + "\n" + \
# 		       "\t" + "Length and location of overlap: " + str(self.overlap_length) + " (" + self.overlap_location + ")\n" + \
# 		       "\t" + "Location of dnaA on original contig: " + str(self.location_of_dnaA) + "\n" + \
# 		       "\t" + "Is dnA on reverse strand: " + str(self.dnaA_on_reverse_strand) + "\n" + \
# 			   "\t" + "Location of gene on plasmid that was used as break point (0 if not relevant): " + str(self.location_of_gene_on_plasmid) + "\n" 
# 			   