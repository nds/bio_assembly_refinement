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
		self.new_name = ''
		self.coverage = 0
		self.overlap_length = 0
		self.overlap_location = ''
		self.location_of_dnaA = 0
		self.dnaA_on_reverse_strand = False
		self.location_of_gene_on_plasmid = 0

					
	def pretty_text(self):
		return "Contig: " + self.original_id + "\n" + \
			   "\t" + "New name: " + self.new_name + "\n" + \
		       "\t" + "Original length: " + str(self.original_length) + "\n" + \
		       "\t" + "Coverage (not implemented yet): " + str(self.coverage) + "x" + "\n" + \
		       "\t" + "Length of overlap: " + str(self.overlap_length) + "(" + self.overlap_location + ")\n" + \
		       "\t" + "Location of dnaA on original contig: " + str(self.location_of_dnaA) + "\n" + \
		       "\t" + "Is dnA on reverse strand: " + str(self.dnaA_on_reverse_strand) + "\n" + \
			   "\t" + "Location of gene on plasmid (0 if not relevant): " + str(self.location_of_gene_on_plasmid) + "\n" 
			   