import unittest
import filecmp
import os
from bio_assembly_refinement import contig_break_finder 

modules_dir = os.path.dirname(os.path.abspath(contig_break_finder.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestContigBreakFinder(unittest.TestCase):
	def test_finding_dnaA(self):
		
		tests = [
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "Providencia_rustigianii.fa"), 
							      gene_file = os.path.join(data_dir, "all_dnaA.fa"), 
							      ), 
							     'BREAKING_Providencia_rustigianni.fa' ],
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "Providencia_rustigianii.fa"),
                                                              gene_file = os.path.join(data_dir, "all_dnaA.fa"),
							      skip = os.path.join(data_dir, "BREAKING_skip_ids.txt"),		
                                                              ),
                                                             'BREAKING_Providencia_rustigianni_2.fa' ],
			]
				
		for t in tests:
			t[0].run()
#			self.assertTrue(os.path.isfile(expected_output_filename))
#			self.assertTrue(os.path.isfile(expected_summary_filename))
#			self.assertTrue(filecmp.cmp(t[0].output_file, os.path.join(data_dir, t[1]) , shallow=False)) 
#			os.remove(t[0].output_file)
#			os.remove(t[0].summary_file)		
		
		
	
	
