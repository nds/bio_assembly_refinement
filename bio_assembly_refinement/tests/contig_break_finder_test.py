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
							      											 choose_random_gene=False, 	
							      											), 
			'BREAKING_Providencia_rustigianni.fa' ],
 			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "Providencia_rustigianii.fa"),
                                                                gene_file = os.path.join(data_dir, "all_dnaA.fa"),
							      								skip = os.path.join(data_dir, "BREAKING_skip_ids.txt"),	
							      								choose_random_gene=False, 	
                                                              ),
            'BREAKING_Providencia_rustigianni_2.fa' ],
 
			]
				
		for t in tests:
			t[0].run()
			self.assertTrue(os.path.isfile(t[0].output_file))
			self.assertTrue(os.path.isfile(t[0].summary_file))
#			self.assertTrue(filecmp.cmp(t[0].output_file, os.path.join(data_dir, t[1]) , shallow=False)) 
			os.remove(t[0].output_file)
			os.remove(t[0].summary_file)	

	def test_skipping_all(self):
			
		breaker = contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "Providencia_rustigianii.fa"),
                                                            gene_file = os.path.join(data_dir, "all_dnaA.fa"),
							      							skip = os.path.join(data_dir, "BREAKING_skip_ids_2.txt"),	
							      							choose_random_gene=False, 
							      							rename = False,	
                                                            )
		breaker.run()
		
		expected_summary_file = os.path.join(data_dir, "BREAKING_summary_file.txt")            
		self.assertTrue(os.path.isfile(breaker.output_file))
		self.assertTrue(os.path.isfile(breaker.summary_file))
# 		self.assertTrue(filecmp.cmp(breaker.summary_file, expected_summary_file , shallow=False)) 
# 		os.remove(breaker.summary_file)     
# 		os.remove(breaker.output_file)       
		
		