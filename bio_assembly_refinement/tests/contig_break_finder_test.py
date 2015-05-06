import unittest
import filecmp
import os
from bio_assembly_refinement import contig_break_finder 
from pyfastaq import tasks


modules_dir = os.path.dirname(os.path.abspath(contig_break_finder.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestContigBreakFinder(unittest.TestCase):
	def test_finding_dnaA(self):
		
		tests = [
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKING_input_0.fa"), 
							      											 gene_file = os.path.join(data_dir, "test_dnaA_1.fa"), 
							      											 choose_random_gene=False, 
							      											 rename = False,	
							      											 skip = os.path.join(data_dir, "BREAKING_skip_ids_0.txt")
							      											), 
			'BREAKING_input_0.fa' ], #skip 1, dnaA normal - nothing should change
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKING_input_0.fa"), 
							      											 gene_file = os.path.join(data_dir, "test_dnaA_1.fa"), 
							      											 choose_random_gene=False, 	
							      											 rename = False,
							      											 skip = os.path.join(data_dir, "BREAKING_skip_ids_0_all.txt")
							      											), 
			'BREAKING_input_0.fa' ], #skip all - - nothing should change
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKING_input_0.fa"), 
							      											 gene_file = os.path.join(data_dir, "test_dnaA_1.fa"),
							      											 rename = False, 
							      											 choose_random_gene=False, 	
							      											), 
			'BREAKING_output_0_all.fa' ], #skip none - test1 contig should be circularised
			
			]
				
		for t in tests:
			t[0].run()
			self.assertTrue(os.path.isfile(t[0].output_file))
			self.assertTrue(os.path.isfile(t[0].summary_file))
			# Read expected output file and compare sequences
			expected_contigs = {}
			tasks.file_to_dict(os.path.join(data_dir, t[1]), expected_contigs)
			for id in expected_contigs.keys():
# 				print(id + "\n")
# 				print("Expected: " + expected_contigs[id].seq)
# 				print("Got : " + t[0].contigs[id].seq)
				self.assertTrue(expected_contigs[id] == t[0].contigs[id])
			os.remove(t[0].output_file)
			os.remove(t[0].summary_file)	

# 	def test_skipping_all(self):
# 		''' Run prodigal in various situations '''

		
