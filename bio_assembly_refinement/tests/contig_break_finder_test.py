import unittest
import filecmp
import os
from bio_assembly_refinement import contig_break_finder 

modules_dir = os.path.dirname(os.path.abspath(contig_break_finder.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestContigBreakFinder(unittest.TestCase):
	def test_finding_dnaA(self):
		test_dnaA_file = os.path.join(data_dir, "test_dnaA.fa")
		test_fasta_file = os.path.join(data_dir, "test_fasta_file.fa")
		break_finder = contig_break_finder.ContigBreakFinder(fasta_file = test_fasta_file,
														 gene_file = test_dnaA_file								       
												        )	
		break_finder.run()
	
	