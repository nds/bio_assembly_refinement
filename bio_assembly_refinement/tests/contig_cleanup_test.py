import unittest
import filecmp
import os
from bio_assembly_refinement import contig_cleanup 
from pymummer import alignment

modules_dir = os.path.dirname(os.path.abspath(contig_cleanup.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestContigCleanup(unittest.TestCase):
	def test_contig_cleanup(self):
		'''Test steps of contig cleanup'''
		input_file = os.path.join(data_dir, 'CLEANUP_input_1.fa')
		skip_ids_file = os.path.join(data_dir, 'CLEANUP_ids_to_skip.txt')
		expected_output_filename = os.path.join(os.getcwd(),"filtered_CLEANUP_input_1.fa")
		expected_summary_filename = os.path.join(os.getcwd(),"contig_cleanup_summary.txt")
	
		tests = [
				[contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=6), 'CLEANUP_output_1.fa' ],
 				[contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=6, skip=skip_ids_file), 'CLEANUP_output_1_b.fa' ],
 				[contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=6, skip=['TEST_CONTIG_11']), 'CLEANUP_output_1_b.fa' ],
				]
				
		for t in tests:
			t[0].run()
			self.assertTrue(os.path.isfile(expected_output_filename))
			self.assertTrue(os.path.isfile(expected_summary_filename))
			self.assertTrue(filecmp.cmp(t[0].output_file, os.path.join(data_dir, t[1]) , shallow=False)) 
			os.remove(t[0].output_file)
			os.remove(t[0].summary_file)
		
        
		