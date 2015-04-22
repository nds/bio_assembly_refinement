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
 				[contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=6, skip={'TEST_CONTIG_11'}), 'CLEANUP_output_1_b.fa' ],
				]
				
		for t in tests:
			t[0].run()
			self.assertTrue(os.path.isfile(expected_output_filename))
			self.assertTrue(os.path.isfile(expected_summary_filename))
			self.assertTrue(filecmp.cmp(t[0].output_file, os.path.join(data_dir, t[1]) , shallow=False)) 
			os.remove(t[0].output_file)
			os.remove(t[0].summary_file)
			
	def test_skipping_all(self):
		'''Test skipping all contigs'''
		
		input_file = os.path.join(data_dir, 'CLEANUP_input_1.fa')
		skip_ids_file_2 = os.path.join(data_dir, 'CLEANUP_ids_to_skip_2.txt')
		expected_summary_file = os.path.join(data_dir, 'CLEANUP_summary_file.txt')
		
		cleaner = contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=6, skip=skip_ids_file_2)
		cleaner.run()
		self.assertTrue(not os.path.isfile(cleaner.output_file))
		self.assertTrue(os.path.isfile(cleaner.summary_file))
		self.assertTrue(filecmp.cmp(cleaner.summary_file, expected_summary_file , shallow=False)) 
		os.remove(cleaner.summary_file)		
		
        
		