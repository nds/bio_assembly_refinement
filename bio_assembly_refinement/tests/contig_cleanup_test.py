import unittest
import filecmp
import os
from bio_assembly_refinement import contig_cleanup 
from pymummer import alignment

modules_dir = os.path.dirname(os.path.abspath(contig_cleanup.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestContigCleanup(unittest.TestCase):
	def test_run(self):
		'''Test steps of contig cleanup '''
		input_file = os.path.join(data_dir, 'test_fasta_file.fa')
		output_file = os.path.join(os.getcwd(), 'filtered_test_fasta_file.fa')
		summary_file = os.path.join(os.getcwd(), 'contig_filtration_summary.txt')
		expected_filtered_file = os.path.join(data_dir, 'test_fasta_file_filtered.fa')
		expected_summary_file = os.path.join(data_dir, 'filtration_summary_file.txt')
				
		ccleaner = contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=6)
		ccleaner.run()
		self.assertTrue(os.path.isfile(output_file))# Does output file exist and is it named right? 
		self.assertTrue(filecmp.cmp(output_file, expected_filtered_file, shallow=False)) 
		self.assertTrue(os.path.isfile(summary_file))# Does summary file exist?
		os.remove(output_file)
		os.remove(summary_file)


	def test_run_with_ids_to_keep(self):
			'''Test running with a list of ids to keep '''
			input_file = os.path.join(data_dir, 'test_fasta_file.fa')
			output_file = os.path.join(os.getcwd(), 'filtered_test_fasta_file.fa')
			expected_filtered_file = os.path.join(data_dir, 'test_fasta_file_filtered_two.fa')
			summary_file = os.path.join(os.getcwd(), 'contig_filtration_summary.txt')
			ids_to_keep = ['TEST_CONTIG_11']
			ids_to_keep_file = os.path.join(data_dir, 'ids_to_keep.txt')
			ccleaner = contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=6, ids=ids_to_keep_file)
			ccleaner.run()
			self.assertTrue(os.path.isfile(output_file))# Does output file exist and is it named right? 
			self.assertTrue(filecmp.cmp(output_file, expected_filtered_file, shallow=False)) 
			os.remove(output_file)
			os.remove(summary_file)



