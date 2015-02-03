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
		
		# Expected values
		small_contigs = ['TEST_CONTIG_0']
		contained_contigs = ['TEST_CONTIG_11', 'TEST_CONTIG_12']
		expected_coords = [
			'\t'.join(['1',	'180',	'1',	'180',	'180',	'180',	'100.00',	'200',	'180',	'1',	'1',	'TEST_CONTIG_1',	'TEST_CONTIG_11']),
			'\t'.join(['1',	'180',	'1',	'180',	'180',	'180',	'100.00',	'180',	'200',	'1',	'1',	'TEST_CONTIG_11',	'TEST_CONTIG_1']),
			'\t'.join(['3',	'200',	'3',	'200',	'198',	'198',	'100.00',	'200',	'200',	'1',	'1',	'TEST_CONTIG_12',	'TEST_CONTIG_2']),
			'\t'.join(['3',	'200',	'3',	'200',	'198',	'198',	'100.00',	'200',	'200',	'1',	'1',	'TEST_CONTIG_2',	'TEST_CONTIG_12']),
		]
		expected_alignments = [alignment.Alignment(coord) for coord in expected_coords]
		expected_filtered_file = os.path.join(data_dir, 'test_fasta_file_filtered.fa')
				
		ccleaner = contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=6)
		self.assertEqual(ccleaner._find_small_contigs(), small_contigs) 
		self.assertEqual(ccleaner._find_contained_contigs(input_file), contained_contigs)
		ccleaner.run()
		self.assertTrue(os.path.isfile(output_file))# Does output file exist and is it named right? 
		self.assertTrue(filecmp.cmp(output_file, expected_filtered_file, shallow=False)) 
		os.remove(output_file)
		

