import unittest
import filecmp
import os
from bio_assembly_refinement.contig_cleanup import contig_cleanup 
from pymummer import alignment

modules_dir = os.path.dirname(os.path.abspath(contig_cleanup.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestContigCleanup(unittest.TestCase):
	def test_find_small_contigs(self):
		'''Test _find_small_contigs'''
		input_file = os.path.join(data_dir, 'test_fasta_file_two.fa')
		expected_ids = ['test5', 'test4']
		ccleaner = contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=9)
		small_contigs_returned = ccleaner._find_small_contigs()
		self.assertEqual(small_contigs_returned.sort(), expected_ids.sort())


	def test_find_hits(self):
		'''Test _find_hits'''		
		expected_coords = [
			'\t'.join(['1',	'180',	'1',	'180',	'180',	'180',	'100.00',	'200',	'180',	'1',	'1',	'TEST_CONTIG_1',	'TEST_CONTIG_11']),
			'\t'.join(['1',	'180',	'1',	'180',	'180',	'180',	'100.00',	'180',	'200',	'1',	'1',	'TEST_CONTIG_11',	'TEST_CONTIG_1']),
			'\t'.join(['3',	'200',	'3',	'200',	'198',	'198',	'100.00',	'200',	'200',	'1',	'1',	'TEST_CONTIG_12',	'TEST_CONTIG_2']),
			'\t'.join(['3',	'200',	'3',	'200',	'198',	'198',	'100.00',	'200',	'200',	'1',	'1',	'TEST_CONTIG_2',	'TEST_CONTIG_12']),
		]
		expected_alignments = [alignment.Alignment(coord) for coord in expected_coords]
		input_file = os.path.join(data_dir, 'test_fasta_file_three.fa') # Has 12 contigs: contig 11 is contained in 1 completely, contig 12 matches 2 except for two bases
		ccleaner = contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=10)
		actual_alignments = ccleaner._find_hits()
		self.assertEqual(actual_alignments, expected_alignments)
		
	def test_find_contained_contigs(self):	
		'''Test _find_contained_contigs'''
		input_file = os.path.join(data_dir, 'test_fasta_file_three.fa')
		expected_ids = ['TEST_CONTIG_11', 'TEST_CONTIG_12']
		ccleaner = contig_cleanup.ContigCleanup(input_file)
		contained_contigs_returned = ccleaner._find_contained_contigs()
		self.assertEqual(contained_contigs_returned.sort(), expected_ids.sort())
		

	def test_run(self):
		'''Test run'''
		input_file = os.path.join(data_dir, 'test_fasta_file_three.fa')
		expected_output = os.path.join(data_dir, 'filtered_test_fasta_file_three.fa')
		ccleaner = contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=9, debug=True)
		ccleaner.run()
		actual_output = "tmp.filtered.fa"
		self.assertTrue(os.path.isfile(actual_output))# Does output file still exist and is it named right? 
		self.assertTrue(filecmp.cmp(actual_output, expected_output, shallow=False)) 
		os.remove("tmp.filtered.fa")

		ccleaner = contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=9)
		ccleaner.run()
		contig_ids_file = input_file + '.contig.ids.remove'
		self.assertFalse(os.path.isfile(contig_ids_file))# Has it cleaned up the ids file?
		
		another_input_file = os.path.join(data_dir, 'test_fasta_file_four.fa')
		ccleaner = contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=9) #Fasta file with contig 13 with deletion and rev comp bases at end
		ccleaner.run()
		self.assertTrue(filecmp.cmp(actual_output, expected_output, shallow=False)) 
		os.remove("tmp.filtered.fa")
