import unittest
import filecmp
import os
from bio_assembly_refinement.overlapping_contig_removal import overlapping_contig_removal 

modules_dir = os.path.dirname(os.path.abspath(overlapping_contig_removal.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestOverlappingContigRemoval(unittest.TestCase):

	def test_overlapping_contig_removal(self):
		'''Test overlapping_contig_removal'''
		input = os.path.join(data_dir, 'test_fasta_file.fa')
		expected_output = os.path.join(data_dir, 'filtered_test_fasta_file.fa')
		actual_output = 'tmp.filtered.fa'
		contig_ids_file = os.path.join(data_dir, 'remove_contig_ids.txt')
		overlapping_contig_removal.remove(input, contig_ids_file, actual_output)
		self.assertTrue(filecmp.cmp(actual_output, expected_output, shallow=False))
		os.unlink(actual_output)
