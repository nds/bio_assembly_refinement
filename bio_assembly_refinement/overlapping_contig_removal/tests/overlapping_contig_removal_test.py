import unittest
import filecmp
import os
from bio_assembly_refinement.overlapping_contig_removal import overlapping_contig_removal 
from pymummer import alignment

modules_dir = os.path.dirname(os.path.abspath(overlapping_contig_removal.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestOverlappingContigRemoval(unittest.TestCase):

	def test_find_hits(self):
		'''Test find_hits'''
		expected_coords = [
            '\t'.join(['1',	'1000',	'1',	'1000',	'1000',	'1000',	'100.00',	'1000',	'1000',	'1',	'1',	'test_contig_1',	'test_contig_1',	'[IDENTITY]' ]),
            '\t'.join(['1',	'840',	'1',	'840',	'840',	'840',	'100.00',	'840',	'840',	'1',	'1',	'test_contig_2',	'test_contig_2',	'[IDENTITY]' ])
		]
		expected_alignments = [alignment.Alignment(coord) for coord in expected_coords]
		input_file = os.path.join(data_dir, 'test_fasta_file_one.fa')
		actual_alignments = overlapping_contig_removal.find_hits(input_file)
		self.assertEqual(actual_alignments, expected_alignments)
		
	def test_find_overlapping_contigs(self):	
		input_file = os.path.join(data_dir, 'test_fasta_file_one.fa')
		alignments = overlapping_contig_removal.find_hits(input_file)
		ids_to_keep = overlapping_contig_removal.find_overlapping_contigs(alignments)
		self.assertEquals(2,len(ids_to_keep))
		

	def test_overlapping_contig_removal(self):
		'''Test overlapping_contig_removal'''
		input = os.path.join(data_dir, 'test_fasta_file_two.fa')
		expected_output = os.path.join(data_dir, 'filtered_test_fasta_file_two.fa')
		actual_output = 'tmp.filtered.fa'
		contig_ids_file = os.path.join(data_dir, 'remove_contig_ids.txt')
		overlapping_contig_removal.remove_overlapping_contigs(input, contig_ids_file, actual_output)
		self.assertTrue(filecmp.cmp(actual_output, expected_output, shallow=False))
		os.unlink(actual_output)
		
		input_file = os.path.join(data_dir, 'test_fasta_file_one.fa')
		expected_output = os.path.join(data_dir, 'empty_fasta_file.fa')
		actual_output = 'tmp.filtered.fa'
		contig_ids_file = 'tmp.contid.ids.file'
		alignments = overlapping_contig_removal.find_hits(input_file)
		ids_to_keep = overlapping_contig_removal.find_overlapping_contigs(alignments)
		with open(contig_ids_file, mode='wt', encoding='utf-8') as ids_file:
			ids_file.write('\n'.join(ids_to_keep))
		overlapping_contig_removal.remove_overlapping_contigs(input_file, contig_ids_file, actual_output)
		self.assertTrue(filecmp.cmp(actual_output, expected_output, shallow=False))
		os.unlink(actual_output)
		os.unlink(contig_ids_file)
