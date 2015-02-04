import unittest
import filecmp
import os
from bio_assembly_refinement import circularisation 
from pymummer import alignment

modules_dir = os.path.dirname(os.path.abspath(circularisation.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestCircularisation(unittest.TestCase):
	def test_circularisation_with_precomputed_values(self):
		'''Test _circularisation_with_precomputed_values'''		
		# Input data
		test_dnaA_file = os.path.join(data_dir, "test_dnaA.fa")
		test_contigs = { "contig1": "TATCGAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTAT", # 4 base overlap between ends of contig 1
					 	 "contig2": "AGTCCACCGGGCACTGCAAGGTAAATTCTTACGCCCACTTTGTAGACCCTACCGTAAAGC"
				   		}				   		
		test_overlap_coords = ['\t'.join(['1', '4', '60', '57', '4', '4', '100.00', '60', '60', '1', '1', 'contig1', 'contig1']),
				 	 		   '\t'.join(['7', '9', '54', '52', '3', '3', '100.00', '60', '60', '1', '1', 'contig2', 'contig2'])
							  ]
		test_overlap_alignments = [alignment.Alignment(coord) for coord in test_overlap_coords]   		
		# Test dnaA: GGACCGCCTC
		test_dnaA_coords = ['\t'.join(['19', '28', '1', '10', '10', '10', '100.00', '56', '10', '1', '1', 'contig1', 'dnaA'])] # Against trimmed contig 1 TCGAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCT
		test_dnaA_alignments = [alignment.Alignment(coord) for coord in test_dnaA_coords]
		
		# Output data
		expected_output = os.path.join(data_dir, "circularised_test_contigs.fa")
		actual_output = os.path.join(os.getcwd(), "circularised_file.fa")
	   	
		circulariser = circularisation.Circularisation(dnaA_sequence = test_dnaA_file,
													   contigs=test_contigs,
												       alignments = test_overlap_alignments,
												       dnaA_alignments = test_dnaA_alignments,
												       overlap_offset = 7,
												       
												      )
		# Does circularisable check work?
		self.assertTrue(circulariser._circularisable("contig1"))
		self.assertFalse(circulariser._circularisable("contig2"))
		# Does trim and circularise work?
		circulariser.run()
		self.assertTrue(filecmp.cmp(actual_output, expected_output, shallow=False)) 
		os.remove(actual_output)
		
		
	def test_circularisation_with_fasta_file(self):	
		# Does constructor take fasta file?
		input_file = os.path.join(data_dir, "Salmonella_pacbio_unitig_0.fa")
		dnaA_file = os.path.join(data_dir, "dnaA.fa")
		expected_output = os.path.join(data_dir, "Salmonella_pacbio_unitig_0_circularised.fa")
		actual_output = os.path.join(os.getcwd(), "circularised_Salmonella_pacbio_unitig_0.fa")
		
		circulariser = circularisation.Circularisation(dnaA_sequence = dnaA_file, fasta_file=input_file)
		circulariser.run()
		self.assertTrue(filecmp.cmp(actual_output, expected_output, shallow=False)) 
		os.remove(actual_output)