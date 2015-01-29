import unittest
import filecmp
import os
from bio_assembly_refinement.circularisation import circularisation 
from pymummer import alignment

modules_dir = os.path.dirname(os.path.abspath(circularisation.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestCircularisation(unittest.TestCase):
	def test_circularisable(self):
		'''Test _circularisable'''
		
		test_dnaA_file = os.path.join(data_dir, "test_dnaA.fa")
		
		# Test data with 4 base overlap between ends of contig1 
		test_contigs = { "contig1": "TATCGAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTAT",
					 	 "contig2": "AGTCCACCGGGCACTGCAAGGTAAATTCTTACGCCCACTTTGTAGACCCTACCGTAAAGC"
				   		}
				   		
		test_hits = ['\t'.join(['1', '4', '60', '57', '4', '4', '100.00', '60', '60', '1', '1', 'contig1', 'contig1']),
				 	 '\t'.join(['7', '9', '54', '52', '3', '3', '100.00', '60', '60', '1', '1', 'contig2', 'contig2'])
		]
		test_alignments = []
		for algn in test_hits:
			test_alignments.append(alignment.Alignment(algn))
   		
		# Test dnaA: GGACCGCCTC
		test_dnaA_hits = ['\t'.join(['22', '31', '1', '10', '10', '10', '100.00', '60', '10', '1', '1', 'contig1', 'dnaA']),
		]
		test_dnaA_alignments = []
		for algn in test_dnaA_hits:
			test_dnaA_alignments.append(alignment.Alignment(algn))
		
		# Does constructor take pre-computed results?	   
		circulariser = circularisation.Circularisation(dnaA_sequence = test_dnaA_file,
													   contigs=test_contigs,
												       alignments = test_alignments,
												       dnaA_alignments = test_dnaA_alignments,
												       output_file="circularised.fa", 
												       offset = 7
												      )

		# Does circularisable check work?
		self.assertTrue(circulariser._circularisable("contig1"))
		self.assertFalse(circulariser._circularisable("contig2"))
		
		# Does trim and circularise work?
		circulariser.run()
		expected_output = os.path.join(data_dir, "circularised_test_contigs.fa")
		self.assertTrue(filecmp.cmp("circularised.fa", expected_output, shallow=False)) 
		
		# Does constructor take fasta file?
# 		input_file = os.path.join(data_dir, "test_contigs.fa")
# 		circulariser = circularisation.Circularisation(dnaA_sequence = test_dnaA_file, fasta_file=input_file, offset=7)
# 		circulariser.run()
		
	

