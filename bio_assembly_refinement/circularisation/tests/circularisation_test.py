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
		# Test data with 4 base overlap between ends of contig1 to simulate circularisable contig 
		test_hits = ['\t'.join(['1', '4', '60', '57', '4', '4', '100.00', '60', '60', '1', '1', 'contig1', 'contig1']),
				 	'\t'.join(['7', '9', '54', '52', '3', '3', '100.00', '60', '60', '1', '1', 'contig2', 'contig2'])
		]
		test_alignments = []
		for algn in test_hits:
			test_alignments.append(alignment.Alignment(algn))

		test_contigs = { "contig1": "TATCGAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTAT",
					 	 "contig2": "AGTCCACCGGGCACTGCAAGGTAAATTCTTACGCCCACTTTGTAGACCCTACCGTAAAGC"
				   		}
		# Does constructor take pre-computed results?	   
		circulariser = circularisation.Circularisation(contigs=test_contigs,
												   alignments = test_alignments,
												   offset = 7
												   )

		# Does circularisable check work?
		self.assertTrue(circulariser._circularisable("contig1"))
		self.assertFalse(circulariser._circularisable("contig2"))
		
		# Does constructor take fasta file?
		input_file = os.path.join(data_dir, "Salmonella_typhi_CT18_pacbio.fa")
		circulariser = circularisation.Circularisation(fasta_file=input_file)
		
	

