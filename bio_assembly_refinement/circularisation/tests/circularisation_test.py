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
		input_file = os.path.join(data_dir, 'Salmonella_typhi_CT18_pacbio_ordered.fa')
# 		expected_ids = ['test5', 'test4']
		circulariser = circularisation.Circularisation(input_file)
		circulariser.run()
		
		# Check circularisable step independently
		# Check if class can take fasta file and/or contigs and nucmer hits


