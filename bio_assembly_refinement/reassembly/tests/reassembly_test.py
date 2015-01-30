import unittest
import os
from bio_assembly_refinement.reassembly import reassembly 


modules_dir = os.path.dirname(os.path.abspath(reassembly.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestReassembly(unittest.TestCase):
	def test_reassembly(self):
		'''Test reassembly'''
		
		test_file = os.path.join(data_dir, "test_contigs.fa")
		
		reassembler = reassembly.Reassembly(input_file=test_file,
											read_data=data_dir,
											pacbio_exec=data_dir + "/dummy_pacbio_script",
											output_directory="reassembled"
											)
		
		reassembler.run()
		self.assertTrue(os.path.isdir("reassembly"))