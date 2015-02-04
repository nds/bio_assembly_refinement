import unittest
import os
from bio_assembly_refinement import reassembly 

modules_dir = os.path.dirname(os.path.abspath(reassembly.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestReassembly(unittest.TestCase):
	def test_reassembly(self):
		'''Test reassembly'''
		
		test_file = os.path.join(data_dir, "test_fasta_file.fa")
		actual_output = os.path.join(os.getcwd(), 'reassembled_test_fasta_file.fa')
		
		reassembler = reassembly.Reassembly(input_file=test_file,
											read_data=data_dir,
											pacbio_exec=data_dir + "/dummy_pacbio_script",
											output_directory="reassembly"
											)
		
		reassembler.run()
		
		self.assertTrue(os.path.exists(actual_output))
		os.remove(actual_output)
		
		
