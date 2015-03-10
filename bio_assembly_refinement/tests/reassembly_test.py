import unittest
import os
import shutil
from bio_assembly_refinement import reassembly 

modules_dir = os.path.dirname(os.path.abspath(reassembly.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestReassembly(unittest.TestCase):
	def test_reassembly(self):
		'''Test reassembly'''
		
		test_file = os.path.join(data_dir, "test_fasta_file.fa")
		output_file = os.path.join(os.getcwd(), 'reassembled_test_fasta_file.fa')
		summary_file = os.path.join(os.getcwd(), "smrtanalysis_summary.txt")
		
		reassembler = reassembly.Reassembly(input_file=test_file,
											read_data=data_dir,
											pacbio_exec=data_dir + "/dummy_pacbio_script",
											output_directory="reassembly"
											)
		
		reassembler.run()
		
		self.assertTrue(os.path.exists(output_file))
		self.assertTrue(os.path.exists(summary_file))
		
		os.remove(output_file)
		os.remove(summary_file)
# 		shutil.rmtree('reassembly')

		
		
