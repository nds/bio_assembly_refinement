import unittest
import os
import shutil
import filecmp
from bio_assembly_refinement import reassembly 

modules_dir = os.path.dirname(os.path.abspath(reassembly.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestReassembly(unittest.TestCase):
	def test_reassembly(self):
		'''Test reassembly'''
		
		test_file = os.path.join(data_dir, "CLEANUP_input_1.fa")
		summary_file = os.path.join(os.getcwd(), "quiver_command_summary.txt")
		
		reassembler = reassembly.Reassembly(input_file=test_file,
											read_data=data_dir,
											pacbio_exec=data_dir + "/dummy_pacbio_script",
											)
		
		reassembler.run()
		self.assertTrue(os.path.exists(summary_file))
		os.remove(summary_file)
		
		another_test_file = os.path.join(data_dir, "empty_file.fa")
#		expected_summary_file = os.path.join(data_dir, "smrtanalysis_summary_empty_file.txt")
		reassembler = reassembly.Reassembly(input_file=another_test_file,
											read_data=data_dir,
											pacbio_exec=data_dir + "/dummy_pacbio_script",
											)
		
		reassembler.run()
#		self.assertTrue(filecmp.cmp(summary_file, expected_summary_file, shallow=False))
		os.remove(summary_file)
		
		shutil.rmtree(os.path.join(os.getcwd(), "improved_assembly"))

		
		
