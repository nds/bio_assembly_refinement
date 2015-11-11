import unittest
import filecmp
import os
from bio_assembly_refinement import main 

modules_dir = os.path.dirname(os.path.abspath(main.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestMain(unittest.TestCase):
	def test_process_assembly(self):
		'''Test steps of contig cleanup '''
		input_file = os.path.join(data_dir, 'CLEANUP_input_1.fa')
		test_dnaA_file = os.path.join(data_dir, "test_dnaA.fa")
		intermediate_file = os.path.join(os.getcwd(), 'circularised_trimmed_filtered_CLEANUP_input_1.fa')
			
		processor = main.Main( fasta_file = input_file, 
						       dnaA_sequence = test_dnaA_file,
						       bax_files = data_dir,
						       pacbio_exec=data_dir + "/dummy_pacbio_script",		
						       debug = False				  
						)
		processor.process_assembly()
		self.assertTrue(os.path.exists(intermediate_file)) # i.e. before running Quiver
		os.remove(intermediate_file)
		
		# check for other summary files and clean them up
		os.remove(os.path.join(os.getcwd(), 'contig_cleanup_summary.txt'))
		

