import unittest
import filecmp
import os
from bio_assembly_refinement import main 
from pymummer import alignment

modules_dir = os.path.dirname(os.path.abspath(main.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestMain(unittest.TestCase):
	def test_process_assembly(self):
		'''Test steps of contig cleanup '''
		input_file = os.path.join(data_dir, 'Salmonella_pacbio_unitig_0.fa')
		test_dnaA_file = os.path.join(data_dir, "dnaA.fa")
		output_file = os.path.join(os.getcwd(), 'reassembled_circularised_filtered_Salmonella_pacbio_unitig_0.fa')
			
		processor = main.Main( fasta_file = input_file, 
						       dnaA_sequence = test_dnaA_file,
						       bax_files = data_dir,
						       pacbio_exec=data_dir + "/dummy_pacbio_script",		
						       debug = False				  
						)
		processor.process_assembly()
		self.assertTrue(os.path.exists(output_file))
		os.remove(output_file)
		
		
