''' 
Class to run pacbio_smrtanalysis RS_sequencing

Attributes:
-----------
input_file : input fasta file name
read_data : path to directory containing bax.h5 files
output_directory : directory for results (default reassembly)
working_directory : path to working directory (default to current working directory)
debug : do not delete temp files if set to true (default false)

Sample usage:
-------------

from bio_assembly_refinement import reassembly

reassembler = reassembly.Reassembly(input_file=my_file.fa,
									read_data=data_dir,
									)
reassembler.run()

'''
import os
import shutil
from pyfastaq import tasks, sequences
from pyfastaq import utils as fastaqutils
from bio_assembly_refinement import utils

class Reassembly:
	def __init__(self, 
				 input_file,
				 read_data,
				 output_directory="reassembly",
				 pacbio_exec="pacbio_smrtanalysis",				 
				 working_directory=None,
				 summary_file="smrtanalysis_summary.txt",
				 debug=False
				 ):
				 
		''' Constructor '''
		self.input_file = input_file
		self.read_data = read_data
		self.output_directory = output_directory
		self.pacbio_exec = pacbio_exec
		self.working_directory = working_directory if working_directory else os.getcwd()
		self.summary_file = summary_file 			
		self.debug = debug
		self.output_file = self._build_final_filename()
		
		
	def _build_final_filename(self):
		input_filename = os.path.basename(self.input_file)
		return os.path.join(self.working_directory, "reassembled_" + input_filename)	
		
		
	def _produce_summary(self, command):
		'''Generate summary text and write to summary file'''
		text = "~~smrtanalysis command~~ \n" + command + "\n"	
		utils.write_text_to_file(text, self.summary_file)
		
		
	def run(self):
		''' Run pacbio_smrtanalysis RS_sequencing'''
		original_dir = os.getcwd()
		os.chdir(self.working_directory)
			
		command = " ".join([self.pacbio_exec,
							'--reference', self.input_file,
							'RS_Resequencing',
							self.output_directory,
							self.read_data + "/*.bax.h5"
							])							
		fastaqutils.syscall(command)
# 		pacbio_smrtanalysis --reference /path/to/reference.fa RS_Resequencing Outputdir *.bax.h5

		# Move results file 
		default_results_file = os.path.join(self.working_directory, self.output_directory, "consensus.fasta")

		# TODO: Make it wait for bsub to finish before moving results file and delete reassembly folder
		if os.path.exists(default_results_file):
 			shutil.move(default_results_file, self.output_file)
			
		self._produce_summary(command)
		os.chdir(original_dir)
		
