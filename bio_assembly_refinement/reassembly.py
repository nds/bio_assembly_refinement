''' 
Class to run pacbio_smrtanalysis RS_sequencing

Attributes:
-----------
input_file : input fasta file name
read_data : path to directory containing bax.h5 files
output_directory : directory for results (default reassembly)
no_bsub : If set, will run quiver as a child process and not bsub it (useful for pipeline)
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
import time
import subprocess
from pyfastaq import tasks, sequences
from pyfastaq import utils as fastaqutils
from bio_assembly_refinement import utils

class Reassembly:
	def __init__(self, 
				 input_file,
				 read_data,
				 output_directory="improved_assembly",
				 pacbio_exec="pacbio_smrtanalysis",		
				 no_bsub=False,		 
				 working_directory=None,
				 summary_file="quiver_command_summary.txt",
				 debug=False
				 ):
				 
		''' Constructor '''
		self.input_file = input_file
		self.read_data = read_data
		self.output_directory = output_directory
		self.pacbio_exec = pacbio_exec
		self.no_bsub = no_bsub
		self.working_directory = working_directory if working_directory else os.getcwd()
		self.summary_file = summary_file 			
		self.debug = debug
		self.output_file = self._build_final_filename()
		
		
	def _build_final_filename(self):
		input_filename = os.path.basename(self.input_file)
		return os.path.join(self.working_directory, "reassembled_" + input_filename)
		
		
	def _build_default_filename(self):
		return os.path.join(self.working_directory, self.output_directory, "consensus.fasta")
		
		
	def _produce_summary(self, command):
		'''Generate summary text and write to summary file'''
		text = "~~smrtanalysis command~~ \n" + command + "\n"	
		utils.write_text_to_file(text, self.summary_file)
		
		
	def run(self):
		''' Run pacbio_smrtanalysis RS_sequencing'''
		original_dir = os.getcwd()
		os.chdir(self.working_directory)
		no_bsub_option = "--no_bsub" if self.no_bsub else ""
# 		pacbio_smrtanalysis --memory 6 --reference /path/to/reference.fa RS_Resequencing Outputdir *.bax.h5
		command = " ".join([self.pacbio_exec,
							"--memory 6",
							no_bsub_option,
							"--reference", self.input_file,
							"RS_Resequencing",
							self.output_directory,
							self.read_data + "/*.bax.h5"
							])	
							
		if(os.path.getsize(self.input_file)):						
			fastaqutils.syscall(command)
#			subprocess.call(command, shell=True)
			self._produce_summary(command)
		else:
			self._produce_summary("File empty: " + self.input_file)

		# Wait for bsub to finish. Look into using process.wait() here
# 		while not os.path.exists(self._build_default_filename()): 
# 			time.sleep(30) 
# 	
# 		# Move results file and produce summary 
# 		shutil.move(self._build_default_filename(), self._build_final_filename())

# 		We have abandoned this effort of trying to wait for the bsub quiver job to finish 
# 		When run in the pipeline, we can carry out this step as a separate task
# 		When run on the command line the user will just have to live with the quiver output (consensus.fasta)
			
		
		
		# Clean up quiver directory
# 		if not self.debug and os.path.exists(os.path.join(self.working_directory, self.output_directory)):
# 			shutil.rmtree(os.path.join(self.working_directory, self.output_directory))
			
		os.chdir(original_dir)
		
