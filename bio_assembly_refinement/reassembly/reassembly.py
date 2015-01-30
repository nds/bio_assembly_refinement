from fastaq import tasks, sequences, utils

class Reassembly:
	def __init__(self, 
				 input_file,
				 read_data,
				 pacbio_exec,				 
				 output_directory,
				 debug=False
				 ):
				 
		''' Constructor '''
		self.input_file = input_file
		self.read_data = read_data
		self.pacbio_exec = pacbio_exec
		self.output_directory = output_directory
		self.debug = debug
	
		
	def run(self):
		''' Run pacbio_smrtanalysis RS_sequencing'''	
		command = " ".join([self.pacbio_exec,
							'--refernce', self.input_file,
							'RS_Resequencing',
							self.output_directory,
							self.read_data + "/*.bax.h5"
							])							
		utils.syscall(command)
	