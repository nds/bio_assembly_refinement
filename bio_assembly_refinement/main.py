'''
Main module to run post assembly refinement on a fasta file
'''

import os
import tempfile
import shutil
from bio_assembly_refinement import contig_cleanup, circularisation, reassembly

class Main:
	def __init__(self, 
				dnaA_sequence, # Path to multi-fasta file with dnaA/refA/refB
				pacbio_exec, # Path to pacbio resequencing tool
				nucmer_exec, # Path to nucmer
    			cutoff_contig_length = 10000, # Contigs smaller than this will be disregarded
    			contained_percent_match = 95, # Percentage identity to determine if a contig is contained in another
				overlap_start_offset = 12, # A match has to start between 0 and this percentage of contig length when looking for overlapping ends
				overlap_max_length = 50, # A match has to end before this percentage of the contig length when looking for overlapping ends (default 50%)
				dnaA_match_percent_identity = 99, # A match to the dnaA/refA/refB has to have this percentage identity inorder to make that site the start of chromosome/plasmid
				debug=False):
				
		''' Constructor '''
		self.cutoff_contig_length = cutoff_contig_length
		self.contained_percent_match = contained_percent_match
		self.overlap_start_offset = overlap_start_offset
		self.overlap_max_length = overlap_max_length
		self.dnaA_match_percent_identity = dnaA_match_percent
		self.dnaA_sequence = dnaA_sequence
		self.pacbio_exec = pacbio_exec
		self.nucmer_exec = nucmer_exec
		self.debug = debug


	def process_assembly(self):
		'''Run three steps: clean contigs, trim & circularise, run pacbio resequencing'''	
		
		tmpdir = tempfile.mkdtemp(prefix='tmp.post_assembly_processing.', dir=os.getcwd())
		original_dir = os.getcwd()
		os.chdir(tmpdir)
		ccleaner = contig_cleanup.ContigCleanup(input_file, cutoff_contig_length=10)
		circulariser = circularisation.Circularisation(dnaA_sequence = test_dnaA_file,
													   contigs=test_contigs,
												       alignments = test_alignments,
												       dnaA_alignments = test_dnaA_alignments,
												       output_file="circularised.fa", 
												       offset = 7
												      )

		reassembler = reassembly.Reassembly(input_file=test_file,
											read_data=data_dir,
											pacbio_exec=data_dir + "/dummy_pacbio_script",
											output_directory="reassembled"
											)
		os.chdir(original_dir)
		shutil.rmtree(tmpdir)

   		 

