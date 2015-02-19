'''
Class to run post assembly refinement on a fasta file (filter, circularise, re-assemble)

Attributes:
-----------
fasta_file : input fasta file
dnaA_sequence : path to fasta file with dnaA/refA/refB 
bax_files : directory containing bax.h5 files
cutoff_contig_length : minimum contig length considered
contained_percent_match : percent identity to determine if contig is contained in another
overlap_offset : offset from the ends where an overlap can begin
overlap_max_length : maximum length of the overlap between ends
overlap_percent_identity : percent identity to use when determining if ends overlap
dnaA_hit_percent_identity : percent identity to use when looking at hits to dnaA
dnaA_hit_length_minimum : minimum length of hit to dnaA
working_directory : working directory (default current working directory) 
pacbio_exec : pacbio resequencing exec (default pacbio_smrtanalysis) 
nucmer_exec : nucmer exec (default nucmer) 
debug : do not delete temp files if set to true (default false)

Sample usage:
-------------

from bio_assembly_refinement import contig_cleanup

ccleaner = contig_cleanup.ContigCleanup("myassembly.fa")
ccleaner.run()
contigs_removed = ccleaner.get_filtered_contigs()

'''

import os
import tempfile
import shutil
from bio_assembly_refinement import contig_cleanup, circularisation, reassembly, utils

class Main:
	def __init__(self,
				fasta_file, 
				dnaA_sequence, 
				bax_files,
				cutoff_contig_length=10000,
				contained_percent_match=95,
				overlap_offset=49, 
				overlap_boundary_max=50, 
				overlap_min_length=2000,
				overlap_percent_identity=85,
				dnaA_hit_percent_identity=80,
				dnaA_hit_length_minimum=95,			
				working_directory=None, 
				pacbio_exec = "pacbio_smrtanalysis", 
				nucmer_exec = "nucmer", 
				reassembly_dir = "reassembly",
				debug = False
				):
		self.fasta_file = fasta_file 
		self.dnaA_sequence = dnaA_sequence
		self.bax_files = bax_files
		self.cutoff_contig_length = cutoff_contig_length
		self.contained_percent_match = contained_percent_match
		self.overlap_offset = overlap_offset 
		self.overlap_boundary_max = overlap_boundary_max
		self.overlap_min_length = overlap_min_length
		self.overlap_percent_identity = overlap_percent_identity
		self.dnaA_hit_percent_identity = dnaA_hit_percent_identity
		self.dnaA_hit_length_minimum = dnaA_hit_length_minimum		 
		self.pacbio_exec = pacbio_exec
		self.nucmer_exec = nucmer_exec 
		self.reassembly_dir = reassembly_dir
		self.debug = debug   		

		if not working_directory:
			self.working_directory = os.getcwd()		
		


	def process_assembly(self):
		'''Run three steps: clean contigs, trim & circularise, run pacbio resequencing'''	
		
		original_dir = os.getcwd()
		os.chdir(self.working_directory)
		
		ccleaner = contig_cleanup.ContigCleanup(fasta_file = self.fasta_file,
												working_directory = self.working_directory, 
												cutoff_contig_length=self.cutoff_contig_length,
												percent_match = self.contained_percent_match,
												debug = self.debug)
		ccleaner.run()
		
		circulariser = circularisation.Circularisation(fasta_file = ccleaner.get_results_file(), # Need the filename to retain naming scheme even though we pass in pre-computed contigs
													   dnaA_sequence = self.dnaA_sequence,
													   working_directory = self.working_directory,
													   contigs = ccleaner.get_filtered_contigs(),
												       alignments = ccleaner.get_alignments(),
												       overlap_offset = self.overlap_offset,
												       overlap_boundary_max = self.overlap_boundary_max,
												       overlap_min_length = self.overlap_min_length,
												       overlap_percent_identity = self.overlap_percent_identity,
												       dnaA_hit_percent_identity = self.dnaA_hit_percent_identity,
												       dnaA_hit_length_minimum = self.dnaA_hit_length_minimum,
													   debug = self.debug
												      )
												      
		circulariser.run()      
				
		 
		reassembler = reassembly.Reassembly(input_file=circulariser.get_results_file(),
											read_data=self.bax_files,
											pacbio_exec=self.pacbio_exec,
											working_directory = self.working_directory,
											output_directory = self.reassembly_dir,
											debug = self.debug
											)
											
		reassembler.run()
		
		if not self.debug:
			utils.delete(ccleaner.get_results_file())
# 			utils.delete(circulariser.get_results_file()) #Only delete once code added to wait for bsub to finish
		
		os.chdir(original_dir)
   		 

