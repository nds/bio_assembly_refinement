'''
Class to run post assembly refinement on a fasta file (filter, circularise, re-assemble)

Attributes:
-----------
fasta_file : input fasta file
dnaA_sequence : path to fasta file with dnaA/refA/refB 
bax_files : directory containing bax.h5 files
cutoff_contig_length : minimum contig length considered
contained_percent_match : percent identity to determine if contig is contained in another
overlap_offset : offset from the ends where an overlap can begin (1000)
overlap_boundary_max : max boundary of overlap expressed as % of length of reference (default 50)
overlap_min_length : minimum length of overlap (default 1KB)
overlap_max_length : minimum length of overlap (default 3KB)
overlap_percent_identity : percent identity to use when determining if ends overlap
min_trim_length : minimum trimmed length of contig over total contig length (default 0.8)
trim : trim overlaps (default true)
trim_reversed_overlaps: trims overlaps even if reversed (default false)
dnaA_hit_percent_identity : percent identity to use when looking at hits to dnaA
dnaA_hit_length_minimum : minimum length of hit to dnaA
no_bsub : If set, will run quiver as a child process and not bsub it (useful for pipeline)
working_directory : working directory (default current working directory) 
pacbio_exec : pacbio resequencing exec (default pacbio_smrtanalysis) 
nucmer_exec : nucmer exec (default nucmer) 
reassembly_dir : directory sent to quiver (default reassembly)
summary_file :  summary file (default pacbio_postprocess_summary.txt)
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
from bio_assembly_refinement import contig_cleanup, contig_break_finder, contig_overlap_trimmer, reassembly, utils

class Main:
	def __init__(self,
				fasta_file, 
				dnaA_sequence, 
				bax_files,  
				# clean up arguments
				cutoff_contig_length=10000,
				contained_percent_match=95,
				# trimming arguments
				overlap_offset=1000, 
				overlap_boundary_max=50, 
				overlap_min_length=1000,
				overlap_max_length=3000,
				overlap_percent_identity=85,
				min_trim_length = 0.8,
				trim = True,
				trim_reversed_overlaps = False,
				# contig breaker arguments
				dnaA_hit_percent_identity=80,
				dnaA_hit_length_minimum=65,	
				# quiver arguments	
				no_bsub = False,	
				pacbio_exec = "pacbio_smrtanalysis", 
				# general arguments
				working_directory=None,
				reassembly_dir = "improved_assembly",
				summary_file = "assembly_refinement_summary.txt",
				debug = False
				):
		self.fasta_file = fasta_file 
		self.dnaA_sequence = dnaA_sequence
		self.bax_files = bax_files
		#clean up arguments
		self.cutoff_contig_length = cutoff_contig_length
		self.contained_percent_match = contained_percent_match
		# trimming arguments
		self.overlap_offset = overlap_offset 
		self.overlap_boundary_max = overlap_boundary_max
		self.overlap_min_length = overlap_min_length
		self.overlap_max_length = overlap_max_length
		self.overlap_percent_identity = overlap_percent_identity
		self.min_trim_length = min_trim_length
		self.trim = trim
		self.trim_reversed_overlaps = trim_reversed_overlaps
		# contig break arguments
		self.dnaA_hit_percent_identity = dnaA_hit_percent_identity
		self.dnaA_hit_length_minimum = dnaA_hit_length_minimum	
		# quiver arguments
		self.no_bsub = no_bsub
		self.pacbio_exec = pacbio_exec
		# general arguments
		self.working_directory = working_directory if working_directory else os.getcwd()	 
		self.reassembly_dir = reassembly_dir
		self.summary_file = summary_file
		self.debug = debug   		
		

	def process_assembly(self):
		'''Run three steps: clean contigs, trim & circularise, run pacbio resequencing'''	
		
		original_dir = os.getcwd()
		os.chdir(self.working_directory)
		
		# Step 1: Clean
		ccleaner = contig_cleanup.ContigCleanup(fasta_file = self.fasta_file,
												working_directory = self.working_directory, 
												cutoff_contig_length=self.cutoff_contig_length,
												percent_match = self.contained_percent_match,
												debug = self.debug)
		ccleaner.run()
		
		# Step 2: Trim
		contig_trimmer = contig_overlap_trimmer.ContigOverlapTrimmer(fasta_file = ccleaner.output_file, 
													   working_directory = self.working_directory,
												       trim = self.trim,
				 									   trim_reversed_overlaps = self.trim_reversed_overlaps,
												       overlap_offset = self.overlap_offset,
												       overlap_boundary_max = self.overlap_boundary_max,
												       overlap_min_length = self.overlap_min_length,
												       overlap_max_length = self.overlap_max_length,
												       overlap_percent_identity = self.overlap_percent_identity,
												       min_trim_length = self.min_trim_length,
													   debug = self.debug
												      )												      
		contig_trimmer.run()      
		
		# Step 3: Break contigs at same point		
		if os.path.exists(contig_trimmer.output_file):
			contig_breaker = contig_break_finder.ContigBreakFinder(fasta_file = contig_trimmer.output_file, 
																   gene_file = self.dnaA_sequence,
																   hit_percent_id = self.dnaA_hit_percent_identity,
																   match_length_percent = self.dnaA_hit_length_minimum,
																   working_directory = self.working_directory,
																   rename = False,	
																   debug = self.debug
																  )
			contig_breaker.run()	
			
		# Step 4: Run quiver
		if os.path.exists(contig_breaker.output_file):
			reassembler = reassembly.Reassembly(input_file=contig_breaker.output_file,
												read_data=self.bax_files,
												pacbio_exec=self.pacbio_exec,
												no_bsub = self.no_bsub,
												working_directory = self.working_directory,
												debug = self.debug,
												)
											
			reassembler.run()
			
		
		if not self.debug:
			utils.delete(ccleaner.output_file)
			utils.delete(contig_trimmer.output_file)
#			utils.delete(contig_breaker.output_file) # delete after bsub quiver has finished (to implement)

		
		os.chdir(original_dir)
   		 

