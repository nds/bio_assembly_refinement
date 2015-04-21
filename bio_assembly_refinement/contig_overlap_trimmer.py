'''
Class to find and trim overlapping ends of a contig

Attributes:
-----------
fasta_file : input fasta file
working_directory : path to working directory (default to current working directory)
contigs : dict of contigs (instead of fasta file)
trim : trim overlaps (default true)
trim_reversed_overlaps: trims overlaps even if reversed (default false)
alignments : pre-computed alignments (if available from previous step)
overlap_offset: offset from edge that the overlap can start (default 1000)
overlap_boundary_max : max boundary of overlap expressed as % of length of reference (default 50)
overlap_min_length : minimum length of overlap (default 1KB)
overlap_max_length : maximum length of overlap (default 3KB)
overlap_percent_identity : percent identity of match between ends (default 85)
min_trim_length : minimum trimmed length of contig over total contig length (default 0.8)
summary_file :  summary file (default contig_overlap_summary.txt)
debug : do not delete temp files if set to true (default false)
			  
Sample usage:
-------------


'''

import os
import re
from pyfastaq import tasks, sequences
from pyfastaq import utils as fastaqutils
from pymummer import alignment
from bio_assembly_refinement import utils, contig_history

class ContigOverlapTrimmer:
	def __init__(self, 
				 fasta_file='', 
				 working_directory=None, 
				 contigs={},
				 alignments=[],
				 trim = True,
				 trim_reversed_overlaps = False,
				 overlap_offset=1000, 
				 overlap_boundary_max=50, 
				 overlap_min_length=1000,
				 overlap_max_length=3000,
				 overlap_percent_identity=85,
				 min_trim_length=0.89,
				 summary_file = "circularisation_summary_file.txt",			  
				 debug=False):

		''' Constructor '''
		self.fasta_file = fasta_file
		self.working_directory = working_directory if working_directory else os.getcwd()		
		self.contigs = contigs
		self.alignments = alignments
		self.trim = trim
		self.trim_reversed_overlaps = trim_reversed_overlaps
		self.overlap_offset = overlap_offset
		self.overlap_boundary_max = overlap_boundary_max * 0.01
		self.overlap_min_length = overlap_min_length
		self.overlap_max_length = overlap_max_length
		self.overlap_percent_identity = overlap_percent_identity
		self.min_trim_length = min_trim_length
		self.summary_file = summary_file
		self.output_file = self._build_final_filename()
		self.debug = debug
		
		# Extract contigs and generate nucmer hits if not provided
		if not self.contigs:
			self.contigs = {}
			tasks.file_to_dict(self.fasta_file, self.contigs) 
		
		# Run nucmer
		if not self.alignments:
			self.alignments = utils.run_nucmer(self.fasta_file, self.fasta_file, self._build_alignments_filename(), min_percent_id=self.overlap_percent_identity)
		
		
	def _find_best_overlap(self, contig_id):
		''' Look for the (best) overlap'''		
		best_overlap = None
		boundary = self.overlap_boundary_max * len(self.contigs[contig_id])
		for algn in self.alignments:
			if algn.qry_name == contig_id and \
			   algn.ref_name == contig_id and \
			   algn.ref_start < self.overlap_offset and \
			   algn.ref_end < boundary and \
			   algn.qry_start > boundary and \
			   algn.qry_end > (algn.qry_length - self.overlap_offset) and \
			   algn.hit_length_ref >= self.overlap_min_length and \
			   algn.hit_length_ref <= self.overlap_max_length and \
			   algn.percent_identity > self.overlap_percent_identity:
				if not best_overlap or \
				   (algn.ref_start <= best_overlap.ref_start and \
					algn.qry_end > best_overlap.qry_end ):
				   best_overlap = algn
		return best_overlap
		
		
	def _trim(self, contig_id, best_overlap):
		''' trim overlap off the start of contig '''
		original_sequence = self.contigs[contig_id]
		trim_start = best_overlap.ref_end+1
		trim_end = best_overlap.qry_end+1
		trim_status = ''
		if not best_overlap.on_same_strand:
			if not self.trim_reversed_overlaps:
				trim_status = "overlap reversed, not trimming"
				return trim_status
			else:
				trim_start = min(best_overlap.ref_start, best_overlap.ref_end) + 1
				trim_end = max(best_overlap.qry_start, best_overlap.qry_end) + 1
				trim_status = "overlap reversed, trimming"
		trimmed_sequence = original_sequence[trim_start:trim_end]		
		if(len(trimmed_sequence)/len(original_sequence) < self.min_trim_length):
			trim_status = "trimmed length would be too short, not trimming"
			return trim_status
		else:
			self.contigs[contig_id] = trimmed_sequence
			trim_status = "trimmed length " + str(len(trimmed_sequence))
		return trim_status
		
		
	def _write_summary(self, contig_id, best_overlap, trim_status):
		'''Write summary'''
		if (not os.path.exists(self.summary_file)) or os.stat(self.summary_file).st_size == 0:
			header = '\t'.join(['id', 'overlap length', 'overlap location', 'trim status']) +'\n'
			utils.write_text_to_file(header, self.summary_file)
		overlap_length = '-'
		overlap_location = '-'
		if best_overlap:
			overlap_length = str(best_overlap.hit_length_ref)
			overlap_location = str(best_overlap.ref_start) + ',' + str(best_overlap.ref_end) + '-' + \
							   str(best_overlap.qry_start) + ',' + str(best_overlap.qry_start)
		line = "\t".join([contig_id, overlap_length, overlap_location, trim_status]) + "\n"
		utils.write_text_to_file(line, self.summary_file)
						   
						   
	def _build_alignments_filename(self):
		return os.path.join(self.working_directory, "nucmer_all_contigs.coords")
		
		
	def _build_final_filename(self):
		input_filename = os.path.basename(self.fasta_file)
		return os.path.join(self.working_directory, "trimmed_" + input_filename)	
		
				
	def _build_intermediate_filename(self):
		input_filename = os.path.basename(self.fasta_file)
		return os.path.join(self.working_directory, "unsorted_trimmed_" + input_filename)	
			   
			   
	def run(self):	
		original_dir = os.getcwd()
		os.chdir(self.working_directory)	
		output_fw = fastaqutils.open_file_write(self._build_intermediate_filename())
		for contig_id in self.contigs.keys():
			#Look for overlaps, trim if applicable
			best_overlap = self._find_best_overlap(contig_id)
			trim_status = ''			
			if best_overlap and self.trim:
				trim_status = self._trim(contig_id, best_overlap)
			self._write_summary(contig_id, best_overlap, trim_status)
			print(sequences.Fasta(contig_id, self.contigs[contig_id]), file=output_fw)	
					
		tasks.sort_by_size(self._build_intermediate_filename(), self.output_file) # Sort contigs in final file according to size
		
		if not self.debug:
			utils.delete(self._build_alignments_filename())
			utils.delete(self._build_intermediate_filename())

		os.chdir(original_dir)
		
	