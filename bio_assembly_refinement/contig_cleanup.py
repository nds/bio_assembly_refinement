''' 
Class to remove small and contained contigs 

Attributes:
-----------
fasta_file : input fasta file name
working_directory : path to working directory (default to current working directory)
cutoff_contig_length : contigs smaller than this will be disregarded (default 200)
percent_match : percent identity of nucmer hit when deciding if contig is contained in another (default 95)
skip : contig ids to skip i.e. keep no matter what (file or list)
summary_file : summary file
summary_prefix : prefix for lines in summary file
debug : do not delete temp files if set to true (default false)

Sample usage:
-------------

from bio_assembly_refinement import contig_cleanup

ccleaner = contig_cleanup.ContigCleanup("myassembly.fa")
ccleaner.run()
ccleaner.output_file will be the cleaned fasta file
ccleaner.summary_file will be the summary file

'''

import os
from pyfastaq import tasks
from pymummer import alignment
from bio_assembly_refinement import utils
from pyfastaq import utils as fastaqutils
from pyfastaq import sequences

class ContigCleanup:
	def __init__(self, 
				 fasta_file, 
				 working_directory=None, 
				 cutoff_contig_length=2000, 
				 percent_match=95, 
				 skip = None,
				 summary_file="contig_cleanup_summary.txt",
				 summary_prefix="[contig cleanup]",
				 debug=False):
				 
		''' Constructor '''
		self.fasta_file = fasta_file
		self.working_directory = working_directory if working_directory else os.getcwd()			
		self.cutoff_contig_length = cutoff_contig_length
		self.percent_match = percent_match
		self.summary_file = summary_file
		self.summary_prefix = summary_prefix
		self.debug = debug		
		self.contigs = {}
		tasks.file_to_dict(self.fasta_file, self.contigs) #Read contig ids and sequences into dict
		
		self.ids_to_skip = set()		
		if skip:
			if type(skip) == set:
				self.ids_to_skip = set(skip) # Assumes ids is a list
			else:
				fh = fastaqutils.open_file_read(skip)
				for line in fh:
					self.ids_to_skip.add(line.rstrip())
				fastaqutils.close(fh)
		self.output_file = self._build_final_filename()		
	
	
	def _write_summary(self, small_contigs, contained_contigs):
		'''Write summary'''
		header = self.summary_prefix + " contig_id\tclean_log\n"
		utils.write_text_to_file(header, self.summary_file)
		for contig_id in self.contigs.keys():
			if contig_id in small_contigs:
				line = self.summary_prefix + " " + contig_id + " small (removed)\n"
			elif contig_id in contained_contigs:
				line = self.summary_prefix + " " + contig_id + " contained (removed)\n"
			elif contig_id in self.ids_to_skip:
				line = self.summary_prefix + " " + contig_id + " skipped\n"
			else:
				line = self.summary_prefix + " " + contig_id + " kept\n"
			utils.write_text_to_file(line, self.summary_file)
				
		
	def _build_final_filename(self):
		'''Build output filename'''
		input_filename = os.path.basename(self.fasta_file)
		return os.path.join(self.working_directory, "filtered_" + input_filename)	
	
	
	def _build_nucmer_filename(self):
		'''Build temp nucmer filename'''
		return os.path.join(self.working_directory, "nucmer_all_contigs.coords")
		
	
	def run(self):
		'''Produce a filtered fasta file.'''	
		original_dir = os.getcwd()
		os.chdir(self.working_directory)
		small_contigs = set()
		contained_contigs = set()
		if len(self.contigs) > len(self.ids_to_skip):
			alignments = utils.run_nucmer(self.fasta_file, self.fasta_file, self._build_nucmer_filename(), min_percent_id=self.percent_match, run_promer=False)
			for id in self.contigs.keys():
				if not id in self.ids_to_skip:
					if len(self.contigs[id]) < self.cutoff_contig_length:
						small_contigs.add(id)
					else:
						for algn in alignments:
							if (not algn.is_self_hit()) \
							   and algn.qry_name == id \
							   and algn.ref_name != algn.qry_name \
							   and not algn.ref_name in contained_contigs \
							   and (algn.hit_length_qry/algn.qry_length) * 100 >= self.percent_match:
								contained_contigs.add(id)
					
			discard = small_contigs.union(contained_contigs)
			ids_file = utils.write_ids_to_file(discard, "contig.ids.discard")  
			tasks.filter(self.fasta_file, self.output_file, ids_file=ids_file, invert=True)	
								
			if not self.debug:
				utils.delete(ids_file)
				utils.delete(self._build_nucmer_filename())
		else:
			output_fw = fastaqutils.open_file_write(self.output_file)
			for contig_id in self.contigs:
				print(sequences.Fasta(contig_id, self.contigs[contig_id]), file=output_fw)
			fastaqutils.close(output_fw)
		
		self._write_summary(small_contigs, contained_contigs)	
		os.chdir(original_dir)
