''' 
Class to remove small and contained contigs 

Attributes:
-----------
fasta_file : input fasta file name
working_directory : path to working directory (default to current working directory)
cutoff_contig_length : contigs smaller than this will be disregarded (default 10,000)
percent_match : percent identity of nucmer hit when deciding if contig is contained in another
summary_file : summary file
debug : do not delete temp files if set to true (default false)

Sample usage:
-------------

from bio_assembly_refinement import contig_cleanup

ccleaner = contig_cleanup.ContigCleanup("myassembly.fa")
ccleaner.run()
contigs_removed = ccleaner.get_filtered_contigs()

'''

import os
from pyfastaq import tasks
from pymummer import alignment
from bio_assembly_refinement import utils

class ContigCleanup:
	def __init__(self, 
				 fasta_file, 
				 working_directory=None, 
				 cutoff_contig_length=2000, 
				 percent_match=95, 
				 summary_file="contig_filtration_summary.txt",
				 debug=False):
		''' Constructor '''
		self.fasta_file = fasta_file
		self.working_directory = working_directory if working_directory else os.getcwd()			
		self.cutoff_contig_length = cutoff_contig_length
		self.percent_match = percent_match
		self.summary_file = summary_file
		self.debug = debug		
		self.contigs = {}
		tasks.file_to_dict(self.fasta_file, self.contigs) #Read contig ids and sequences into dict
		self.alignments = []
		self.filtered_contigs = []
		self.output_file = self._build_final_filename()

	
	def _find_small_contigs(self):
		'''	Remove contigs smaller than cutoff length. '''
		small_contigs = []
		for id in self.contigs.keys():
			if len(self.contigs[id]) < self.cutoff_contig_length:
				small_contigs.append(id)
		return small_contigs		
	
	
	def _find_contained_contigs(self, filename):
		'''Parse alignments and identify contained contigs'''
		self.alignments = utils.run_nucmer(filename, filename, self._build_alignments_filename())
		contained_contigs = []
		sorted_contig_ids = sorted(self.contigs.keys()) #Sorting so that results are consistent each time
		for contig_id in sorted_contig_ids:
			for algn in self.alignments:
				if (not algn.is_self_hit()) \
				   and algn.qry_name == contig_id \
				   and algn.ref_name != algn.qry_name \
				   and not algn.ref_name in contained_contigs \
				   and (algn.hit_length_qry/algn.qry_length) * 100 >= self.percent_match: #Does it have a very large hit to another contig?
					contained_contigs.append(contig_id)
		return contained_contigs
		
		
	def _produce_summary(self, small_contigs, contained_contigs):
		'''Generate summary text and write to summary file'''
		text = 	 "~~Contig filtration~~ \n" + \
					"Small contigs removed: \n" + \
					"\n".join(small_contigs) + "\n" \
					"Contained contigs removed: \n" + \
					"\n".join(contained_contigs) + "\n"
		utils.write_text_to_file(text, self.summary_file)
								  		
		
	def _build_intermediate_filename(self):
		return os.path.join(self.working_directory, "intermediate.fa")
		
	
	def _build_alignments_filename(self):
		return os.path.join(self.working_directory, "nucmer_all_contigs.coords")
		
		
	def _build_final_filename(self):
		input_filename = os.path.basename(self.fasta_file)
		return os.path.join(self.working_directory, "filtered_" + input_filename)	
	
	
	def run(self):
		'''Produce a filtered fasta file.'''	
		original_dir = os.getcwd()
		os.chdir(self.working_directory)
		intermediate_file = self._build_intermediate_filename()
				
		small_contigs = self._find_small_contigs()
		contig_ids_file_1 = utils.write_ids_to_file(small_contigs, "contig.ids.too.small")  
		tasks.filter(self.fasta_file, intermediate_file, ids_file=contig_ids_file_1, invert=True)
			
		contained_contigs = self._find_contained_contigs(intermediate_file) #Run nucmer after filtering small contigs to save holding unnecessary contigs and hits in memory
		contig_ids_file_2 = utils.write_ids_to_file(contained_contigs, "contig.ids.contained")
		tasks.filter(intermediate_file, self.output_file, ids_file=contig_ids_file_2, invert=True)
		
		for id in small_contigs + contained_contigs:
			del self.contigs[id] #No longer care about contigs thrown away			
	
		self._produce_summary(small_contigs, contained_contigs)
		
		if not self.debug:
			utils.delete(contig_ids_file_1)
			utils.delete(contig_ids_file_2)
			utils.delete(intermediate_file)
			utils.delete(self._build_alignments_filename())
		
		os.chdir(original_dir)