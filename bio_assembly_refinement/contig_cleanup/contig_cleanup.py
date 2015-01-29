import os
from fastaq import tasks
from pymummer import coords_file, alignment, nucmer

class ContigCleanup:
	def __init__(self, fasta_file, output_file="tmp.filtered.fa", cutoff_contig_length=10000, percent_match=95, debug=False):
		''' Constructor takes
		fasta file
		an output fasta filename
		a cutoff length to determine which contigs are too small to keep (default 10,000)
		a percent match to identify contained contigs (default 95)
		debug (do not delete temp files if set to true)
		'''
		self.fasta_file = fasta_file
		self.output_file = output_file
		self.cutoff_contig_length = cutoff_contig_length
		self.percent_match = percent_match
		self.debug = debug
		self.contigs = {}
		tasks.file_to_dict(self.fasta_file, self.contigs) #Read ids and sequences into hash

	
	def _find_small_contigs(self):
		'''	
		Remove contigs smaller than cutoff length 
		This is useful in bacterial assemblies to remove bits of sequence resulting from phages and plasmids
		'''
		small_contigs = []
		for id in self.contigs.keys():
			if len(self.contigs[id]) < self.cutoff_contig_length:
				small_contigs.append(id)
		return small_contigs		
				
				
	def _find_hits(self):
		'''Runs nucmer and returns a list of alignment objects (self hits removed)'''
		results_file = str(self.fasta_file) + ".coords"
		runner = nucmer.Runner(self.fasta_file, self.fasta_file, results_file, coords_header=False, maxmatch=True) # nucmer default breaklength is 200
		runner.run()
		file_reader = coords_file.reader(results_file)
		alignments = [coord for coord in file_reader if not coord.is_self_hit()] #Remove self hits (?) Keep these actually...
		return alignments

			
	def _find_contained_contigs(self):
		'''Parse alignments and identify contained contigs'''
		contained_contigs = []
		alignments = self._find_hits()		
		sorted_contig_ids = sorted(self.contigs.keys()) #Sorting so that results are consistent each time
		for contig_id in sorted_contig_ids:
			for algn in alignments:
				if algn.qry_name == contig_id \
				  and algn.ref_name != algn.qry_name \
				  and not algn.ref_name in contained_contigs \
				  and (algn.hit_length_qry/algn.qry_length) * 100 >= self.percent_match: #Does it have a very large hit?
					contained_contigs.append(contig_id)
		return contained_contigs
	
	
	def run(self):
		'''
		Remove small & contained contigs from fasta file to produce new fasta file
		Delete temp files unless debug is true
		'''
		contigs_to_remove = []
		contigs_to_remove.extend(self._find_small_contigs()) #TODO: Remove small contigs first, before running nucmer etc to save holding unnecessary contigs in memory
		contigs_to_remove.extend(self._find_contained_contigs())
		contig_ids_file = self.fasta_file + '.contig.ids.remove'
		with open(contig_ids_file, mode='wt') as ids_file:
			ids_file.write('\n'.join(contigs_to_remove))
		tasks.filter(self.fasta_file, self.output_file, ids_file=contig_ids_file, invert=True) #Invert has to be true here inorder for filter to remove contigs
		if not self.debug:
			os.remove(contig_ids_file)
