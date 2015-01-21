# Run nucmer on given fasta file (against itself), analyse results and remove any overlapping
# contigs
# The rules to define what 'overlapping' means are encoded below
# In the script, leave an option to retain all intermediate steps
# If not, delete the temp directory
# TODO: Catch errors

import os
from fastaq import tasks
from pymummer import coords_file, alignment, nucmer

class ContigCleanup:
	def __init__(self, fastafile, outfile="temp.filtered.fa", cutoff_contig_length=10000):
		''' 
		Constructor
		Reads in contig data from fasta file
		'''
		self.fastafile = fastafile
		self.outfile = outfile
		self.cutoff_contig_length = cutoff_contig_length
		self.contigs = {}
		tasks.file_to_dict(fastafile, self.contigs) #Read ids and sequences into hash

	
	def _find_small_contigs(self):
		'''
		Remove contigs smaller than set cutoff length (default 10000 ?)
		This is useful in bacterial assemblies to remove bits of sequence resulting from phages and plasmids
		'''
		small_contigs = []
		for id in self.contigs.keys():
			if len(self.contigs[id]) < self.cutoff_contig_length:
				small_contigs.append(id)
		return small_contigs		
				
				
	def _find_hits(self):
		'''Runs nucmer and returns a list of alignment objects (self hits removed)'''
		resultsfile = str(self.fastafile) + ".coords"
		runner = nucmer.Runner(self.fastafile, self.fastafile, resultsfile, coords_header=False, maxmatch=True) # nucmer default breaklength is 200
		runner.run()
		filereader = coords_file.reader(resultsfile)
		alignments = [coord for coord in filereader if not coord.is_self_hit()] #Remove self hits
		return alignments

			
	def _find_contained_contigs(self):
		'''
		Parse alignments and identify contained contigs
		Returns list of contig ids
		'''
		contained_contigs = ['TEST_CONTIG_11', 'TEST_CONTIG_12']
# 		alignments = self._find_hits()
		# For each contig, see if it has a reasonably large hit			
		return contained_contigs
	
	
	def run(self, debug=False):
		''' Get list of small contigs
			Get list of contained contigs
			Remove above contigs from fasta file to produce new fasta file
			Delete temp files unless debug is true
		'''
		contigs_to_remove = []
		# When filter can take a set of contigs, probably best to remove small contigs first before running nucmer etc to save holding any small contigs in memory 
		contigs_to_remove.extend(self._find_small_contigs())
		contigs_to_remove.extend(self._find_contained_contigs())
		contig_ids_file = 'tmp.contig.ids.file'
		with open(contig_ids_file, mode='wt', encoding='utf-8') as ids_file:
			ids_file.write('\n'.join(contigs_to_remove))
		tasks.filter(self.fastafile, self.outfile, ids_file=contig_ids_file, invert=True) #Invert has to be true here inorder for filter to remove contigs in list
		#Delete ids file
