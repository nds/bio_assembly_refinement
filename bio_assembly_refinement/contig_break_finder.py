''' 
Find the point at which to break a circular contig (at origin of replication or random gene)

Attributes:
-----------
fasta_file : input fasta file name
gene_file : file with genes (e.g. dnaA)
skip : list or file of contig ids to skip (i.e. do not break at gene location)
hit_percent_id : min percent identity of matches to gene
match_length_percent : min length of ref match expressed as % of gene length (default 80%)
choose_random_gene : if genes in file cannot be found, run prodigal and find random gene (default True)
rename : rename contigs (default True)
working_directory : path to working directory (default to current working directory)
summary_file : summary file
debug : do not delete temp files if set to true (default false)

Sample usage:
-------------

from bio_assembly_refinement import contig_break_finder

break_finder = contig_break_finder.ContigBreakFinder(fasta_file = myfasta_file,
													 gene_file = mydnaA_file,													     								       
												     )
break_finder.run()
break_finder.output_file will be the cleaned fasta file
break_finder.summary_file will be the summary file
'''

import os
import shutil
from bio_assembly_refinement import utils
from pyfastaq import sequences
from pyfastaq import utils as fastaqutils
from pyfastaq import tasks
from pymummer import alignment


class ContigBreakFinder:
	def __init__(self, 
			     fasta_file, 
			     gene_file, 
			     skip=None, #Avoid circularising contigs with these ids
			     hit_percent_id=80, 
			     match_length_percent=100, 
			     choose_random_gene=True, 
			     rename=True,
			     working_directory=None,
			     summary_file = "contig_breaks_summary.txt",			  
				 debug=False):				 
		''' Attributes '''
		self.fasta_file = fasta_file
		self.gene_file = gene_file
		self.hit_percent_id = hit_percent_id
		self.match_length_percent = match_length_percent
		self.choose_random_gene = choose_random_gene
		self.rename = rename
		self.working_directory = working_directory if working_directory else os.getcwd()
		self.summary_file = summary_file
		self.output_file = self._build_final_filename()
		self.debug = debug
		self.contigs = {}
		tasks.file_to_dict(self.fasta_file, self.contigs) #Read contig ids and sequences into dict
		# run promer
		self.dnaA_alignments = utils.run_nucmer(self.fasta_file, self.gene_file, self._build_promer_filename(), min_percent_id=self.hit_percent_id, run_promer=True)
		self.random_gene_starts = {}
		if self.choose_random_gene:
			self.random_gene_starts = self._run_prodigal_and_get_gene_starts()
		
		self.ids_to_skip = set()		
		if skip:
			if type(skip) == set:
				self.ids_to_skip = set(skip) # Assumes ids is a list
			else:
				fh = fastaqutils.open_file_read(skip)
				for line in fh:
					self.ids_to_skip.add(line.rstrip())
				fastaqutils.close(fh)
	
		
	def _run_prodigal_and_get_gene_starts(self):
		'''Run prodigal and find gene starts''' 
		gene_starts = {}
		fastaqutils.syscall("prodigal -i " + self.fasta_file + " -o " + self._build_prodigal_filename() +  " -f gff -c -q")	# run on whole fasta as prodgal works better with larger sequences
		fh = fastaqutils.open_file_read(self._build_prodigal_filename())
		for line in fh:
			if not line.startswith("#"):
				columns = line.split('\t')
				start_location = int(columns[3])
				contig_id = columns[0]
				boundary_start = round(0.4 * len(self.contigs[contig_id]))
				boundary_end = round(0.6 * len(self.contigs[contig_id]))
				if start_location > boundary_start and start_location < boundary_end:
					gene_starts[contig_id] = start_location - 1 #Interbase				
		fastaqutils.close(fh)
		return gene_starts	
		
		
	def _build_final_filename(self):
		'''Build output filename'''
		input_filename = os.path.basename(self.fasta_file)
		return os.path.join(self.working_directory, "circularised_" + input_filename)
		
		
	def _build_promer_filename(self):
		'''Build temp promer filename'''
		return os.path.join(self.working_directory, "promer_dnaA_hits.coords")
		
		
	def _build_prodigal_filename(self):
		'''Build temp prodigal filename'''
		return os.path.join(self.working_directory, "prodigal_genes.gff")
	
		
	def _write_summary(self, contig_id, break_point, gene_name, gene_reversed, new_name):
		'''Write summary'''
		if (not os.path.exists(self.summary_file)) or os.stat(self.summary_file).st_size == 0:
			header = '\t'.join(['id', 'break_point', 'gene_name', 'gene_reversed', 'new_name']) +'\n'
			utils.write_text_to_file(header, self.summary_file)
		name_to_print = new_name if self.rename else '-'			
		text = '\t'.join(map(str, [contig_id, break_point, gene_name, gene_reversed, name_to_print])) + '\n'
		utils.write_text_to_file(text, self.summary_file)		
	

	def run(self):
		'''Look for break point in contigs and rename if needed'''	
		
		chromosome_count = 1
		plasmid_count = 1
		output_fw = fastaqutils.open_file_write(self.output_file)
		for contig_id in self.contigs:
			contig_sequence = self.contigs[contig_id]
			if contig_id not in self.ids_to_skip:		
#				print("Working on contig: " + contig_id)
				dnaA_found = False
				gene_name = '-'
				gene_on_reverse_strand = False
				new_name = contig_id #Stick with old name if no new name comes along
				break_point = 0
				
				for algn in self.dnaA_alignments:			
					if algn.ref_name == contig_id and \
					   algn.hit_length_qry >= (algn.qry_length * self.match_length_percent/100) and \
					   algn.percent_identity >= self.hit_percent_id and \
					   algn.qry_start == 0:	     
#						print("dnaA found")
						dnaA_found = True
						gene_name = algn.qry_name
						if algn.on_same_strand():
							break_point = algn.ref_start						
						else:
							# Reverse complement sequence, circularise using new start of dnaA in the right orientation
#							print("dnaA reversed")
#							sequence_tmp = str(original_sequence)
#							original_sequence = sequence_tmp.translate(str.maketrans("ATCGatcg","TAGCtagc"))[::-1]
							contig_sequence.revcomp()
							break_point = (algn.ref_length - algn.ref_start) - 1 #interbase
							gene_on_reverse_strand = True
						new_name = 'chromosome_' + str(chromosome_count)
						chromosome_count += 1		
						break;
					
				if not dnaA_found:
					if len(self.contigs[contig_id]) < 200000: # Only rename if it's roughly plasmid size
						new_name = 'plasmid_' + str(plasmid_count)
						plasmid_count += 1
					if self.choose_random_gene and self.random_gene_starts[contig_id]:			
						break_point = self.random_gene_starts[contig_id]
						gene_name = 'prodigal'
				
				if break_point > 0:
#					print("Reorganising sequence")
					contig_sequence = contig_sequence[break_point:] + contig_sequence[0:break_point]
#					print(original_sequence)
				
				contig_name = new_name if self.rename else contig_id
				print(sequences.Fasta(contig_name, contig_sequence), file=output_fw)
				self._write_summary(contig_id, break_point, gene_name, gene_on_reverse_strand, contig_name)
				
		output_fw.close()
		
		if not self.debug:
			utils.delete(self._build_promer_filename())
			utils.delete(self._build_prodigal_filename())
