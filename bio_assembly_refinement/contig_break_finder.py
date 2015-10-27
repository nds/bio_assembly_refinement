''' 
Find the best point at which to break a circular contig (at origin of replication or random gene)
Uses promer and prodigal

Attributes:
-----------
fasta_file : input fasta file name
gene_file : file with genes (e.g. dnaA)
skip : list or file of contig ids to skip 
hit_percent_id : min percent identity of matches to gene
match_length_percent : min length of ref match expressed as % of gene length (default 80%)
choose_random_gene : if genes in file cannot be found, run prodigal and find random gene (default True)
rename : rename contigs to be either chromosomes or plasmids depending on length and presence/absence of dnaA (default True)
working_directory : path to working directory (default to current working directory)
summary_file : summary file
summary prefix : the prefix for each line in summary file
debug : do not delete temp files if set to true (default false)

Sample usage:
-------------

from bio_assembly_refinement import contig_break_finder

break_finder = contig_break_finder.ContigBreakFinder(fasta_file = myfasta_file,
													 gene_file = mydnaA_file,													     								       
												     )
break_finder.run()
break_finder.output_file will be the resulting fasta file
break_finder.summary_file will be the summary file

'''

import os
import shutil
from bio_assembly_refinement import utils, prodigal_hit
from pyfastaq import sequences, tasks, intervals
from pyfastaq import utils as fastaqutils
from pymummer import alignment
import subprocess



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
			     summary_prefix="[contig break finder]",			  
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
		self.summary_prefix = summary_prefix
		self.output_file = self._build_final_filename()
		self.debug = debug
		self.contigs = {}
		tasks.file_to_dict(self.fasta_file, self.contigs) #Read contig ids and sequences into dict
		self.random_gene_starts = {}

		self.ids_to_skip = set()		
		if skip:
			if type(skip) == set:
				self.ids_to_skip = set(skip) # Assumes ids is a list
			else:
				fh = fastaqutils.open_file_read(skip)
				for line in fh:
					self.ids_to_skip.add(line.rstrip())
				fastaqutils.close(fh)
				
	
	def _get_length_of_fasta_file(self):
		''' Sum up lengths of all contigs'''
		d = {k: len(v) for k, v in self.contigs.items()}
		return sum(d.values())
	
		
	def _run_prodigal_and_get_gene_starts(self):
		'''Run prodigal and find the start of a gene around the middle of each contig''' 
		gene_starts = {}
		# run prodigal
		prodigal_output = utils.run_prodigal(self.fasta_file, self._build_prodigal_filename(), self._get_length_of_fasta_file())
		prodigal_genes = {}
		if(prodigal_output):
			fh = fastaqutils.open_file_read(self._build_prodigal_filename())
			for line in fh:
				if not line.startswith("#"):
					columns = line.split('\t')
					start_location = int(columns[3])
					end_location = int(columns[4])
					contig_id = columns[0]
					strand = columns[6]	
					middle = abs((len(self.contigs[contig_id])/2))
					p = prodigal_hit.ProdigalHit(start_location, end_location, strand, middle)				
					prodigal_genes.setdefault(contig_id, []).append(p)
			fastaqutils.close(fh)
			# look for best distance
			for id in self.contigs.keys():
				best_gene = None
				if id in prodigal_genes.keys():
					all_prodigal_hits = prodigal_genes[id]
					min_distance = abs(len(self.contigs[contig_id])/2)
					for p in all_prodigal_hits:
						if p.distance <= min_distance:
							best_gene = p
							min_distance = p.distance
				if best_gene:
					gene_starts[id] = best_gene
				else:
					gene_starts[id] = None # Could not find a gene			
		return gene_starts	
		
		
	def _best_promer_hit(self, promer_hits, contig_id):
		'''Look for the best promer hit'''
		complete_dnaA_found = False
		break_point = None
		on_reverse_strand = False
		dnaA_name = None
		if not promer_hits:
			return complete_dnaA_found, break_point, on_reverse_strand, dnaA_name
		
		best_hit = promer_hits[0]		
		partial_hits_at_start = []
		partial_hits_at_end = []	
			
		for algn in promer_hits:
			# look for a complete match	
			if self._is_complete_match_to_dnaa(algn, best_hit, contig_id):
				complete_dnaA_found = True
				best_hit = algn
			else:
				# look for partial matches at the start and end
				if self._is_partial_match_at_start(algn, contig_id):
					partial_hits_at_start.append(algn)
				if self._is_partial_match_at_end(algn, contig_id):
					partial_hits_at_end.append(algn)
					
		if complete_dnaA_found:
			dnaA_name = best_hit.qry_name
			if best_hit.on_same_strand():
				break_point = best_hit.ref_start						
			else:
				break_point = (best_hit.ref_length - best_hit.ref_start) - 1 #interbase
				on_reverse_strand = True
		else:
			# try partial hits
			if(partial_hits_at_start and partial_hits_at_end):
				complete_dnaA_found, break_point, on_reverse_strand, dnaA_name = self._best_promer_partial_hits(partial_hits_at_start, partial_hits_at_end)
			
		return complete_dnaA_found, break_point, on_reverse_strand, dnaA_name
		
		
	def _is_complete_match_to_dnaa(self, hit, best_hit, contig_id):
		''' Is hit a complete match to dnaA '''
		if hit.ref_name == contig_id and \
		   hit.qry_start == 0 and \
		   hit.hit_length_qry >= (hit.qry_length * self.match_length_percent/100) and \
		   hit.percent_identity >= self.hit_percent_id and \
		   hit.hit_length_qry >= best_hit.hit_length_qry and \
		   hit.percent_identity >= best_hit.percent_identity:	 
		   return True
		return False 
	
		
	def _is_partial_match_at_start(self, hit, contig_id):
		''' Is hit a partial match to dnaA at the start of contig '''
		if hit.on_same_strand():
			if hit.ref_name == contig_id and \
		   	   hit.qry_end == hit.qry_length - 1 and \
			   hit.ref_start == 0 and \
			   hit.percent_identity >= self.hit_percent_id:	
			   return True
		else:
			if hit.ref_name == contig_id and \
		   	   hit.qry_end == 0 and \
			   hit.ref_start == 0 and \
			   hit.percent_identity >= self.hit_percent_id:	
			   return True	
		return False
		
	
	def _is_partial_match_at_end(self, hit, contig_id):
		''' Is hit a partial match to dnaA at the end of contig '''
		if hit.on_same_strand():
			if hit.ref_name == contig_id and \
			   hit.qry_start == 0 and \
			   hit.ref_end == hit.ref_length - 1 and \
			   hit.percent_identity >= self.hit_percent_id:
			   return True
		else:
			if hit.ref_name == contig_id and \
			   hit.qry_start == hit.qry_length - 1 and \
			   hit.ref_end == hit.ref_length - 1 and \
			   hit.percent_identity >= self.hit_percent_id:
			   return True
		return False
	
		
	def _best_promer_partial_hits(self, partial_hits_at_start, partial_hits_at_end):
		''' Look through partial hits for cases where dnaA is split across the start and end of contig'''
		complete_dnaA_found = False
		break_point = None
		on_reverse_strand = False
		dnaA_name = None
		best_start_hit = partial_hits_at_start[0]
		best_end_hit = partial_hits_at_end[0]
				
		for hit_at_start in partial_hits_at_start:
			for hit_at_end in partial_hits_at_end:
				if  hit_at_start.qry_name == hit_at_end.qry_name and \
					(hit_at_start.hit_length_qry + hit_at_end.hit_length_qry) == hit_at_start.qry_length and \
					hit_at_start.qry_length == hit_at_end.qry_length and \
					(hit_at_start.on_same_strand() is hit_at_end.on_same_strand()) and \
 					hit_at_start.percent_identity >= best_start_hit.percent_identity and \
 					hit_at_end.percent_identity >= best_end_hit.percent_identity:					
					best_start_hit = hit_at_start
					best_end_hit = hit_at_end
					complete_dnaA_found = True
		
		if complete_dnaA_found:

			dnaA_name = best_start_hit.qry_name
			if best_end_hit.on_same_strand():
				break_point = best_end_hit.ref_start						
			else:
				break_point = (best_end_hit.ref_length - best_end_hit.ref_start) - 1 #interbase
				on_reverse_strand = True
		return complete_dnaA_found, break_point, on_reverse_strand, dnaA_name
									
		
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
	
		
	def _write_summary(self, contig_id, break_point, gene_name, gene_reversed, new_name, skipped):
		'''Write summary'''
		if (not os.path.exists(self.summary_file)) or os.stat(self.summary_file).st_size == 0:
			header = self.summary_prefix + "\t" + '\t'.join(['id', 'break_point', 'gene_name', 'gene_reversed', 'new_name', 'skipped']) +'\n'
			utils.write_text_to_file(header, self.summary_file)
			
		breakpoint_to_print = break_point if break_point else '-'
		gene_name_to_print = gene_name if gene_name else '-'
		gene_reversed_to_print = '-'
		if gene_reversed:
			gene_reversed_to_print = 'yes'
		else:
			gene_reversed_to_print = 'no' if gene_name else '-'
		new_name_to_print = '-'
		if new_name and self.rename:
			new_name_to_print = new_name
		skipped_print = 'skipped' if skipped else '-'			
		text = self.summary_prefix + "\t" + '\t'.join(map(str, [contig_id, breakpoint_to_print, gene_name_to_print, gene_reversed_to_print, new_name_to_print, skipped_print])) + '\n'
		utils.write_text_to_file(text, self.summary_file)		
	

	def run(self):
		'''Look for break points in contigs'''
		contigs_in_file = set(self.contigs.keys())		
		if contigs_in_file != self.ids_to_skip:
			# run promer and prodigal only if needed
			dnaA_alignments = utils.run_nucmer(self.fasta_file, self.gene_file, self._build_promer_filename(), min_percent_id=self.hit_percent_id, run_promer=True)
			if self.choose_random_gene:
				self.random_gene_starts = self._run_prodigal_and_get_gene_starts()
		
		chromosome_count = 1
		plasmid_count = 1
		
		output_fw = fastaqutils.open_file_write(self.output_file)
		for contig_id in self.contigs:
			contig_sequence = self.contigs[contig_id]	
			gene_name = None
			new_name = contig_id 
			skipped = False
			break_point = None
			on_reverse_strand = False
			
			if contig_id not in self.ids_to_skip:				
				dnaA_found, break_point, on_reverse_strand, gene_name  = self._best_promer_hit(dnaA_alignments, contig_id)
				if dnaA_found:
					new_name = 'chromosome_' + str(chromosome_count)
					chromosome_count += 1		
				
				# Or look for a gene in prodigal results
				if not dnaA_found and self.choose_random_gene:
					if contig_id in self.random_gene_starts and self.random_gene_starts[contig_id]:
						gene_name = 'prodigal'
						if self.random_gene_starts[contig_id].strand == '+':			
							break_point = self.random_gene_starts[contig_id].start
						else:		
							break_point = (len(self.contigs[contig_id]) - self.random_gene_starts[contig_id].start) - 1 #interbase
							gene_on_reverse_strand = True
						new_name = 'plasmid_' + str(plasmid_count)
				
				# circularise the contig				
				if break_point:
					if on_reverse_strand:
						contig_sequence.revcomp()
					contig_sequence = contig_sequence[break_point:] + contig_sequence[0:break_point]
					self.contigs[contig_id].seq = contig_sequence
				
			else: # Skipped, just write contig as it is
				skipped = True
	
			# write the contig out			
			contig_name = new_name if self.rename else contig_id
			print(sequences.Fasta(contig_name, contig_sequence), file=output_fw)
			self._write_summary(contig_id, break_point, gene_name, on_reverse_strand, new_name, skipped)
			
		fastaqutils.close(output_fw)

		# clean up
		if not self.debug:
			utils.delete(self._build_promer_filename())
			utils.delete(self._build_prodigal_filename())
