''' 
Find the point at which to break a circular contig (at origin of replication or random gene)

Attributes:
-----------

Sample usage:
-------------

'''

import os
import shutil
from bio_assembly_refinement import utils
from pyfastaq import sequences
from pyfastaq import utils as fastaqutils
from pymummer import alignment


class ContigBreakFinder:
	def __init__(self, 
			     fasta_file, 
			     gene_file, 
			     avoid=None, #Avoid circularising contigs with these ids
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
		
		self.ids_to_avoid = set()		
		if avoid:
			if isinstance(avoid, str) and os.path.isfile(avoid):			
				f = fastaqutils.open_file_read(avoid)
				for line in f:
					self.ids_to_avoid.add(line.rstrip())
				fastaqutils.close(f)
			else:
				self.ids_to_avoid = set(avoid) # Assumes ids_avoid is a list
		
		


	def _run_prodigal_and_get_start_of_a_gene(self, sequence):
		'''Run prodigal and find start of a gene''' 
		output_fw = fastaqutils.open_file_write("tmp_seq.fa")
		print(sequences.Fasta("contig", sequence), file=output_fw)
		output_fw.close()
		fastaqutils.syscall("prodigal -i tmp_seq.fa -o tmp_genes.gff -f gff -c")	
		boundary_start = round(0.3 * len(sequence)) # Look for a gene that starts after 30% of the sequence length 
		gene_start = 0
		fh = fastaqutils.open_file_read('tmp_genes.gff')
		for line in fh:
			if not line.startswith("#"):
				columns = line.split('\t')
				start_location = int(columns[3])
				if start_location > boundary_start:
					gene_start = start_location - 1 #Interbase
					break; 
		if not self.debug:   		
			utils.delete('tmp_genes.gff')
			utils.delete('tmp_seq.fa')
		return gene_start	
		
		
	def _build_final_filename(self):
		'''Build output filename'''
		input_filename = os.path.basename(self.fasta_file)
		return os.path.join(self.working_directory, "circularised_" + input_filename)
		
		
	def _write_summary(self, contig_id, break_point, gene_name, gene_reversed, new_name):
		'''Write summary'''
		if (not os.path.exists(self.summary_file)) or os.stat(self.summary_file).st_size == 0:
			header = '\t'.join(['id', 'break_point', 'gene_name', 'gene_reversed', 'new_name']) +'\n'
			utils.write_text_to_file(header, self.summary_file)
		name_to_print = new_name if self.rename else '-'			
		text = '\t'.join(map(str, [contig_id, break_point, gene_name, gene_reversed, name_to_print])) + '\n'
		utils.write_text_to_file(text, self.summary_file)		
	

	def run(self):
		'''Look for break point in contigs. If found, circularise and rename contig. Write to a log.'''	
		self.dnaA_alignments = utils.run_nucmer(self.fasta_file, self.gene_file, "gene.hits", min_percent_id=self.hit_percent_id, run_promer=True)
		plasmid_count = 1
		chromosome_count = 1
		seq_reader = sequences.file_reader(self.fasta_file)
		output_fw = fastaqutils.open_file_write(self.output_file)
		print("Ids to not circiularise")
		print(self.ids_to_avoid)
		for seq in seq_reader:
			if seq.id not in self.ids_to_avoid:		
				plasmid = True
				gene_name = '-'
				gene_on_reverse_strand = False
				new_name = ''
				break_point = 0
				
				for algn in self.dnaA_alignments:			
					if algn.ref_name == seq.id and \
					   algn.hit_length_qry >= (algn.qry_length * self.match_length_percent/100) and \
					   algn.percent_identity >= self.hit_percent_id and \
					   algn.qry_start == 0:	     
						plasmid = False
						gene_name = algn.qry_name
						if algn.on_same_strand():
							break_point = algn.ref_start						
						else:
							# Reverse complement sequence, circularise using new start of dnaA in the right orientation
							seq.seq = seq.seq.translate(str.maketrans("ATCGatcg","TAGCtagc"))[::-1]
							break_point = (algn.ref_length - algn.ref_start) - 1 #interbase
							dnaA_on_reverse_strand = True
					
						seq.seq = seq.seq[break_point:] + seq.seq[0:break_point]		
						new_name = 'chromosome' + str(chromosome_count)
						chromosome_count += 1		
						break;
					
				if plasmid:
					if len(seq.seq) < 200000: # Only rename if it's roughly plasmid size
						new_name = 'plasmid' + str(plasmid_count)
						plasmid_count += 1
					if self.choose_random_gene and len(seq.seq) > 20000:
						# If suitable, choose random gene in plasmid, and circularise. Prodigal only works for contigs > 20000 bases				
						gene_start = self._run_prodigal_and_get_start_of_a_gene(seq.seq)
						if gene_start:
							seq.seq = seq.seq[gene_start:] + seq.seq[0:gene_start]	
							break_point = gene_start
							gene_name = 'predicted_gene'
					
				contig_name = new_name if self.rename else seq.id
				print(sequences.Fasta(contig_name, seq.seq), file=output_fw)
				self._write_summary(seq.id, break_point, gene_name, gene_on_reverse_strand, new_name)
				
		output_fw.close()
		
		if not self.debug:
			utils.delete("nucmer_gene.hits")
