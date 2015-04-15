'''
Helper functions like deleting files, running nucmer and finding break point
'''

import os
import sys
from pymummer import coords_file, alignment, nucmer
from pyfastaq import utils as fastaqutils, sequences

class Error (Exception): pass

def delete(filename):
	''' Delete a file if it exists '''
	if os.path.exists(filename):
		try: 
			os.remove(filename)
		except OSError as e:
			raise Error("Error deleting file '" + e.filename + "'")
			
			
def run_nucmer(ref, query, output, min_percent_id=95, run_promer=False):
	'''Run nucmer and return a list of alignment objects'''
	runner = nucmer.Runner(ref, query, output, min_id=min_percent_id, coords_header=False, maxmatch=True, simplify=False, promer=True) # nucmer default break length is 200
	runner.run()
	file_reader = coords_file.reader(output)
	alignments = [coord for coord in file_reader]
	return alignments
	
	
def write_ids_to_file(ids, filename):
	'''Write contig ids to a file'''
	with open(filename, mode='wt') as ids_file:
		ids_file.write('\n'.join(ids))
	ids_file.close()
	return filename
	
	
def write_text_to_file(text, filename):
	'''Append text to file'''
	with open(filename, mode='a') as text_file:
		text_file.write(text)
	text_file.close()
	return filename
	
	
def run_prodigal_and_get_start_of_a_gene(sequence):
	# Write plasmid to a file
	# Run prodigal
	# Return a gene start location (first, or middle?)
	# Delete plasmid file and the gff file output by prodigal
	
	output_fw = fastaqutils.open_file_write("tmp_seq.fa")
	print(sequences.Fasta("contig", sequence), file=output_fw)
	output_fw.close()
	fastaqutils.syscall("prodigal -i tmp_seq.fa -o tmp_genes.gff -f gff -c")
	
	boundary_start = round(0.3 * len(sequence)) # Look for a gene that starts after 30% of the sequence length (i.e. be sure to avoid region of overlap)
	gene_start = 0
	
	fh = fastaqutils.open_file_read('tmp_genes.gff')
	for line in fh:
		if not line.startswith("#"):
			columns = line.split('\t')
			start_location = int(columns[3])
			if start_location > boundary_start:
				gene_start = start_location - 1 #Interbase
				break;    		
	delete('tmp_genes.gff')
	delete('tmp_seq.fa')
	return gene_start		
	

