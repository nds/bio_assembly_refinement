'''
Helper functions like deleting files, running nucmer and finding break point
'''

import os
import sys
from pymummer import coords_file, alignment, nucmer
from pyfastaq import utils as fastaqutils, sequences
import re
import subprocess
from distutils.version import LooseVersion

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
	runner = nucmer.Runner(ref, query, output, min_id=min_percent_id, coords_header=False, maxmatch=True, simplify=False, promer=run_promer) # nucmer default break length is 200
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
	
	
def parse_file_or_set(s):
	'''Parse a file or set and return set of items in it '''
	items = set()	
	if s:
		if type(s) == set:
			items = s
		else:
			fh = fastaqutils.open_file_read(s) #Will just fail is file not found. Handle properly
			for line in fh:
				items.add(line.rstrip())
			fastaqutils.close(fh)
	return items
	
	
def get_prodigal_version():
	'''Get prodigal version''' #There must be a better way to do this!
	cmd = "prodigal -v"
	cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	for l in cmd_output:
		l = l.decode()
		m = re.findall(r'[V](\d\.\d)', l)
	if m:
		return m[0]
		
		
def run_prodigal(input_fasta, prodigal_output_file, length):
	'''Run prodigal - check version, check errors'''
	try:
		p_option = ""
		if (length < 20000):
			# prodigal needs -p meta option for sequences less than 20000
			# annoyingly newer version of prodigal has different -p option!
			version = get_prodigal_version()
			if LooseVersion(version) >= LooseVersion('3.0'):
				p_option = "-p anon"
			else:
				p_option = "-p meta"
		cmd = "prodigal -i " + input_fasta + " -o " + prodigal_output_file + " -f gff -c -q " + p_option
		output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
		return prodigal_output_file
	except  subprocess.CalledProcessError as error:
		print('Error running prodigal, so will not use gene predictions \n')
		print(error.output.decode())
		return None
        

	

