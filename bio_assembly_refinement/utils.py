'''
Helper functions like deleting files and running nucmer
'''

import os
import sys
from pymummer import coords_file, alignment, nucmer

class Error (Exception): pass

def delete(filename):
	''' Delete a file if it exists '''
	if os.path.exists(filename):
		try: 
			os.remove(filename)
		except OSError as e:
			raise Error("Error deleting file '" + e.filename + "'")
			
			
def run_nucmer(ref, query, output):
	'''Run nucmer and return a list of alignment objects'''
	runner = nucmer.Runner(ref, query, output, coords_header=False, maxmatch=True, simplify=False) # nucmer default break length is 200
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

	