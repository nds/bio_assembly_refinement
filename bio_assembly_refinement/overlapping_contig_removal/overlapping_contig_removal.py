from fastaq import tasks

def remove(infile, ids_file, outfile):
#   if not outfile:
# 	    outfile = 'filtered_' + infile
	tasks.filter(infile, outfile, ids_file=ids_file, invert=True)