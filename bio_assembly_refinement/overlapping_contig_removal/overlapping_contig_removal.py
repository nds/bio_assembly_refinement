# Run nucmer on given fasta file (against itself), analyse results and remove any overlapping
# contigs
# The rules to define what 'overlapping' means are encoded below
# In the script, leave an option to retain all intermediate steps
# If not, delete the temp directory
# TODO: Catch errors

from fastaq import tasks
from pymummer import coords_file, alignment, nucmer

def find_hits(fastafile):
	'''Currently runs nucmer and returns a list of alignment objects'''
	resultsfile = str(fastafile) + ".coords"
	runner = nucmer.Runner(fastafile, fastafile, resultsfile, coords_header=False)
	runner.run()
	filereader = coords_file.reader(resultsfile)
	alignments = [coord for coord in filereader]
	return alignments

# Parse alignments and identify overlapping hits
def find_overlapping_contigs(alignments):
	ids_to_remove = []
	for algn in alignments:
		if algn.is_self_hit():
			ids_to_remove.append(algn.ref_name)
			
	return ids_to_remove
			
			
	
	
# 	
# 	def _contig_contained_in_nucmer_hits(self, hits, contig, min_percent):
#         assert contig in self.contigs
#         contig_length = len(self.contigs[contig])
#         coords = []
#         for hit in [hit for hit in hits if contig in [hit.qry_name, hit.ref_name] and hit.qry_name != hit.ref_name]:
#             start = min(hit.qry_start, hit.qry_end)
#             end = max(hit.qry_start, hit.qry_end)
#             coords.append(fastaq.intervals.Interval(start, end))
# 
#         if len(coords) == 0:
#             return False
# 
#         fastaq.intervals.merge_overlapping_in_list(coords)
#         total_bases_matched = fastaq.intervals.length_sum_from_list(coords)
#         return min_percent <= 100.0 * total_bases_matched / len(self.contigs[contig])
	

def remove_overlapping_contigs(infile, ids_file, outfile):
	tasks.filter(infile, outfile, ids_file=ids_file, invert=True)