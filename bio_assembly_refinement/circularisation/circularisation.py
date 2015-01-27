import os
from fastaq import tasks
from pymummer import coords_file, alignment, nucmer

class Circularisation:
	def __init__(self, 
				 fasta_file, 
				 output_file="tmp.circularised.fa", 
				 contigs={},
				 alignments=[],
				 acceptable_offset=12000, #Make this a % of length?
				 acceptable_length=0.5, 
				 percent_identity=99,
				 dnaa_sequence="file",
				 blast_program="farm_blast",
				 reassembling_program="quiver",				  
				 debug=False):
				 
		''' Constructor
		(complete)
		'''
		self.fasta_file = fasta_file
		self.output_file = output_file
		
		# If a set of contigs is not provided, read fasta file and get them
		self.contigs = contigs
		if not self.contigs:
			self.contigs = {}
			tasks.file_to_dict(self.fasta_file, self.contigs) #Read ids and sequences into hash

		# If a set of nucmer_hits if not provided, run nucmer and get alignments
		self.alignments = alignments
		if not self.alignments:
			results_file = str(self.fasta_file) + ".coords"
			runner = nucmer.Runner(self.fasta_file, self.fasta_file, results_file, coords_header=False, maxmatch=True) # nucmer default breaklength is 200
			runner.run()
			file_reader = coords_file.reader(results_file)
			self.alignments = [coord for coord in file_reader]
		
		self.acceptable_offset = acceptable_offset
		self.acceptable_length = acceptable_length
		self.percent_identity = percent_identity
		self.dnaa_sequence = dnaa_sequence,
		self.blast_program = blast_program,
		self.reassembling_program = reassembling_program,
		self.debug = debug
		
		
	def _circularisable(self, contig_id):
		'''
		Check if it's possible to circularise given contig
		TODO: Optimise. We go through each alignment in contigcleanup too
		How can that be reused?
		TODO: Move this check to alignment class
		
		1	2182	4783104	4780922	2182	2183	99.82	4791129	4791129	1	unitig_0|quiver	unitig_0|quiver
		'''
		
		for algn in self.alignments:		
			if algn.qry_name == contig_id and \
			   algn.ref_name == contig_id and \
			   algn.ref_start < self.acceptable_offset and \
			   algn.ref_end < (algn.ref_length * self.acceptable_length) and \
			   algn.qry_start > (algn.qry_length * self.acceptable_length) and \
			   algn.qry_end > (algn.qry_length - self.acceptable_offset) and \
			   algn.percent_identity > self.percent_identity:
				print(algn)
				return True  
		return False
				
				
	def run(self):
		'''
		'''
		sorted_contig_ids = sorted(self.contigs.keys()) #Sorting so that results are consistent each time
		possible_chromosomes = []
		for contig_id in sorted_contig_ids:
			if self._circularisable(contig_id):
 				print (contig_id)
 				# Blast using dnaa, cut and re-assemble

			   
		