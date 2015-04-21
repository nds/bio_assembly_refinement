import unittest
import filecmp
import os
from bio_assembly_refinement import contig_overlap_trimmer 
from pymummer import alignment

modules_dir = os.path.dirname(os.path.abspath(contig_overlap_trimmer.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestContigOverlapTrimmer(unittest.TestCase):
	def test_contig_overlap_trimmer(self):
		'''Test contig overlap trimming'''		
		# test data
		test_fasta_file = os.path.join(data_dir, "TRIMMING_input_1.fa")
		test_overlap_coords = ['\t'.join(['1', '4', '57', '60', '4', '4', '100.00', '60', '60', '1', '1', 'contig1', 'contig1']),
							   '\t'.join(['1', '4', '57', '60', '4', '4', '100.00', '60', '60', '1', '1', 'contig2', 'contig2']),
				 	 		   '\t'.join(['2', '4', '57', '59', '3', '3', '100.00', '60', '60', '1', '1', 'contig3', 'contig3']),	
				 	 		   '\t'.join(['2', '4', '54', '56', '3', '3', '100.00', '60', '60', '1', '1', 'contig4', 'contig4']),
				 	 		   '\t'.join(['1', '3', '58', '60', '3', '3', '100.00', '60', '60', '1', '1', 'contig4', 'contig4']),		 	 		   
				 	 		   '\t'.join(['1', '4', '57', '60', '4', '4', '100.00', '60', '60', '1', '1', 'contig5', 'contig5']),				 	 		   
				 	 		   '\t'.join(['1', '2', '59', '60', '2', '2', '100.00', '60', '60', '1', '1', 'contig6', 'contig6']),				 	 		   
				 	 		   '\t'.join(['1', '12', '49', '60', '12', '12', '100.00', '60', '60', '1', '1', 'contig7', 'contig7']),				 	 		   
				 	 		   '\t'.join(['13', '15', '57', '54', '4', '4', '100.00', '60', '60', '1', '-1', 'contig8', 'contig8']),
				 	 		   # No overlap for contig 9
				 	 		   '\t'.join(['4', '7', '36', '38', '4', '4', '100.00', '60', '60', '1', '-1', 'contig10', 'contig10']),
							  ]							  
		test_overlap_alignments = [alignment.Alignment(coord) for coord in test_overlap_coords] 
		overlap_trimmer = contig_overlap_trimmer.ContigOverlapTrimmer(fasta_file = test_fasta_file,
														     		  alignments = test_overlap_alignments,
														     		  overlap_offset = 10,
												       				  overlap_min_length=3,
												                      overlap_max_length=12,										     								       
												                     )	
		overlap_trimmer.run()
		