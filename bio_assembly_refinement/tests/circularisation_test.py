import unittest
import filecmp
import os
from bio_assembly_refinement import circularisation 
from pymummer import alignment

modules_dir = os.path.dirname(os.path.abspath(circularisation.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestCircularisation(unittest.TestCase):
	def test_circularisation_with_precomputed_values(self):
		'''Test _circularisation_with_precomputed_values'''		
		# Input data
		test_dnaA_file = os.path.join(data_dir, "test_dnaA.fa")
		test_contigs = { 
						 # Chromosomes 	
						 "contig1": "TATCGAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTAT", # 4 base overlap at the very ends, dnaA present
						 "contig2": "TATCGAGTATATTATCAACTGAGGCGGTCCCGACGCATATAATTATGAAAATGGCTCTAT", # 4 base overlap at the very ends, dnaA reverse orientation
						 "contig3": "TATCGAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTAG", # 3 base overlap but not starting at the beginning, dnaA present	
						 "contig4": "TATCGAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGCTACTAT", # Two 3 base overlaps, dnaA present						 
						 # Plasmids
						 "contig5": "TATCGAGTATATTATCAACTGTACCCCCTCCGACGCATATAATTATGAAAATGGCTCTAT", # 4 base overlap at the very ends, dnaA absent
						 # Do not circularise
					 	 "contig6": "AGTCCACCGGGCACTGCAAGGTAAATTCTTACGCCCACTTTGTAGACCCTACCGTAAAGC", # Overlap smaller than min length of 3 
					 	 "contig7": "AAGTCGGTAGCTTGGCATACGTTACGTATTACGCCCACTTTGTAGAAATCGATGGCTGAA", # 6 base overlap, no dnaA (trimmed length would be too small)					 	
					 	 "contig8": "GCTTATCGAGTATAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTATCTA", # 4 base overlap but reversed (ignore)
					 	 "contig9": "AGTCCACCGGGCACTGCAAGGTAAATTCTTACGCCCACTTTGTAGACCCTACCGTAAAGC", # No overlap
					 	 "contig10": "AGTCCACCGGGCACTGCAAGGTAAATTCTTACGCCCACTTTGTAGACCCTACCGTAAAGC", # Overlap start beyond offset				 	  
					   }	
					   			   		
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
		
		# Trimmed sequences
		#"contig1": "GAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTAT", # 4 base overlap at the very ends, dnaA present
		#"contig2": "GAGTATATTATCAACTGAGGCGGTCCCGACGCATATAATTATGAAAATGGCTCTAT", # 4 base overlap at the very ends, dnaA reverse orientation
		#"contig3": "GAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTA", # 3 base overlap but not starting at the beginning, dnaA present
		#"contig4": "CGAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGCTACTAT", # Two 3 base overlaps, dnaA present						 
		# Plasmids
		#"contig5": "GAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTAT", # 4 base overlap at the very ends, dnaA absent
	
		# Test dnaA: GGACCGCCTC (revcomp GAGGCGGTCC)
		# Simulated hits to dnaA using trimmed sequences
		test_dnaA_coords = ['\t'.join(['17', '26', '1', '10', '10', '10', '100.00', '56', '10', '1', '1', 'contig1', 'dnaA']),
							'\t'.join(['26', '17', '1', '10', '10', '10', '100.00', '56', '10', '-1', '1', 'contig2', 'dnaA']),
							'\t'.join(['17', '26', '1', '10', '10', '10', '100.00', '55', '10', '1', '1', 'contig3', 'dnaA']),
							'\t'.join(['18', '27', '1', '10', '10', '10', '100.00', '57', '10', '1', '1', 'contig4', 'dnaA'])
		
		]  
		test_dnaA_alignments = [alignment.Alignment(coord) for coord in test_dnaA_coords]
		# Circularised sequences
		circularised_contigs = {"contig1": "GGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTATGAGTATATTATCAACT", # chromosome1		 
								"contig2": "GGACCGCCTCAGTTGATAATATACTCATAGAGCCATTTTCATAATTATATGCGTCG", # chromosome2
								"contig3": "GGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTAGAGTATATTATCAACT", # chromosome3
								"contig4": "GGACCGCCTCCGACGCATATAATTATGAAAATGCTACTATCGAGTATATTATCAACT", # chromosome4
								"contig5": "GAGTATATTATCAACTGTACCCCCTCCGACGCATATAATTATGAAAATGGCTCTAT", # plasmid1	
								
								# Following contigs should not have changed
								"contig6": "AGTCCACCGGGCACTGCAAGGTAAATTCTTACGCCCACTTTGTAGACCCTACCGTAAAGC", # Overlap smaller than min length of 3 
					 	 		"contig7": "AAGTCGGTAGCTTGGCATACGTTACGTATTACGCCCACTTTGTAGAAATCGATGGCTGAA", # 6 base overlap, no dnaA (trimmed length would be too small)					 	
					 			"contig8": "GCTTATCGAGTATAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTATCTA", # 4 base overlap but reversed (ignore)
					 	 		"contig9": "AGTCCACCGGGCACTGCAAGGTAAATTCTTACGCCCACTTTGTAGACCCTACCGTAAAGC", # No overlap
					 	 		"contig10": "AGTCCACCGGGCACTGCAAGGTAAATTCTTACGCCCACTTTGTAGACCCTACCGTAAAGC", # Overlap start beyond offset				 	  
								}
		# Output data
		expected_output = os.path.join(data_dir, "circularised_test_contigs.fa")
		actual_output = os.path.join(os.getcwd(), "circularised_file.fa")
		summary_file = os.path.join(os.getcwd(), 'circularisation_summary_file.txt')
	   	
		circulariser = circularisation.Circularisation(dnaA_sequence = test_dnaA_file,
													   contigs=test_contigs,
												       alignments = test_overlap_alignments,
												       dnaA_alignments = test_dnaA_alignments,
												       overlap_offset = 10,
												       overlap_min_length=3,
												       overlap_max_length=12,
												       debug = False											       
												      )
	
		circulariser.run()
		output_contigs = circulariser.contigs
		
		for id in circularised_contigs.keys():		
			self.assertEqual(circularised_contigs[id], output_contigs[id])
		self.assertTrue(os.path.isfile(actual_output))# Does output file exist and is it named right? 
		self.assertTrue(os.path.isfile(summary_file))# Does summary file exist and is it named right? 
		os.remove(actual_output)
		os.remove(summary_file)
		