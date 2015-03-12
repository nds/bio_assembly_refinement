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
		test_contigs = { "contig1": "TATCGAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTAT", # 4 base overlap and dnaA present
					 	 "contig2": "AGTCCACCGGGCACTGCAAGGTAAATTCTTACGCCCACTTTGTAGACCCTACCGTAAAGC", # Overlap smaller than cutoff (should be ignored)
					 	 "contig3": "AAGTCGGTAGCTTGGCATACGTTACGTATTACGCCCACTTTGTAGAAACTACCGGCTGAA", # 6 base overlap but no dnaA (plasmid)
					 	 "contig4": "TATCGAGTATATTATCAACTGAGGCGGTCCCGACGCATATAATTATGAAAATGGCTCTAT", # 4 base overlap and dnaA on negative strand (rev comp sequence)
					 	 "contig5": "GCTTATCGAGTATAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTATCTA", # 4 base overlap but not right at the end, dnaA present
					 	  
					   }	
					   			   		
		test_overlap_coords = ['\t'.join(['1', '4', '57', '60', '4', '4', '100.00', '60', '60', '1', '1', 'contig1', 'contig1']),
				 	 		   '\t'.join(['7', '9', '52', '54', '3', '3', '100.00', '60', '60', '1', '1', 'contig2', 'contig2']),
				 	 		   '\t'.join(['1', '6', '55', '60', '6', '6', '100.00', '60', '60', '1', '1', 'contig3', 'contig3']),
				 	 		   '\t'.join(['1', '4', '57', '60', '4', '4', '100.00', '60', '60', '1', '1', 'contig4', 'contig4']),
				 	 		   '\t'.join(['4', '7', '57', '54', '4', '4', '100.00', '60', '60', '1', '1', 'contig5', 'contig5'])
							  ]							  
		test_overlap_alignments = [alignment.Alignment(coord) for coord in test_overlap_coords] 
		
		# Trimmed sequences
		#"contig1": "GAGTATATTATCAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTAT", # had 4 base overlap, dnaA present
		#"contig3": "GTAGCTTGGCATACGTTACGTATTACGCCCACTTTGTAGAAACTACCGGCTGAA", # had 6 base overlap, dnaA absent (plasmid)
		#"contig4": "GAGTATATTATCAACTGAGGCGGTCCCGACGCATATAATTATGAAAATGGCTCTAT", # had 4 base overlap, dnaA on negative strand
		#           "ATAGAGCCATTTTCATAATTATATGCGTCGGGACCGCCTCAGTTGATAATATACTC" # Revcom contig 4 sequence
		#"contig5": "GAGTATAACTGGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTAT", # had 4 base overlap but not right at the end, dnaA present	
		  		
		# Test dnaA: GGACCGCCTC (revcomp GAGGCGGTCC)
		# Simulated hits to dnaA using trimmed sequences
		test_dnaA_coords = ['\t'.join(['17', '26', '1', '10', '10', '10', '100.00', '56', '10', '1', '1', 'contig1', 'dnaA']),
							'\t'.join(['26', '17', '1', '10', '10', '10', '100.00', '56', '10', '-1', '1', 'contig4', 'dnaA']),
							'\t'.join(['11', '20', '1', '10', '10', '10', '100.00', '56', '10', '1', '1', 'contig5', 'dnaA'])
		
		] 
		test_dnaA_alignments = [alignment.Alignment(coord) for coord in test_dnaA_coords]
		
		# Circularised sequences
		circularised_contigs = {"contig1": "GGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTATGAGTATATTATCAACT", # chromosome1
								"contig3": "GTAGCTTGGCATACGTTACGTATTACGCCCACTTTGTAGAAACTACCGGCTGAA", # plasmid1
								"contig4": "GGACCGCCTCAGTTGATAATATACTCATAGAGCCATTTTCATAATTATATGCGTCG", # chromosome2
								"contig5": "GGACCGCCTCCGACGCATATAATTATGAAAATGGCTCTATGAGTATAACT", # chromosome3	
								}
		# Output data
		expected_output = os.path.join(data_dir, "circularised_test_contigs.fa")
		actual_output = os.path.join(os.getcwd(), "circularised_file.fa")
		summary_file = os.path.join(os.getcwd(), 'circularisation_summary_file.txt')
	   	
		circulariser = circularisation.Circularisation(dnaA_sequence = test_dnaA_file,
													   contigs=test_contigs,
												       alignments = test_overlap_alignments,
												       dnaA_alignments = test_dnaA_alignments,
												       overlap_offset = 9,
												       overlap_min_length=3,
												       debug = False											       
												      )
	
		circulariser.run()
		output_contigs = circulariser.contigs
		
				
# 		self.assertTrue(filecmp.cmp(actual_output, expected_output, shallow=False)) 
# 		If the contigs are of equal size there is no guarantee as to which order they will be put in by fastaq sort. So, checking contig sequences one by one.
		
		for id in circularised_contigs.keys():
			self.assertEqual(circularised_contigs[id], output_contigs[id])
		self.assertTrue(os.path.isfile(actual_output))# Does output file exist and is it named right? 
		self.assertTrue(os.path.isfile(summary_file))# Does summary file exist and is it named right? 
		os.remove(actual_output)
		os.remove(summary_file)
		
		
	def test_circularisation_with_fasta_file(self):	
		# Does constructor take fasta file?
		input_file = os.path.join(data_dir, "Salmonella_pacbio_unitig_0.fa")
		dnaA_file = os.path.join(data_dir, "dnaA.fa")
		output_file = os.path.join(os.getcwd(), "circularised_Salmonella_pacbio_unitig_0.fa")	
		summary_file = os.path.join(os.getcwd(), 'circularisation_summary_file.txt')	
		circulariser = circularisation.Circularisation(dnaA_sequence = dnaA_file, fasta_file=input_file)
		circulariser.run()
		self.assertTrue(os.path.isfile(output_file))# Does output file exist and is it named right? 
		self.assertTrue(os.path.isfile(summary_file))# Does summary file exist and is it named right? 
		os.remove(output_file)
		os.remove(circulariser.summary_file)