import unittest
import filecmp
import os
from bio_assembly_refinement import contig_break_finder 
from pyfastaq import tasks


modules_dir = os.path.dirname(os.path.abspath(contig_break_finder.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestContigBreakFinder(unittest.TestCase):
	def test_finding_dnaA_in_various_positions(self):
		
		tests = [
			#dnaa at start - return identical sequence
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKFINDER_input_dnaa_at_start.fa"), 
							      											 gene_file = os.path.join(data_dir, "BREAKFINDER_test_dnaA.fa"), 
							      											 choose_random_gene=False, 
							      											 rename = False,	
							      											), 
													'BREAKFINDER_output_dnaa_at_start.fa' ], 
			
			# dnaa in the middle
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKFINDER_input_dnaa_in_middle.fa"), 
							      											 gene_file = os.path.join(data_dir, "BREAKFINDER_test_dnaA.fa"), 
							      											 choose_random_gene=False, 
							      											 rename = False,	
							      											), 
													'BREAKFINDER_output_dnaa_in_middle.fa' ], 
			
			# dnaa at the end
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKFINDER_input_dnaa_at_end.fa"), 
							      											 gene_file = os.path.join(data_dir, "BREAKFINDER_test_dnaA.fa"), 
							      											 choose_random_gene=False, 
							      											 rename = False,	
							      											), 
													'BREAKFINDER_output_dnaa_at_end.fa' ], 
													
			# dnaa split across start and end
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKFINDER_input_dnaa_split.fa"), 
							      											 gene_file = os.path.join(data_dir, "BREAKFINDER_test_dnaA.fa"), 
							      											 choose_random_gene=False, 
							      											 rename = False,	
							      											), 
													'BREAKFINDER_output_dnaa_split.fa' ], 
			# dnaa in middle of contig  but revcom
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKFINDER_input_dnaa_in_middle_revcom.fa"), 
							      											 gene_file = os.path.join(data_dir, "BREAKFINDER_test_dnaA.fa"), 
							      											 choose_random_gene=False, 
							      											 rename = False,	
							      											), 
													'BREAKFINDER_output_dnaa_in_middle.fa' ], 											
			# dnaa split across start and end but revcom
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKFINDER_input_dnaa_split_revcom.fa"), 
							      											 gene_file = os.path.join(data_dir, "BREAKFINDER_test_dnaA.fa"), 
							      											 choose_random_gene=False, 
							      											 rename = False,	
							      											), 
													'BREAKFINDER_output_dnaa_split_revcom.fa' ], 
			# no dnaa 
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKFINDER_input_no_dnaa.fa"), 
							      											 gene_file = os.path.join(data_dir, "BREAKFINDER_test_dnaA.fa"), 
							      											 choose_random_gene=False, 
							      											 rename = False,	
							      											), 
													'BREAKFINDER_input_no_dnaa.fa' ], #do not change the contig
			# best dnaa hit not first 
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKFINDER_input_multiple_dnaa.fa"), 
							      											 gene_file = os.path.join(data_dir, "BREAKFINDER_test_multiple_dnaA.fa"), 
							      											 choose_random_gene=False, 
							      											 rename = False,	
							      											), 
													'BREAKFINDER_output_multiple_dnaa.fa' ], 			

 			]
				
		for t in tests:
			t[0].run()
			self.assertTrue(os.path.isfile(t[0].output_file))
			self.assertTrue(os.path.isfile(t[0].summary_file))
			expected_contigs = {}
			tasks.file_to_dict(os.path.join(data_dir, t[1]), expected_contigs)
			for id in expected_contigs.keys():
				self.assertTrue(expected_contigs[id].seq == t[0].contigs[id].seq)
			os.remove(t[0].output_file)
			os.remove(t[0].summary_file)	



	def test_options(self):	
		tests_for_options = [
			#rename genes
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKFINDER_input_dnaa_at_start.fa"), 
							      											 gene_file = os.path.join(data_dir, "BREAKFINDER_test_dnaA.fa"), 
							      											 choose_random_gene=False, 
							      											 rename = True,	
							      											), 
													'BREAKFINDER_output_dnaa_at_start.fa' ], 	
			# no dnaa, but use prodigal	
 			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKFINDER_input_no_dnaa_use_prodigal.fa"), 
 							      											 gene_file = os.path.join(data_dir, "BREAKFINDER_real_dnaa.fa"), 
 							      											 choose_random_gene=True, 
 							      											 rename = False,	
 							      											), 
 													'BREAKFINDER_output_no_dnaa_use_prodigal.fa' ], 										
			# skip one contig
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKFINDER_input_multiple_contigs.fa"), 
							      											 gene_file = os.path.join(data_dir, "BREAKFINDER_test_dnaA.fa"), 
							      											 choose_random_gene=False, 
							      											 rename = False,	
							      											 skip = os.path.join(data_dir, "BREAKFINDER_skip_one_id.txt")
							      											), 
													'BREAKFINDER_output_skip_contig.fa' ], 				
			# skip all contigs	
			[contig_break_finder.ContigBreakFinder(fasta_file = os.path.join(data_dir, "BREAKFINDER_input_multiple_contigs.fa"), 
							      											 gene_file = os.path.join(data_dir, "BREAKFINDER_test_dnaA.fa"), 
							      											 choose_random_gene=False, 
							      											 rename = True,	
							      											 skip = os.path.join(data_dir, "BREAKFINDER_skip_all.txt")
							      											), 
													'BREAKFINDER_input_multiple_contigs.fa' ], 	# do not change anything							
			]
			
		for t in tests_for_options:
			t[0].run()
			expected_contigs = {}
			tasks.file_to_dict(os.path.join(data_dir, t[1]), expected_contigs)
			for id in expected_contigs.keys():
				self.assertTrue(expected_contigs[id].seq == t[0].contigs[id].seq)
			os.remove(t[0].output_file)
			os.remove(t[0].summary_file)	
	
	
