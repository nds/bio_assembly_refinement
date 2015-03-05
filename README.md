bio\_assembly\_refinement
=======================

Modules to filter, circularise and re-assemble contigs, mostly useful for (but not limited to) bacterial assemblies

Description
-----------

Given a fasta file, these modules can be used to:

1. Filter out contigs smaller than a set length, and those completely contained within other contigs
2. Circularise contigs by trimming overlaps and setting the start of the conserved dnaA gene (or other genes in the case of plasmids) as the start of the sequence
3. Run quiver (pacbio_smrtanalysis script) on the new circularised contigs to refine them, particularly around the new joins

There is also a script that can be used to invoke this functionality on the command line

Installation
------------

###Pre-requisites###

__1.	MUMmer__

Instructions to install MUMmer can be found [here](http://mummer.sourceforge.net/manual/#installation)
    
__2.	pacbio\_smrtanalysis RS\_sequencing__
	
[Installation instructions to be completed]

__3.	pyfastaq__ 
		
Install: 
	
	pip3 install pyfastaq
		
__4.	prodigal__ 
		
PRODIGAL gene finding software. [Installation instructions](https://github.com/hyattpd/prodigal/wiki/Installation) for PRODIGAL 


Usage (for developers)
----------------------

Sample usage of main module:

	from bio_assembly_refinement import main 
	processor = main.Main(fasta_file = input_file, 
					      dnaA_sequence = dnaA_file,
					      bax_files = data_dir
					   	 )
	processor.process_assembly()


Attributes of Main.py:
----------------------

**fasta\_file**: input fasta file

**dnaA\_sequence**: fasta file with dnaA/refA/refB sequences
 
**bax\_files**: directory containing bax.h5 files

**cutoff\_contig\_length**: minimum contig length (default 10000)

**contained\_percent\_match**: minimum percent identity in nucmer hits to determine if contig is contained in another (default 95)

**overlap\_offset**: offset from the ends of a contig where an overlap region can begin (expressed as % of length of contig) (default 49)

**overlap\_boundary\_max**: maximum boundary of the overlap between ends (expressed as % of contig length) (default 50)

**overlap\_min\_length**: minimum length of overlap (default 2000 bases)

**overlap\_percent\_identity**: minimum percent identity in nucmer hits to use when determining if ends overlap (default 85)

**dnaA\_hit\_percent\_identity**: minimum percent identity to consider when looking at hits to dnaA/refA/refB (default 80)

**dnaA\_hit\_length\_minimum**: minimum length of hit to dnaA/refA/refB expressed as % of length of dnaA/refA/refB (default 95)

**working\_directory**: working directory (default current working directory) 

**pacbio\_exec**: pacbio resequencing exec (default pacbio_smrtanalysis) 

**nucmer\_exec**: nucmer exec (default nucmer) 

**reassembly\_dir**: directory sent to quiver (default reassembly)

**summary\_file**: summary file (default pacbio\_postprocess\_summary.txt) 

**debug**: do not delete temp files if set to true (default false)

More documentation available in the code.


Contact
-------

Author: Nishadi De Silva

Affiliation: Wellcome Trust Sanger Institute, Hinxton, UK

Email: path-help@sanger.ac.uk
      
