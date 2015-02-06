bio\_assembly\_refinement
=======================

Modules to filter, circularise and re-assemble contigs, mostly useful for (but not limited to) bacterial assemblies

Description
-----------

Given a fasta file, these modules can be used to:

1. Filter out contigs smaller than a set length, and those completely contained within other contigs
2. Identify contigs that can be circularised (by checking if their ends overlap). These contigs can then be trimmed and circularised, setting the site of the dnaA/refA/refB gene as the start
3. Run pacbio smrtanalysis software on the new assembly to smooth out the new joins in the contigs

Each step can be run individually, or as a collection. There is also a script that can be used to invoke this functionality on the command line

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

**overlap\_offset**: offset from the ends of a contig where an overlap region can begin (expressed as % of length of contig) (default 12)

**overlap\_max\_length**: maximum length of the overlap between ends (expressed as % of contig length) (default 50)

**overlap\_percent\_identity**: minimum percent identity in nucmer hits to use when determining if ends overlap (default 99)

**dnaA\_hit\_percent\_identity**: minimum percent identity to consider when looking at hits to dnaA/refA/refB (default 99)

**dnaA\_hit\_length\_minimum**: minimum length of hit to dnaA/refA/refB expressed as % of length of dnaA/refA/refB (default 95)

**working\_directory**: working directory (default current working directory) 

**pacbio\_exec**: pacbio resequencing exec (default pacbio_smrtanalysis) 

**nucmer\_exec**: nucmer exec (default nucmer) 

**debug**: do not delete temp files if set to true (default false)

More documentation available in the code.


Contact
-------

Author: Nishadi De Silva

Affiliation: Wellcome Trust Sanger Institute, Hinxton, UK

Email: path-help@sanger.ac.uk
      
