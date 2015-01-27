#!/usr/bin/perl

use warnings;
use strict;

while (<>) {
    unless ($_ =~ /^(\s*)\d/){
	next
    }
    $_ =~ s/\|//g;
   
    my @f = split;
    # create crude match score = ((length_of_match * %identity)-(length_of_match * (100 - %identity))) /20
    my $crude_plus_score=($f[4]*$f[6]);
    my $crude_minus_score=($f[4]*(100-$f[6])); 
    my $crude_score=  int(($crude_plus_score  - $crude_minus_score) / 20);      
    
    # reorganise columns and print crunch format to stdout
    # score        %id   S1    E1    seq1  S2    E2    seq2  (description)
    
    print " $crude_score $f[6] $f[0] $f[1] $f[7] $f[2] $f[3] $f[8] nucmer comparison coordinates\n"
}

