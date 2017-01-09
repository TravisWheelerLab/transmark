#!/usr/bin/env perl
use strict;

#my $protfile = "rmarkinsertAA.msa";
#my $dnafile = "rmarkinsertORFandDNA.msa";

my ($protfile, $dnafile) = @ARGV;
 
if (not defined $protfile) {
  die "Need amino acid MSA as first input\n";
}

if (not defined $dnafile) {
  die "Need DNA MSA as second input\n";
}
 
my $seeds = "../Pfam-A.v27.seed";
unlink $protfile;

#get the names of the alignments in the DNA multiple alignment file
#these are  training (query) sequences
my @names = `esl-alistat $dnafile | grep "Alignment name" | awk '{print \$3}'`;
chomp @names;

foreach my $ali (@names) {
  print "$ali\n";
  #get the amino acid multiple alignment from the Pfam multiple alignments file
  do_cmd ("esl-afetch $seeds $ali > $ali.pfam.sto");
  #get the names of the sequences in the amino acid multiple alignment
  do_cmd ("esl-alistat --list $ali.pfam.names $ali.pfam.sto");

  #get the DNA multiple alignment 
  do_cmd ("esl-afetch  $dnafile $ali > $ali.dna.sto");
  #get the names of the sequences in the DNA multiple alignment; these are the query sequence
  #names for the homology search that will be run
  do_cmd ("esl-alistat --list $ali.dna.names $ali.dna.sto");
  #sometimes, the DNA names are formatted like "001|PUR5_SYMTH" (why!?).
  # clean this up:
#  do_cmd (qq[awk 'BEGIN{FS="|"} {print \$2}' $ali.dna.names > $ali.dna.names2]);
  do_cmd (qq[awk 'BEGIN{FS="|"} {if ( \$2 == \'\\n\' )  print \$1; else print \$2}' $ali.dna.names > $ali.dna.names2]);

  #get a list of (target) sequence names that are in the DNA multiple alignment
  #and also in the amino acid multiple alignment
  do_cmd ("grep -f $ali.dna.names2 $ali.pfam.names > $ali.filter.names");

  #remove all the amino acid sequences from the amino acid multiple alignment
  #that are not in the list of names of DNA query sequences 
  #from the DNA multple alignment
  #Now we have an amino acid multiple alignment that has all the sequences that are 
  #in the DNA multiple alignment
  #The query sequences will be used to search the benchmark for the target
  #sequences planted in the background sequence by rmark-create
  do_cmd ("esl-alimanip --seq-k $ali.filter.names $ali.pfam.sto >> $protfile");

  #`rm $ali.*`;
}


sub do_cmd {
   my $cmd = $_[0];
   #print "$cmd\n"; 
   `$cmd`;
}

