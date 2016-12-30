#!/usr/bin/env perl
use strict;

my $protfile = "rmarkinsertAA.msa";
my $dnafile = "rmarkinsertORFandDNA.msa";
my $seeds = "../Pfam-A.v27.seed";
unlink $protfile;

#get the names of the alignments in the DNA multiple alignment file
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
  #get the names of the sequences in the DNA multiple alignment; these are the taget sequences
  #for the homology search that will be run
  do_cmd ("esl-alistat --list $ali.dna.names $ali.dna.sto");
  #sometimes, the DNA names are formatted like "001|PUR5_SYMTH" (why!?).
  # clean this up:
#  do_cmd (qq[awk 'BEGIN{FS="|"} {print \$2}' $ali.dna.names > $ali.dna.names2]);
  do_cmd (qq[awk 'BEGIN{FS="|"} {if ( \$2 == \'\\n\' )  print \$1; else print \$2}' $ali.dna.names > $ali.dna.names2]);

  #get a list of (target) sequence names that are in the DNA multiple alignment
  #and also in the amino acid multiple alignment
  do_cmd ("grep -f $ali.dna.names2 $ali.pfam.names > $ali.filter.names");

  #remove all the amino acid sequences from the amino acid multiple alignement
  #that except those sequences that are in the list of names of (target) sequences 
  #from the DNA multple alignement
  #Now we have an amino acid multiple alignement that has sequences that are 
  #the query sequences. These sequences are not in the DNA multiple alignement
  #which contains the target sequences
  do_cmd ("esl-alimanip --seq-k $ali.filter.names $ali.pfam.sto >> $protfile");

  #`rm $ali.*`;
}


sub do_cmd {
   my $cmd = $_[0];
   #print "$cmd\n"; 
   `$cmd`;
}

