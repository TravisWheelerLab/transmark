#!/usr/bin/env perl
use strict;

my $protfile = "rmarkinsertAA.msa";
my $dnafile = "rmarkinsertORFandDNA.msa";
my $seeds = "../Pfam-A.v27.seed";
unlink $protfile;

my @names = `esl-alistat $dnafile | grep "Alignment name" | awk '{print \$3}'`;
chomp @names;

foreach my $ali (@names) {
  print "$ali\n";
  do_cmd ("esl-afetch $seeds $ali > $ali.pfam.sto");
  do_cmd ("esl-alistat --list $ali.pfam.names $ali.pfam.sto");

  do_cmd ("esl-afetch  $dnafile $ali > $ali.dna.sto");
  do_cmd ("esl-alistat --list $ali.dna.names $ali.dna.sto");
  #sometimes, the DNA names are formatted like "001|PUR5_SYMTH" (why!?).
  # clean this up:
#  do_cmd (qq[awk 'BEGIN{FS="|"} {print \$2}' $ali.dna.names > $ali.dna.names2]);
  do_cmd (qq[awk 'BEGIN{FS="|"} {if ( \$2 == \'\\n\' )  print \$1; else print \$2}' $ali.dna.names > $ali.dna.names2]);

  do_cmd ("grep -f $ali.dna.names2 $ali.pfam.names > $ali.filter.names");

  do_cmd ("esl-alimanip --seq-k $ali.filter.names $ali.pfam.sto >> $protfile");

  #`rm $ali.*`;
}


sub do_cmd {
   my $cmd = $_[0];
   #print "$cmd\n"; 
   `$cmd`;
}

