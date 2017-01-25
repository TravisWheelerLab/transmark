#!/usr/bin/env perl
use strict;
#use warnings

my ($dnamsafilteredfile, $dnamsafile) = @ARGV;
 
if (not defined $dnamsafilteredfile) {
  die "Need DNA filtered  MSA file name as first input\n";
}

if (not defined $dnamsafile) {
  die "Need DNA MSA as second input\n";
}


unlink $dnamsafilteredfile;

#get the names of the alignments in the DNA multiple alignment file
#these are  training (query) sequences
my @names = `esl-alistat $dnamsafile | grep "Alignment name" | awk '{print \$3}'`;
chomp @names;


foreach my $ali (@names) {
  print "$ali\n";

  #get the DNA multiple alignment (MSA) file 
  do_cmd ("esl-afetch  $dnamsafile $ali > $ali.dna_orf_check.sto");
  #get the names of the sequences in the DNA multiple alignment
  do_cmd ("esl-alistat --list $ali.dna_seq_orf.names $ali.dna_orf_check.sto");
  # sometimes, the DNA names are formatted like "001|PUR5_SYMTH" (why!?).
  # clean this up:
  do_cmd (qq[awk 'BEGIN{FS="|"} {if ( \$2 == \'\\n\' )  print \$1; else print \$2}' $ali.dna_seq_orf.names > $ali.dna_seq_orf_check.names2]);

  #convert the MSA to a file of sequences
  do_cmd ("esl-reformat -o $ali.dna_seqs.fa fasta $ali.dna_orf_check.sto");
 
  print "Opening $ali.dna_seq_orf_check.names2 to read DNA sequence names\n";

#  my dnaseqnamefile = "$ali.dna_seq_orf_check.names2";
  open my $handle, '<', "$ali.dna_seq_orf_check.names2" or die "Cannot open $ali.dna_seq_orf_check.names2: $!";
  chomp(my @dnaseqnames = <$handle>);

  #get a sequence from the DNA MSA and make sure there is an ORF of the same length
  #if there is add the sequence name to a list of names sequences to keep
  #
  my @dna_seqs_to_keep;
  my $seq_len = 0;

  my $dna_seqs_to_keep = "dna_seqs_with_same_length_ORFs.txt";
  print "Open a file named $dna_seqs_to_keep; die if there's an error\n";
  open my $fh, '>', $dna_seqs_to_keep or die "Cannot open dna_seqs_with_same_length_ORFs.txt: $!";

  print "Loop through each DNA sequence in the MSA\n";
  foreach my $dnaseqname (@dnaseqnames) {
    print "checking sequence $dnaseqname\n";

    print "Indexing the MSA sequence file $ali.dna_seqs.fa\n";
    do_cmd ("esl-sfetch --index $ali.dna_seqs.fa");

    #get a sequence from the MSA and write it to a file
    print "getting the sequence from the MSA sequence file\n";
    do_cmd ("esl-sfetch $ali.dna_seqs.fa $dnaseqname > $ali.dna_seq_orf.fa");
    #get the length of the sequence
    my $seq_len = 0;
#    $seq_len = do_cmd (qq[esl-seqstat $ali.dna_seq_orf.fa | awk ' {seq_len_line = match(\$0, /Total # residues: /); if (seq_len_line > 0) {\$seq_len = \$4; print \$4}  } ' ]);
#awk ' {seq_len_line = match($0, /Total # residues:/); if (seq_len_line > 0) { print $4}  }'

    $seq_len = `esl-seqstat $ali.dna_seq_orf.fa | awk ' {seq_len_line = match(\$0, /Total # residues: /); if (seq_len_line > 0) { print \$4 } } '`;
    $seq_len =~ s/\r|\n//g;

#    my $dna_seq_orf = "dna_seq_orf.fa";
    print "The sequence length is $seq_len\n";
    #Try to find an ORF that is as long as the DNA sequence
    my $output = do_cmd ("esl-translate -l $seq_len $ali.dna_seq_orf.fa");
    print "output of esl-translate is $output\n";
    #if there is one then save the name of the sequence
    if ($output != "") {
      print "Found an ORF the same length\n";
      print $fh "$dnaseqname\n"; 
#      push @dna_seqs_to_keep, $dnaseqname;
    } else {
      print "Did not find an ORF of the same length\n";
    }
  }
  close $fh;

  #keep only the DNA sequences in the MSA that have an ORF as long as the sequence
  print "Saving only sequences with an equivalent length ORF to MSA and writing MSA to file $dnamsafilteredfile\n";
  do_cmd ("esl-alimanip --seq-k $dna_seqs_to_keep $ali.dna_orf_check.sto >> $dnamsafilteredfile");


  close($handle);

  #`rm $ali.*`;
}


sub do_cmd {
   my $cmd = $_[0];
   print "$cmd\n"; 
   `$cmd`;
}

