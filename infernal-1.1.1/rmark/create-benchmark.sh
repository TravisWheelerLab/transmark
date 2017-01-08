#!/bin/bash

#The all_alignments.stk file was created by Travis from Pfam DB
#and contains all the multiple sequence alignments in the Pfam DB
#to get the number of MSAs in the file
#grep "\/\/" ../all_alignments.stk   | wc
#  14724   14724   44172

#First get the names of all the alignments
esl-alistat all_alignments.stk  | awk '/Alignment name/ { print $3}' > all_ali_names.lst

#Just subsample 50% of all MSAs so there are not too many families and benchmark is smaller
#and get the names of the alignments to use
#(using "--random-source all_ali_names.lst" to set the seed)
sort -R --random-source all_ali_names.lst all_ali_names.lst | head -n 7362 > some_ali_names.lst

#now get the named  MSAs from the file with all the MSAs
esl-afetch -f  all_alignments.stk some_ali_names.lst > 7362_alignments.stk


mkdir transmark_benchmark_data
cd transmark_benchmark_data

#put my_msub in path
#my_msub should be in repository

#commands to generate test benchmark:
#generate the DNA background benchmark with decoy shuffled ORFs inserted into the background
/home/um/wshands/TravisWheelerLabTransMark/transmark/infernal-1.1.1/rmark/rmark-create  -N 10 -L 100000000 -R 10 -E 10 --maxtrain 30 --maxtest 20  -D ../Pfam-A.v27.seed rmarkinsertORFandDNA ../7362_alignments.stk  /home/um/wshands/TravisWheelerLabTransMark/transmark/infernal-1.1.1/rmark/rmark3-bg.hmm

#create a DB for tblastn to use
#used ncbi-blast-2.2.31+-x64-linux.tar.gz
makeblastdb -dbtype nucl -in rmarkinsertORFandDNA.fa

#create the file that has the amino acid MSAs that have the sequences will be used as query sequences
#against the background sequences
../build_protein_training_seeds.pl

#create the HMMs to use as queries from the amino acid MSAs
hmmbuild transmarkAminoAcidTest.hmm transmarkAminoAcidTest.msa


#do the searches 
#NOTE the result directories, e.g. tbn.w3.e100.fpw must be created before the command below is called
#otherwise the script will try to write to the directory possibly before it is created and you will loose
#result data
perl ../rmark/rmark-master.pl -F -N 16 -C transmarkAminoAcidTest  $H_PATH ../rmark ../rmark ptr.std.e100 ../rmark_opts/phmmert.e100.opts rmarkinsertORFandDNA ../rmark/x-phmmert  1000000

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

perl ../rmark/rmark-master.pl -G transmarkAminoAcidTest  -F -N 16 $H_PATH ../rmark ../rmark tbn.w3.e100.cons ../rmark_opts/tblastn-w3-e100.opts rmarkinsertORFandDNA ../rmark/x-tblastn-cons 1000000

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

perl ../rmark/rmark-master.pl -G transmarkAminoAcidTest  -F -N 16 $H_PATH ../rmark ../rmark tbn.w3.e100.fpw ../rmark_opts/tblastn-w3-e100.opts rmarkinsertORFandDNA ../rmark/x-tblastn-fpw 1000000

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

perl ../rmark/rmark-master.pl -G transmarkAminoAcidTest  -F -N 16 $H_PATH ../rmark ../rmark exonerate.fpw ../rmark_opts/exonerate.opts rmarkinsertORFandDNA ../rmark/x-exonerate-fpw  1000000


#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

#gather statistics for how many positive embedded squences were found by the search tools
my_msub gather "../rmark/rmark-pp.sh rmarkinsertORFandDNA tbn.w3.e100.cons 1" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

my_msub gather "../rmark/rmark-pp.sh rmarkinsertORFandDNA tbn.w3.e100.fpw 1" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

my_msub gather "../rmark/rmark-pp.sh rmarkinsertORFandDNA  ptr.std.e100 1" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

my_msub gather "../rmark/rmark-pp.sh pfamtransDNAmark exonerate.fpw 1" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

my_msub gather "../rmark/rmark-pp.sh rmarkinsertORFandDNA tbn.w3.e100.fpw 1 .orf" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

my_msub gather "../rmark/rmark-pp.sh rmarkinsertORFandDNA tbn.w3.e100.cons 1 .orf" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

my_msub gather "../rmark/rmark-pp.sh rmarkinsertORFandDNA  ptr.std.e100 1 .orf" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done


