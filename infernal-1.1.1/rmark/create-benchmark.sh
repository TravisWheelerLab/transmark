#!/bin/bash

#stop the script if it attempts to use any unset variables 
set -o nounset
#stop the script if any command fails (returns non zero)
set -o errexit

#The all_alignments.stk file was created by Travis from Pfam DB
#and contains all the multiple sequence alignments in the Pfam DB
#to get the number of MSAs in the file
#grep "\/\/" ../all_alignments.stk   | wc
#  14724   14724   44172

phmmertpath=/home/um/wshands/pulloftranslatedsearch/hmmer/src
export PATH=${phmmertpath}:$PATH

transmarkpath=/home/um/wshands/TravisWheelerLabTransMark/transmark/infernal-1.1.1/


#First get the names of all the alignments
echo "getting the names of the protein MSAs"
esl-alistat all_alignments.stk  | awk '/Alignment name/ {split($3, basename, ".");  print basename[1]}' > all_ali_names.lst


#Just subsample 50% of all MSAs so there are not too many families and benchmark is smaller
#and get the names of the alignments to use
#(using "--random-source all_ali_names.lst" to set the seed)
echo "getting the names of only half of the protein MSAs"
sort -R --random-source all_ali_names.lst all_ali_names.lst | head -n 7362 > 7362_ali_names.lst

#now get the named  MSAs from the file with all the MSAs
echo "getting the protein MSAs"
esl-afetch -f  all_alignments.stk 7362_ali_names.lst > 7362_alignments.stk

echo "making the benchmark directory"
mkdir transmark_benchmark_data
echo "cd'ing into the benchmark directory"
cd transmark_benchmark_data

#put my_msub in path
#my_msub should be in repository

echo "generating the DNA background benchmark with decoy shuffled ORFs inserted into the background"
#remove the seed index file if it exists
rm -f  ../Pfam-A.v27.seed.ssi

#${transmarkpath}/rmark/rmark-create  -N 10 -L 100000000 -R 10 -E 10 --maxtrain 30 --maxtest 20  -D ../Pfam-A.v27.seed transmarkORFandDNA ../7362_alignments.stk  ${transmarkpath}/rmark/rmark3-bg.hmm

#small test background sequence
${transmarkpath}/rmark/rmark-create -X 0.2  -N 1 -L 100000000  -R 10 -E 10 --maxtrain 30 --maxtest 20  -D ../Pfam-A.v27.seed transmarkORFandDNA ../7362_alignments.stk  ${transmarkpath}/rmark/rmark3-bg.hmm

echo "downloading NCBI stand alone BLAST"
#curl -O ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
mkdir ncbi-blast
cd ncbi-blast
curl -O ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.6.0+-x64-linux.tar.gz
export PATH=$(pwd)/ncbi-blast/ncbi-blast-2.6.0+/bin/:$PATH
cd ..

echo "creating a DB for tblastn to use"
makeblastdb -dbtype nucl -in transmarkORFandDNA.fa

echo "creating the file that has the amino acid MSAs that have the sequences will be used as query sequences against the background sequences"
${transmarkpath}/../build_protein_training_seeds.pl transmarkAminoAcidTest.msa transmarkORFandDNA.msa ../Pfam-A.v27.seed

echo "creating the HMMs to use as queries from the amino acid MSAs"
hmmbuild transmarkAminoAcidTest.hmm transmarkAminoAcidTest.msa

echo "making the phmmert result directory"
mkdir ptr.std.e100
#do the searches 
#NOTE the result directories, e.g. tbn.w3.e100.fpw must be created before the command below is called
#otherwise the script will try to write to the directory possibly before it is created and you will loose
#result data
echo "running the phmmert search against the benchmark"
perl ${transmarkpath}/rmark/rmark-master.pl -F -N 16 -C transmarkAminoAcidTest  $H_PATH ${transmarkpath}/rmark ${transmarkpath}/rmark ptr.std.e100 ${transmarkpath}/rmark_opts/phmmert.e100.opts transmarkORFandDNA ${transmarkpath}/rmark/x-phmmert  1000000

#wait until the running jobs have finished (there is no output from qstat)
while [[ $(qstat -u wshands) ]]
do
  echo "Waiting for phmmert to finish; press [CTRL+C] to stop.."
  sleep 1
done

mkdir tbn.w3.e100.cons
perl ${transmarkpath}/rmark/rmark-master.pl -G transmarkAminoAcidTest  -F -N 16 $H_PATH ${transmarkpath}/rmark ${transmarkpath}/rmark tbn.w3.e100.cons ${transmarkpath}/rmark_opts/tblastn-w3-e100.opts transmarkORFandDNA ${transmarkpath}/rmark/x-tblastn-cons 1000000

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

#end the script here for debugging
exit

mkdir tbn.w3.e100.fpw
perl ${transmarkpath}/rmark/rmark-master.pl -G transmarkAminoAcidTest  -F -N 16 $H_PATH ${transmarkpath}/rmark ${transmarkpath}/rmark tbn.w3.e100.fpw ${transmarkpath}/rmark_opts/tblastn-w3-e100.opts transmarkORFandDNA ${transmarkpath}/rmark/x-tblastn-fpw 1000000

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

mkdir exonerate.fpw
perl ${transmarkpath}/rmark/rmark-master.pl -G transmarkAminoAcidTest  -F -N 16 $H_PATH ${transmarkpath}/rmark ${transmarkpath}/rmark exonerate.fpw ${transmarkpath}/rmark_opts/exonerate.opts transmarkORFandDNA ${transmarkpath}/rmark/x-exonerate-fpw  1000000


#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

mkdir exonerate.cons
perl ${transmarkpath}/rmark/rmark-master.pl -G transmarkAminoAcidTest  -F -N 16 $H_PATH ${transmarkpath}/rmark ${transmarkpath}/rmark exonerate.cons ${transmarkpath}/rmark_opts/exonerate.opts transmarkORFandDNA ${transmarkpath}/rmark/x-exonerate-fpw  1000000


#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done


#gather statistics for how many positive embedded squences were found by the search tools
my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA tbn.w3.e100.cons 1" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA tbn.w3.e100.fpw 1" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA  ptr.std.e100 1" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA exonerate.fpw 1" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA tbn.w3.e100.fpw 1 .orf" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA tbn.w3.e100.cons 1 .orf" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done

my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA  ptr.std.e100 1 .orf" 1

#wait until the running jobs have finished (there is no output from qstat)
while [[ qstat -u wshands ]]
do
  echo "Press [CTRL+C] to stop.."
  sleep 1
done


