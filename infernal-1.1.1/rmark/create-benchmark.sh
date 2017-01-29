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




#phmmertpath=/home/um/wshands/pulloftranslatedsearch/hmmer/src
phmmertpath=/home/um/wshands/gitroot/hmmer/src
#export PATH=${phmmertpath}:$PATH
#export H_PATH=${phmmertpath}


transmarkpath=/home/um/wshands/TravisWheelerLabTransMark/transmark/infernal-1.1.1/

my_sub_path=/home/um/wshands/TravisWheelerLabTransMark/transmark/
#export PATH=${my_sub_path}:$PATH

if [ $# -gt 0 ]; then
    all_DNA_MSA_file=$1
else
    echo "Your command line contains no arguments, the first argument must be the file of all DNA MSAs"
    exit 1
fi

echo before comment
: <<'COMMENT'

#In some organisms only UGA is decoded as a stop codon, while UAG and UAA are 
#reassigned as sense codons. So what's happened is that we've asked easel to translate the DNA into proteins using the default codon table, but that table doesn't apply here. I haven't actually looked, but I'd bet money the the ORFs called by easel's translate code are all being stopped by stop codons "UAG" and "UAA" ... which aren't actually stop codons in Paramecium.

#How to overcome?  In this case, we could ask the translate code to use the ciliate code (you'd use "-c 6"). 

#But that's not really the solution to our problem. The thing is: we're building a benchmark "genome" with protein-coding sequences from all over the tree of life. That's not realistic, and it's getting us in trouble. Normally, when using phmmert, you'd know which kind of organism you were working with, so could just name the codon usage table at runtime with -c. But we can't do that here since each inserted sequence is coming from a different genome.  

#I think the solution is to restrict which sequences we put into the benchmark, ensuring that they all work with the standard translation table.  The way I'd do that is to take the DNA alignment file and run esl-translate on every sequence. Only keep sequences for which the translation finds a full-length ORF. After paring the list of sequences per alignment in this way, you could then go back, and run the benchmark-creation script on the alignment. 
#echo "Filtering the DNA MSAs so that only they only contain sequences with ORFs as long as the DNA sequence"⏎
#${transmarkpath}/../filter_sequences_with_same_size_ORF.pl all_filtered_ORF_alignments.stk $all_DNA_MSA_file⏎

#First get the names of all the alignments
echo "getting the names of the protein MSAs"
#esl-alistat all_filtered_ORF_alignments.stk  | awk '/Alignment name/ {split($3, basename, ".");  print basename[1]}' > all_ali_names.lst
esl-alistat all_filtered_ORF_alignments.stk  | awk '/Alignment name/ { print $3}' > all_ali_names.lst


#Just subsample 50% of all MSAs so there are not too many families and benchmark is smaller
#and get the names of the alignments to use
#(using "--random-source all_ali_names.lst" to set the seed)
echo "getting the names of only half of the protein MSAs"
sort -R --random-source all_ali_names.lst all_ali_names.lst | head -n 7362 > 7362_ali_names.lst

#now get the named  MSAs from the file with all the MSAs
echo "getting the protein MSAs"
#esl-afetch -f  all_filtered_ORF_alignments.stk 7362_ali_names.lst > 7362_alignments.stk
#also remove the '.stk.' suffix from the alignment names in the new 7362_alignments.stk file
#this is becuase the all_filtered_ORF_alignments.stk file erroneously has the '.stk' suffix at the end of the alignment names
esl-afetch -f  all_filtered_ORF_alignments.stk 7362_ali_names.lst | awk ' {where = match($0, /#=GF ID /); if (where !=0) {split($3, basename, "."); $3 =  basename[1]; print} else {print}}' > 7362_alignments.stk

echo "making the benchmark directory"
mkdir transmark_benchmark_data2
echo "cd'ing into the benchmark directory"
cd transmark_benchmark_data2

#put my_msub in path
#my_msub should be in repository

echo "generating the DNA background benchmark with decoy shuffled ORFs inserted into the background"
#remove the seed index file if it exists
rm -f  ../Pfam-A.v27.seed.ssi

${transmarkpath}/rmark/rmark-create  -N 10 -L 100000000 -R 10 -E 10 --maxtrain 30 --maxtest 20  -D ../Pfam-A.v27.seed transmarkORFandDNA ../7362_alignments.stk  ${transmarkpath}/rmark/rmark3-bg.hmm
#small test background sequence
#${transmarkpath}/rmark/rmark-create -X 0.2  -N 1 -L 100000000  -R 10 -E 10 --maxtrain 30 --maxtest 20  -D ../Pfam-A.v27.seed transmarkORFandDNA ../7362_alignments.stk  ${transmarkpath}/rmark/rmark3-bg.hmm


echo "downloading NCBI stand alone BLAST"
#curl -O ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
mkdir ncbi-blast
cd ncbi-blast
curl -O ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.6.0+-x64-linux.tar.gz
cd ..
$tblastn_path=$(pwd)/ncbi-blast/ncbi-blast-2.6.0+/bin/

echo "creating a DB for tblastn to use"
${tblastn_path}/makeblastdb -dbtype nucl -in transmarkORFandDNA.fa

echo "creating the file that has the amino acid MSAs that have the sequences will be used as query sequences against the background sequences"
${transmarkpath}/../build_protein_training_seeds.pl transmarkAminoAcidTest.msa transmarkORFandDNA.msa ../Pfam-A.v27.seed

echo "creating the HMMs to use as queries from the amino acid MSAs"
${phmmertpath}/hmmbuild transmarkAminoAcidTest.hmm transmarkAminoAcidTest.msa


echo "making the phmmert result directory"
mkdir ptr.std.e100
#do the searches 
#NOTE the result directories, e.g. tbn.w3.e100.fpw must be created before the command below is called
#otherwise the script will try to write to the directory possibly before it is created and you will loose
#result data
echo "running the phmmert search against the benchmark"
perl ${transmarkpath}/rmark/rmark-master.pl -F -N 16 -C transmarkAminoAcidTest.hmm  $H_PATH ${transmarkpath}/rmark ${transmarkpath}/rmark ptr.std.e100 ${transmarkpath}/rmark_opts/phmmert.e100.opts transmarkORFandDNA ${transmarkpath}/rmark/x-phmmert  1000000



#wait until the running jobs have finished (there is no output from qstat)
echo "Waiting for phmmert to finish; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done

COMMENT

tblastn_path=$(pwd)/ncbi-blast/ncbi-blast-2.6.0+/bin/


mkdir tbn.w3.e100.cons
perl ${transmarkpath}/rmark/rmark-master.pl -G transmarkAminoAcidTest  -F -N 16 $H_PATH ${transmarkpath}/rmark ${tblastn_path} tbn.w3.e100.cons ${transmarkpath}/rmark_opts/tblastn-w3-e100.opts transmarkORFandDNA ${transmarkpath}/rmark/x-tblastn-cons 1000000

#wait until the running jobs have finished (there is no output from qstat)
echo "Waiting for tblastn cons to finish; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done

mkdir tbn.w3.e100.fpw
perl ${transmarkpath}/rmark/rmark-master.pl -G transmarkAminoAcidTest  -F -N 16 $H_PATH ${transmarkpath}/rmark ${tblastn_path} tbn.w3.e100.fpw ${transmarkpath}/rmark_opts/tblastn-w3-e100.opts transmarkORFandDNA ${transmarkpath}/rmark/x-tblastn-fpw 1000000

#wait until the running jobs have finished (there is no output from qstat)
echo "Waiting for tblastn fpw to finish; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done

mkdir ptr.std.e100.cons
perl ${transmarkpath}/rmark/rmark-master.pl -F -N 16 -G transmarkAminoAcidTest  $H_PATH ${transmarkpath}/rmark ${transmarkpath}/rmark ptr.std.e100.cons ${transmarkpath}/rmark_opts/phmmert.e100.opts transmarkORFandDNA ${transmarkpath}/rmark/x-phmmert-cons  1000000

#wait until the running jobs have finished (there is no output from qstat)
echo "Waiting for phmmert cons to finish; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done


echo "downloading Exonerate"
#http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
#http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0.tar.gz
#http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0-x86_64.tar.gz
mkdir exonerate
cd exonerate
curl -O http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0-x86_64.tar.gz
tar -xvzf exonerate-2.2.0-x86_64.tar.gz
cd ..
exonerate_path=$(pwd)/exonerate/exonerate-2.2.0-x86_64/bin/


mkdir exonerate.fpw
perl ${transmarkpath}/rmark/rmark-master.pl -G transmarkAminoAcidTest  -F -N 16 $H_PATH ${transmarkpath}/rmark ${exonerate_path} exonerate.fpw ${transmarkpath}/rmark_opts/exonerate.opts transmarkORFandDNA ${transmarkpath}/rmark/x-exonerate-fpw  1000000


#wait until the running jobs have finished (there is no output from qstat)
echo "Waiting for exonerate fpw to finish; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done

mkdir exonerate.cons
perl ${transmarkpath}/rmark/rmark-master.pl -G transmarkAminoAcidTest  -F -N 16 $H_PATH ${transmarkpath}/rmark ${exonerate_path} exonerate.cons ${transmarkpath}/rmark_opts/exonerate.opts transmarkORFandDNA ${transmarkpath}/rmark/x-exonerate-fpw  1000000


#wait until the running jobs have finished (there is no output from qstat)
echo "Waiting for exonerate cons to finish; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done


#gather statistics for how many positive embedded squences were found by the search tools
my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA tbn.w3.e100.cons 1" 1

#wait until the running jobs have finished (there is no output from qstat)
echo "Gathering tblastn cons statistics; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done

my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA tbn.w3.e100.fpw 1" 1

#wait until the running jobs have finished (there is no output from qstat)
echo "Gathering tblastn fpw statistics; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done

my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA  ptr.std.e100 1" 1

#wait until the running jobs have finished (there is no output from qstat)
echo "Gathering phmmert statistics; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done

my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA  ptr.std.e100.cons 1" 1

#wait until the running jobs have finished (there is no output from qstat)
echo "Gathering phmmert cons statistics; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done


my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA exonerate.fpw 1" 1

#wait until the running jobs have finished (there is no output from qstat)
echo "Gathering exonerate statistics; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done

my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA tbn.w3.e100.fpw 1 .orf" 1

#wait until the running jobs have finished (there is no output from qstat)
echo "Gathering tblastn fpw ORF statistics; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done

my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA tbn.w3.e100.cons 1 .orf" 1

#wait until the running jobs have finished (there is no output from qstat)
echo "Gathering tblastn cons ORF  statistics; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done

my_msub gather "${transmarkpath}/rmark/rmark-pp.sh transmarkORFandDNA  ptr.std.e100 1 .orf" 1

#wait until the running jobs have finished (there is no output from qstat)
echo "Gathering phmmert ORF statistics; press [CTRL+C] to stop.."
while [[ $(qstat -u wshands) ]]
do
  sleep 1
done


