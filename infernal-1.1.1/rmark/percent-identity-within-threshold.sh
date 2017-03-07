#!/bin/bash

#stop the script if it attempts to use any unset variables 
set -o nounset
#stop the script if any command fails (returns non zero)
set -o errexit

threshold=$3

tblastn_path=../transmark_benchmark_data2/ncbi-blast/ncbi-blast-2.6.0+/bin

echo "creating a DB for input target sequence for tblastn to use"
${tblastn_path}/makeblastdb -dbtype nucl -in ${1}

#find the percent identity of the two sequences
#if it is greater than the threshold return it 
echo "using tblastx to calculate percent identity for target sequence ${2} and query sequence ${1}"
percent_identity=0
percent_identity="$(${tblastn_path}/tblastx -word_size 3 -evalue 1.0e-50 -db ${1} -query ${2} -outfmt '7 pident' | awk -v threshold="${threshold}" ' /^#/ {next};  $1 > threshold  { ret=$1; print ret; exit; }')"

echo "percent identity of target sequence "${2}" and query sequence "${1}" is "${percent_identity}""

if [[ "${percent_identity}" > 0 ]] 
then
   echo "returning 1" 
   exit 1
else
   echo "returning 0"
   exit 0
fi


