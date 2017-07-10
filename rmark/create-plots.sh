#stop the script if it attempts to use any unset variables 
set -o nounset
#stop the script if any command fails (returns non zero)
set -o errexit



#create the table of E-Values for the hits for Positive inserted sequences by the tools
#sort by phmmert E-Values (hmm)
plot_src_path=/home/um/wshands/gitroot/transmark/rmark

$plot_src_path/compute-glob-table.pl hmm > hmmout

echo "creating file with difference between phmmert E-Value exponents and tblastn fpw E-Value exponents"
awk  'BEGIN { print "                               search\ttblastnfpw\ttblastncons\tphmmert\tphmmertcons"} FNR>1 {printf("%40s\t%1.2e\t%1.2e\t%1.2e\t%1.2e\n", $1, $2 == -1 ? 100 : $2, $3 == -1 ? 100 : $3 , $4 == -1 ? 100 : $4 , $5 == -1 ? 100 : $5  ) }' hmmout | sort -g -b -t$'\t' -k 4,8  > hmmadjustedout

echo "plotting the E-Values for positive  hits for phmmert and tblastn fpw"
Rscript --vanilla  ${plot_src_path}/plot.POShits.R fpw
echo "plotting the E-Values for positive  hits for phmmert and tblastn cons"
Rscript --vanilla  ${plot_src_path}/plot.POShits.R cons
echo "plotting the E-Values for positive  hits for phmmert cons and tblastn cons"
Rscript --vanilla  ${plot_src_path}/plot.POShits.R phmmertcons


echo "creating the table of E-Values for the hits for inserted decoy shuffled ORF sequences by the tools sorting by phmmert E-Values (hmm)"
${plot_src_path}/compute-glob-table.pl hmm orf > hmmorfout

echo "creating file with difference between phmmert E-Value exponents and tblastn fpw E-Value exponents"
awk  'BEGIN { print "                               search\ttblastnfpw\ttblastncons\tphmmert"} FNR>1 {printf("%40s\t%1.2e\t%1.2e\t%1.2e\n", $1, $2 == -1 ? 100 : $2, $3 == -1 ? 100 : $3 , $4 == -1 ? 100 : $4 ) }' hmmorfout | sort -g -b -t$'\t' -k 4,8  > hmmorfadjustedout

echo "plotting the E-Values for shuffled ORF hits for phmmert and tblastn fpw"
Rscript --vanilla  ${plot_src_path}/plot.shuffledORFhits.R fpw
echo "plot the E-Values for shuffled ORF hits for phmmert and tblastn cons"
Rscript --vanilla  ${plot_src_path}/plot.shuffledORFhits.R cons



echo "creating file with difference between phmmert E-Value exponents and tblastn fpw E-Value exponents"
awk  'BEGIN { print "                               search\ttblastn\tphmmert\tlogevaluediff"} FNR>1 {if (($2 > -1) && ($4 > -1)) printf("%40s\t%1.2e\t%1.2e\t%1.2e\n", $1,$2,$4, (log($2)/log(10))-(log($4)/log(10)) ) }' hmmout | sort -g -b -t$'\t' -k 4,8  > hmmfpwlogevaluediff

cp hmmfpwlogevaluediff logevaluediff
echo "creating positive hits histogram for phmmert vs tblastn fpw"
Rscript --vanilla  ${plot_src_path}/plot.POShitsHistogram.R fpw

echo "creating file with difference between phmmert E-Value exponents and tblastn cons E-Value exponents"
awk  'BEGIN { print "                               search\ttblastn\tphmmert\tlogevaluediff"} FNR>1 {if (($3 > -1) && ($4 > -1)) printf("%40s\t%1.2e\t%1.2e\t%1.2e\n", $1,$3,$4, (log($3)/log(10))-(log($4)/log(10)) ) }' hmmout | sort -g -b -t$'\t' -k 4,8  > hmmconslogevaluediff

echo "creating positive hits histogram for phmmert vs tblastn cons"
cp hmmconslogevaluediff logevaluediff
Rscript --vanilla  ${plot_src_path}/plot.POShitsHistogram.R cons

echo "getting results where hit for phmmert exists and hit for tblastn cons did not exist"
#just use the log of the E-Value for phmmert, don't include log E-Value for tlbastn 
awk  'BEGIN { print "                               search\ttblastn\tphmmert\tlogevaluediff"} FNR>1 {if (($3 == -1) && ($4 > -1)) printf("%40s\t%1.2e\t%1.2e\t%1.2e\n", $1,$3,$4, log($4)/log(10) ) }' hmmout | sort -g -b -t$'\t' -k 4,8  > hmmbutnoconslogevaluediff

echo "creating positive hits histogram where hit for phmmert exists and hit for tblastn cons did not"
cp hmmbutnoconslogevaluediff logevaluediff
Rscript --vanilla  ${plot_src_path}/plot.POShitsHistogram.R cons no


echo "getting results where hit for phmmert exists and hit for tblastn fpw did not exist"
#just use the log of the E-Value for phmmert, don't include log E-Value for tlbastn 
awk  'BEGIN { print "                               search\ttblastn\tphmmert\tlogevaluediff"} FNR>1 {if (($2 == -1) && ($4 > -1)) printf("%40s\t%1.2e\t%1.2e\t%1.2e\n", $1,$2,$4, log($4)/log(10) ) }' hmmout | sort -g -b -t$'\t' -k 4,8  > hmmbutnofpwlogevaluediff

echo "creating positive hits histogram where hit for phmmert exists and hit for tblastn fpw did not"
cp hmmbutnofpwlogevaluediff logevaluediff
Rscript --vanilla  ${plot_src_path}/plot.POShitsHistogram.R fpw no

echo "creating sensitivity ROC plot"
perl ${plot_src_path}/transmark-Rroc.pl -R ${plot_src_path}/../listfiles/transmarkORFandDNA.ROC.list phmmert_ROC_plot.pdf 0 "Sensitivity ROC plot"







