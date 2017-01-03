
#create file with difference between phmmert E-Value exponents and tblastn fpw E-Value exponents
awk  'BEGIN { print "                               search\ttblastn\tphmmert\tlogevaluediff"} FNR>1 {if (($2 < 100) && ($4 < 100)) printf("%40s\t%1.2e\t%1.2e\t%1.2e\n", $1,$2,$4, (log($2)/log(10))-(log($4)/log(10)) ) }' hmmout | sort -g -b -t$'\t' -k 4,8  > hmmfpwlogevaluediff

cp hmmfpwlogevaluediff logevaluediff
Rscript --vanilla  ~/TravisWheelerLabTransMark/transmark/infernal-1.1.1/rmark/plot.POShitsHistogram.R fpw

#create file with difference between phmmert E-Value exponents and tblastn fpw E-Value exponents
awk  'BEGIN { print "                               search\ttblastn\tphmmert\tlogevaluediff"} FNR>1 {if (($3 < 100) && ($4 < 100)) printf("%40s\t%1.2e\t%1.2e\t%1.2e\n", $1,$3,$4, (log($3)/log(10))-(log($4)/log(10)) ) }' hmmout | sort -g -b -t$'\t' -k 4,8  > hmmconslogevaluediff

cp hmmconslogevaluediff logevaluediff
Rscript --vanilla  ~/TravisWheelerLabTransMark/transmark/infernal-1.1.1/rmark/plot.POShitsHistogram.R cons

#get results where hit for phmmert exists and hit for tblastn cons did not exist
awk  'BEGIN { print "                               search\ttblastn\tphmmert\tlogevaluediff"} FNR>1 {if (($3 == 100) && ($4 < 100)) printf("%40s\t%1.2e\t%1.2e\t%1.2e\n", $1,$3,$4, (log($3)/log(10))-(log($4)/log(10)) ) }' hmmout | sort -g -b -t$'\t' -k 4,8  > hmmbutnoconslogevaluediff

cp hmmbutnoconslogevaluediff logevaluediff
Rscript --vanilla  ~/TravisWheelerLabTransMark/transmark/infernal-1.1.1/rmark/plot.POShitsHistogram.R cons no


#get results where hit for phmmert exists and hit for tblastn fpw did not exist
awk  'BEGIN { print "                               search\ttblastn\tphmmert\tlogevaluediff"} FNR>1 {if (($2 == 100) && ($4 < 100)) printf("%40s\t%1.2e\t%1.2e\t%1.2e\n", $1,$2,$4, (log($2)/log(10))-(log($4)/log(10)) ) }' hmmout | sort -g -b -t$'\t' -k 4,8  > hmmbutnofpwlogevaluediff

cp hmmbutnofpwlogevaluediff logevaluediff
Rscript --vanilla  ~/TravisWheelerLabTransMark/transmark/infernal-1.1.1/rmark/plot.POShitsHistogram.R fpw no









