#/bin/bash

scriptdir=`dirname $0`

# get running time
ls $2/*.time   | perl $scriptdir/rmark-time.pl > $2/$2.time
# get MER
cat $2/*out$4 | perl $scriptdir/rmark-multiply-evalues.pl $3 | sort -g | perl  $scriptdir/rmark-mer.pl $1.pos$4 $2/$2.time >  $2/$2.mer$4 #$2/$2.em$3.mer
# get ROC
#cat $2/*out$4 | perl $scriptdir/rmark-multiply-evalues.pl $3 | sort -g | $scriptdir/rmark-rocplot -N 10000 --seed 181 $1 - > $2/$2.xy$4 #$2/$2.em$3.xy
# get mer from rmark-rocplot
#cat $2/*out$4 | perl $scriptdir/rmark-multiply-evalues.pl $3 | sort -g | $scriptdir/rmark-rocplot -N 10000 --mer --seed 181 $1 - > $2/$2.bmer$4 # $2/$2.em$3.bmer
# get numbers of false negatives and false positives at E-threshold of 0.1 from rmark-rocplot (after E-value inflation)
#ARG=`echo 0.1/$3|bc -l`
#echo 0.1/$3|bc -l
#cat $2/*out$4 | perl $scriptdir/rmark-multiply-evalues.pl $3 | sort -g | $scriptdir/rmark-rocplot -N 10000 --Ethresh 0.1 --seed 181 $1 - > $2/$2.bEthresh$4 # $2/$2.em$3.bEthresh

cp $2/$2.mer$4 ./
#cp $2/$2.bmer$4 ./
cp $2/$2.time ./
#cp $2/$2.xy$4 ./


# summarize files to stdout
echo -n $2 | awk '{printf("%-70s  ", $0)}' > $2/$2.sum$4 #$2/$2.em$3.sum
echo -n $3 | awk '{printf("%4s  ", $0)}' >> $2/$2.sum$4 #$2/$2.em$3.sum
#cat $2/$2.bmer$4     | awk '{printf("MER:  %5s  %5s  %5s  %5s   ", $3, $7, $11, $15)}' >> $2/$2.sum$4 
#cat $2/$2.bEthresh$4 | awk '{printf("ETHRESH:0.1  %5s  %5s   ", $5, $9)}' >> $2/$2.sum$4 
grep ummary $2/$2.mer$4 | awk '{print $7}' >> $2/$2.sum$4 
cp $2/$2.sum$4 ./
cat $2.sum$4

