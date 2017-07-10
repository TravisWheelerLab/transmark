#/bin/bash
#stop the script if it attempts to use any unset variables 
set -o nounset
#stop the script if any command fails (returns non zero)
set -o errexit

#debug the script
set -x


if [ $# -lt 4 ]; then
    echo "ERROR: Your command line contains less than four arguments,"
    exit 1
fi
#set parameter 5 to the empty string if it is unset or null
orf_suffix=${5:-''}


#if [ -z $5 ]; then
#    orf_suffix=$5
#fi

sort_option='-g'
if [ $4 == "score" ]; then
    echo "using score instead of E-value"
    #when using score sort so that largest score
    #is at top of file
    sort_option='-rg'
fi


scriptdir=`dirname $0`

# get running time
ls $2/*.time   | perl $scriptdir/rmark-time.pl > $2/$2.time
# get MER
cat $2/*out${orf_suffix} | perl $scriptdir/rmark-multiply-evalues.pl $3 | sort ${sort_option} | perl $scriptdir/rmark-mer.pl $1.pos${orf_suffix} $2/$2.time >  $2/$2.mer${orf_suffix} #$2/$2.em$3.mer

#if not processing shuffled ORFs then create xy file
if [ -z $orf_suffix  ]; then
    # get ROC
    cat $2/*out | perl $scriptdir/rmark-multiply-evalues.pl $3 | sort ${sort_option} | $scriptdir/rmark-rocplot -N 10000 --seed 181 $1 - > $2/$2.xy  #$2/$2.em$3.xy
fi

# get mer from rmark-rocplot
#cat $2/*out${orf_suffix} | perl $scriptdir/rmark-multiply-evalues.pl $3 | sort ${sort_option} | $scriptdir/rmark-rocplot -N 10000 --mer --seed 181 $1 - > $2/$2.bmer${orf_suffix} # $2/$2.em$3.bmer
# get numbers of false negatives and false positives at E-threshold of 0.1 from rmark-rocplot (after E-value inflation)
#ARG=`echo 0.1/$3|bc -l`
#echo 0.1/$3|bc -l
#cat $2/*out${orf_suffix} | perl $scriptdir/rmark-multiply-evalues.pl $3 | sort ${sort_option} | $scriptdir/rmark-rocplot -N 10000 --Ethresh 0.1 --seed 181 $1 - > $2/$2.bEthresh${orf_suffix} # $2/$2.em$3.bEthresh

cp $2/$2.mer${orf_suffix} ./
#cp $2/$2.bmer${orf_suffix} ./
cp $2/$2.time ./

if [ -z $orf_suffix  ]; then
    cp $2/$2.xy ./
fi

# summarize files to stdout
echo -n $2 | awk '{printf("%-70s  ", $0)}' > $2/$2.sum${orf_suffix} #$2/$2.em$3.sum
echo -n $3 | awk '{printf("%4s  ", $0)}' >> $2/$2.sum${orf_suffix} #$2/$2.em$3.sum
#cat $2/$2.bmer${orf_suffix}     | awk '{printf("MER:  %5s  %5s  %5s  %5s   ", $3, $7, $11, $15)}' >> $2/$2.sum${orf_suffix} 
#cat $2/$2.bEthresh${orf_suffix} | awk '{printf("ETHRESH:0.1  %5s  %5s   ", ${orf_suffix}, $9)}' >> $2/$2.sum${orf_suffix} 
grep ummary $2/$2.mer${orf_suffix} | awk '{print $7}' >> $2/$2.sum${orf_suffix} 
cp $2/$2.sum${orf_suffix} ./
cat $2.sum${orf_suffix}

