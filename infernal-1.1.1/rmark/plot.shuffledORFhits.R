#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied ('fpw' or 'cons')", call.=FALSE)
} 

tooltechnique <- args[1]
tooltechnique <- sub("[\r\n]", "", tooltechnique)

print(sprintf("tooltechnique:%s",tooltechnique))
#print(sprintf("arg[1]:%s\n",args[1]))


file_suffix <- paste(tooltechnique, ".pdf",sep="")
title <- sprintf("Shuffled ORF E-values phmmert vs tblastn %s", tooltechnique)
output_file <- paste("PlotHitsOnShuffledORFstblastn", file_suffix,sep="")
pdf(output_file,width=5,height=5)

#column_to_read <- paste("tblastn",tooltechnique,sep="")
#print(sprintf("column to read:%s",column_to_read))

myhmmdata <- read.table(file.path("hmmorfout"), header=T, sep="\t")
attach(myhmmdata)

if (identical(tooltechnique,"fpw")) {


plot(tblastnfpw, phmmert, main=title,
       cex=0.8, xlab="tblastn - E-value", ylab="phmmert - E-value", col="blue",
        log="xy", xlim=c(1e+2,1e-5), ylim=c(1e+2,1e-5))

} else {


plot(tblastncons, phmmert, main=title,
       cex=0.8, xlab="tblastn - E-value", ylab="phmmert - E-value", col="blue",
        log="xy", xlim=c(1e+2,1e-5), ylim=c(1e+2,1e-5))




}
segments(1e+2,1e+2, 1e-5, 1e-5, lwd=1, col="pink")



