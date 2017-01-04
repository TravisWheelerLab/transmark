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
title <- sprintf("Positive sequence E-values\nphmmert vs tblastn %s", tooltechnique)
output_file <- paste("PlotHitsOnPositivesphmmertvstblastn", file_suffix,sep="")
pdf(output_file,width=5,height=5)

#pdf("HitsOnPositives.pdf",width=5,height=5)

myhmmdata <- read.table(file.path("hmmadjustedout"), header=T, sep="\t")
attach(myhmmdata)



if (identical(tooltechnique,"fpw")) {

plot(tblastnfpw, phmmert, main=title, 
      cex=0.8, xlab="tblastn - E-value", ylab="phmmert - E-value", col="blue",
      pch=1, log="xy", xlim=c(1e+2,1e-170), ylim=c(1e+2,1e-170))

} else {

plot(tblastncons, phmmert, main=title, 
      cex=0.8, xlab="tblastn - E-value", ylab="phmmert - E-value", col="blue",
      pch=1, log="xy", xlim=c(1e+2,1e-170), ylim=c(1e+2,1e-170))


}

segments(1e+2,1e+2, 1e-170, 1e-170, lwd=1, col="pink")

