#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied ('fpw' or 'cons'), second is tblastn did not find hit ('no')", call.=FALSE)
} 


#else if (length(args)==1) {
  # default output file
#  args[2] = "out.txt"
#}
tooltechnique <- args[1]
tooltechnique <- sub("[\r\n]", "", tooltechnique)

tblastnfoundit <- args[2]
tblastnfoundit <- sub("[\r\n]", "", tblastnfoundit)

##print(sprintf("tooltechnique:%s\n",tooltechnique))

#print(sprintf("tooltechnique:%s",tooltechnique))
#print(sprintf("tblastnfoundit:%s",tblastnfoundit))


#print(sprintf("arg[1]:%s\n",args[1]))
#print(sprintf("arg[2]:%s\n",args[2]))


file_suffix <- paste(tooltechnique, ".pdf",sep="")

#pdf("HistogramHitsOnPositivesNotblastn.pdf",width=5,height=5)

myhmmdata <- read.table(file.path("logevaluediff"), header=T, sep="\t")
attach(myhmmdata)

if (identical(tblastnfoundit,"no")) {

  print("generating plot for phmmert and NO tblastn hit")
  title <- sprintf("log E-Values\n positive hits phmmert\n not found by tblastn %s", tooltechnique)
  output_file <- paste("HistogramHitsOnPositivesNotblastn", file_suffix,sep="")
  pdf(output_file,width=5,height=5)

  hist(logevaluediff,
     main=title,
#     main="Differences in log E-Values\n positive hits phmmert vs. tblastn",
     xlab="log E-Value phmmert",
     border="blue",
     col="green",
     xlim=c(-3.0e+1,3.0e+1),
     ylim=c(0,800),
     breaks=50,
     las=1)


} else {
  
  title <- sprintf("Differences in log E-Values\n positive hits phmmert vs. tblastn %s", tooltechnique)
  output_file <- paste("HistogramHitsOnPositivestblastn", file_suffix,sep="")
  pdf(output_file,width=5,height=5)

  hist(logevaluediff,
     main=title,
#     main="Differences in log E-Values\n positive hits phmmert vs. tblastn",
     xlab="log E-Value tblastn - log E-Value phmmert",
     border="blue",
     col="green",
     xlim=c(-3.0e+1,3.0e+1),
     ylim=c(0,1500),
     breaks=300,
     las=1)


}

#hist(logevaluediff,
#     main="Differences in log E-Values\n positive hits phmmert vs. tblastn",
#     xlab="log E-Value tblastn - log E-Value phmmert",
#     border="blue",
#     col="green",
#     xlim=c(-3.0e+1,3.0e+1),
#     ylim=c(0,800),
#     breaks=50,
#     las=1)

