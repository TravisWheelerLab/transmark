#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied ('fpw' or 'cons' or 'phmmertcons')", call.=FALSE)
} 

tooltechnique <- args[1]
tooltechnique <- sub("[\r\n]", "", tooltechnique)

#print(sprintf("tooltechnique:%s",tooltechnique))
#print(sprintf("arg[1]:%s\n",args[1]))


file_suffix <- paste(tooltechnique, ".pdf",sep="")
if ((length(args))==1 && ((identical(tooltechnique,"fpw")) || (identical(tooltechnique,"cons")))) {
    title <- sprintf("Positive sequence E-values\nphmmert vs tblastn %s", tooltechnique)
    output_file <- paste("PlotHitsOnPositivesphmmertvstblastn", file_suffix,sep="")

} else {
    title <- sprintf("Positive sequence E-values\nphmmert cons vs tblastn cons")
    output_file <- paste("PlotHitsOnPositivesphmmertconsvstblastn", file_suffix,sep="")

}

pdf(output_file,width=5,height=5)

#pdf("HitsOnPositives.pdf",width=5,height=5)

myhmmdata <- read.table(file.path("hmmadjustedout"), header=T, sep="\t")
attach(myhmmdata)



if (identical(tooltechnique,"fpw")) {

# Define the position of tick marks
v1 <- c(0,1.0e-25,1.0e-50,1.0e-75,1.0e-100,1.0e-125,1.0e-150,1.0e-175)

# Define the labels of tick marks
v2 <- c("0","1.0e-25","1.0e-50","1.0e-75","1.0e-100","1.0e-125","1.0e-150","1.0e-175")

# Plot the data
#plot(x,
#     y,
#     xaxt = "n")

# Add axis to the plot 
#axis(side = 1, 
#     at = v1, 
#     labels = v2,
#     tck=-.1,
#     tcl = -0.5,
#    cex.axis=1.05,
#     col.axis="blue",
#     font.axis=5)


plot(tblastnfpw, phmmert, main=title, 
      cex=0.8, xlab="tblastn fpw - E-value", ylab="phmmert - E-value", col="blue",
      pch=1, log="xy", xlim=c(1e+2,1e-170), ylim=c(1e+2,1e-170))


#plot(tblastnfpw, phmmert, main=title, 
#      cex=0.8, xlab="tblastn fpw - E-value", ylab="phmmert - E-value", col="blue",
#      pch=1, log="xy", xlim=c(1e+2,1e-175), ylim=c(1e+2,1e-175),  yaxp=c(1e+2, 1e-175, 5))

} else if (identical(tooltechnique,"cons")) {

plot(tblastncons, phmmert, main=title, 
      cex=0.8, xlab="tblastn cons - E-value", ylab="phmmert - E-value", col="blue",
      pch=1, log="xy", xlim=c(1e+2,1e-170), ylim=c(1e+2,1e-170))

} else {

plot(tblastncons, phmmertcons, main=title, 
      cex=0.8, xlab="tblastn cons - E-value", ylab="phmmert cons - E-value", col="blue",
      pch=1, log="xy", xlim=c(1e+2,1e-170), ylim=c(1e+2,1e-170))

}

segments(1e+2,1e+2, 1e-170, 1e-170, lwd=1, col="pink")

