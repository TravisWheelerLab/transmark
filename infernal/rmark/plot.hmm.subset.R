pdf("hmm.subset.pdf",width=5,height=5)
mydata <- read.table(file.path("hmmout"), header=T, sep="\t")
attach(mydata)
plot(tblastncons, phmmert, main="Sharper E-values with phmmert than tblastn cons", cex=0.8, xlab="tblastn - single sequence", ylab="phmmert - single sequence", 
        log="xy", xlim=c(1e+2,1e-5), ylim=c(1e+2,1e-5))
segments(1e+2,1e+2, 1e-5, 1e-5, lwd=1, col="pink")

