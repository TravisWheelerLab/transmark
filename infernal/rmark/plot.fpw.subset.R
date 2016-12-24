pdf("cons.subset.pdf",width=5,height=5)
mydata <- read.table(file.path("fpwout"), header=T, sep="\t")
attach(mydata)
plot(tblastncons, phmmert, main="Sharper E-values with phmmert than tblastn cons", cex=0.8, xlab="tblastn - E-value", ylab="phmmert - E-value", 
       log="xy", xlim=c(100,1e-5), ylim=c(100,1e-5))
segments( 100,100,1e-5,1e-5, lwd=1, col="pink")

