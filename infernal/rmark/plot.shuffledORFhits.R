pdf("HitsOnShuffledORFs.pdf",width=5,height=5)

myhmmdata <- read.table(file.path("hmmorfout"), header=T, sep="\t")
#myfpwdata <- read.table(file.path("fpwout"), header=T, sep="\t")

#mydata = rbind(myhmmdata, myfpwdata)
#attach(mydata)

attach(myhmmdata)
plot(tblastnfpw, phmmert, main="Shuffled ORF E-values phmmert vs tblastn",
       cex=0.8, xlab="tblastn - E-value", ylab="phmmert - E-value", col="blue",
        log="xy", xlim=c(1e+2,1e-5), ylim=c(1e+2,1e-5))

points(tblastncons, phmmert, col="blue", pch=1)

myfpwdata <- read.table(file.path("fpworfout"), header=T, sep="\t")
attach(myfpwdata)
points(tblastnfpw, phmmert, col="red", pch=2)

myconsdata <- read.table(file.path("consorfout"), header=T, sep="\t")
attach(myconsdata)
points(tblastncons, phmmert, col="orange", pch=3)


segments(1e+2,1e+2, 1e-5, 1e-5, lwd=1, col="pink")

legend(1e+1,1e-4,c("phmmert","tblastn fpw", "tblastn cons"),col=c("blue","red","orange"),pch=c(1,2,3))
