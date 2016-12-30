pdf("FPWHitsOnPositives.pdf",width=5,height=5)

myhmmdata <- read.table(file.path("hmmout"), header=T, sep="\t")

#mydata = rbind(myhmmdata, myfpwdata)
#attach(mydata)

attach(myhmmdata)
plot(tblastnfpw, phmmert, main="Pos seq E-values phmmert vs tblastn fpw", 
      cex=0.8, xlab="tblastn - E-value", ylab="phmmert - E-value", col="blue",
      pch=1, log="xy", xlim=c(1e+2,1e-170), ylim=c(1e+2,1e-170))

#positive hit plots
#myhmmconsdata <- read.table(file.path("hmmout"), header=T, sep="\t")
#attach(myhmmconsdata)
#points(tblastncons, phmmert, col="blue", pch=1)

#myfpwdata <- read.table(file.path("fpwout"), header=T, sep="\t")
#attach(myfpwdata)
#points(tblastnfpw, phmmert, col="red", pch=2)

#myconsdata <- read.table(file.path("consout"), header=T, sep="\t")
#attach(myconsdata)
#points(tblastncons, phmmert, col="orange", pch=3)

segments(1e+2,1e+2, 1e-170, 1e-170, lwd=1, col="pink")

#legend(1e-100,1e-50,c("phmmert", "tblastn cons"),
#        col=c("blue","orange"),
#         pch=c(1,2))

