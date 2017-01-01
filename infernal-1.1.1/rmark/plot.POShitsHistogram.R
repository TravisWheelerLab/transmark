pdf("HistogramHitsOnPositives.pdf",width=5,height=5)

myhmmdata <- read.table(file.path("logevaluediff"), header=T, sep="\t")

#mydata = rbind(myhmmdata, myfpwdata)
#attach(mydata)

attach(myhmmdata)
#plot(tblastncons, phmmert, main="Pos seq E-values phmmert vs tblastn cons", 
#      cex=0.8, xlab="tblastn - E-value", ylab="phmmert - E-value", col="blue",
#      pch=1, log="xy", xlim=c(1e+2,1e-170), ylim=c(1e+2,1e-170))

hist(logevaluediff,
     main="Differences in log E-Values\n positive hits phmmert vs. tblastn",
     xlab="log E-Value phmmert - log E-Value tblastn",
     border="blue",
     col="green",
     xlim=c(-3.0e+1,3.0e+1),
     ylim=c(0,250),
     breaks=300,
     las=1)



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

#segments(1e+2,1e+2, 1e-170, 1e-170, lwd=1, col="pink")

#legend(1e-100,1e-50,c("phmmert", "tblastn cons"),
#        col=c("blue","orange"),
#         pch=c(1,2))

