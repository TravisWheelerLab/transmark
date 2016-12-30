pdf("HitsOnPositives.pdf",width=5,height=5)

myhmmdata <- read.table(file.path("hmmout"), header=T, sep="\t")

#mydata = rbind(myhmmdata, myfpwdata)
#attach(mydata)

attach(myhmmdata)
plot(tblastnfpw, phmmert, main="Positive sequence E-values phmmert vs tblastn", 
      cex=0.8, xlab="tblastn - E-value", ylab="phmmert - E-value", col="blue",
      pch=1, log="xy", xlim=c(1e+2,1e-170), ylim=c(1e+2,1e-170))

#positive hit plots
myhmmconsdata <- read.table(file.path("hmmout"), header=T, sep="\t")
attach(myhmmconsdata)
points(tblastncons, phmmert, col="blue", pch=1)

myfpwdata <- read.table(file.path("fpwout"), header=T, sep="\t")
attach(myfpwdata)
points(tblastnfpw, phmmert, col="red", pch=2)

myconsdata <- read.table(file.path("consout"), header=T, sep="\t")
attach(myconsdata)
points(tblastncons, phmmert, col="orange", pch=3)

#shuffled ORF hit plots
#myhmmorfconsdata <- read.table(file.path("hmmorfout"), header=T, sep="\t")
#attach(myhmmorfconsdata)
#points(tblastncons, phmmert, col="sienna", pch=23)

#myhmmorffpwdata <- read.table(file.path("hmmorfout"), header=T, sep="\t")
#attach(myhmmorffpwdata)
#points(tblastnfpw, phmmert, col="sienna", pch=23)

#myfpworfdata <- read.table(file.path("fpworfout"), header=T, sep="\t")
#attach(myfpworfdata)
#points(tblastnfpw, phmmert, col="yellow", pch=17)

#myconsorfdata <- read.table(file.path("consorfout"), header=T, sep="\t")
#attach(myconsorfdata)
#points(tblastncons, phmmert, col="purple", pch=8)


segments(1e+2,1e+2, 1e-170, 1e-170, lwd=1, col="pink")

legend(1e-100,1e-50,c("phmmert", "tblastn fpw", "tblastn cons"),
        col=c("blue","red","orange"),
         pch=c(1,2,3))

#legend(1e-175,1e-150,c("phmmert", "tblastn fpw", "tblastn cons", "phmmert ORF", "tblastn fpw ORF", "tblastn cons ORF"),
#        col=c("blue","red","orange","sienna","yellow","purple"),
#        pch=c(1,2,3,23,17,8))
