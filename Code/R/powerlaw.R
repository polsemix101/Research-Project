library(statcomp)
library(ggplot2)
library(rstudioapi)
library(tidyverse)
library(prodlim)
library(myLibrary)

setwd(dirname(getActiveDocumentContext()$path ))



set.seed(123)
k1=seq(1,5,by=2)
m1 = length(k1)

n=1000
iterations=10
d=3

pwHist = runif(n) #production of a simple histogram
opd = ordinal_pattern_time_series(pwHist,ndemb=d)[1:(n-d+1)]
pdf("./../../Figures/PDFjpg/powerlaw/histogram.pdf")
barplot(table(opd),xlab="pattern",ylab="frequency",main="Histogram of Ordinal Pattern distribution of random numbers",sub ="d=3 and n=1,000" )
dev.off()

pwTest = powernoise(5,n)[["x"]] #entropy of power law series with k=5
print(global_complexity(pwTest,ndemb=3)[1])

for(i3 in 1:m1){ #producing rejection and NaN plots for power law series
  rejectionList = c()
  NaNlist = c()
  tempK2 = exp(-5:0)
  k2 = append(append(rep(c(k1[i3]),4)-rev(tempK2),k1[i3]),rep(c(k1[i3]),4)+tempK2)
  m2=length(k2)
  for(i1 in 1:m2){
    rejectionCounter = 0
    NaNCounter = 0
    for (i2 in 1:iterations){
      pw1 = powernoise(k1[i3],n)
      pw2=powernoise(k2[i1],n)
      pvalue = pvalM2(pw1[["x"]],pw2[["x"]],d,"S")
      if (identical(pvalue,NaN) | pvalue<0.05){
        rejectionCounter = rejectionCounter+1/iterations
      }
      if(identical(pvalue,NaN)){
        NaNCounter=NaNCounter+1/iterations
      }
    }
    NaNlist = append(NaNlist,NaNCounter)
    rejectionList = append(rejectionList,rejectionCounter)
    print(i1)
  }
  pdf(paste("./../../Figures/PDFjpg/powerlaw/rejectionPlot,k1=",k1[i3],",n=",n,",iterations=",iterations,".pdf",sep=""))
  plot(k2[1:m2],rejectionList,xlab="k value of second power law series",ylab="Rejection frequency"
       ,pch=19,main="Testing of entropic hypothesis",
       sub=paste("k=",k1[i3]," for first power law series, n = ",n,", iterations = ",iterations,sep=""))
  dev.off()
  pdf(paste("./../../Figures/PDFjpg/powerlaw/NaNPlot,k1=",k1[i3],",n=",n,",iterations=",iterations,".pdf",sep=""))
  plot(k2[1:m2],NaNlist,xlab="k value of second power law series",ylab="NaN frequency"
       ,pch=19,main="Frequency of NaN pvalues in above testing",
       sub=paste("k = ",k1[i3]," for first power law series, n = ",n,", iterations = ",iterations,sep=""))
  dev.off()
}
