library(statcomp)
library(ggplot2)
library(ggthemes)
library(scales)
library(reshape2)
library(rstudioapi)
library(cowplot)
library(gtools)
library(prodlim)
library(latex2exp)
library(pdc)
library(gridExtra)
library(tidyverse)
library(gtable)
library(latexpdf)
library(roxygen2)
library(devtools)
library(pkgload)
#library(StatOrdPattHxC)

# devtools::check("./StatOrdPattHxCModified/R/Test/")
# devtools::install("./StatOrdPattHxCModified/R/Test/")

setwd(dirname(getActiveDocumentContext()$path ))

files = list.files("./StatOrdPattHxCModified/R/")

for(i in files){
  source(paste("./StatOrdPattHxCModified/R/",i,sep=""))
}

files = list.files("./StatOrdPattHxC/R/")

for(i in files){
  source(paste("./StatOrdPattHxC/R/",i,sep=""))
}


dataPath = "./../../Data/CSV/weather.csv"

runEntropyTable = F
runConfidenceIntervalPlot = F
runNoise = F
runWhiteNoise = F
runPvalTheoretical = F
runPvalRandom = F

#EntropyTable Code
#####################################################################################
if(runEntropyTable){
  set.seed(123)
  df = data.frame()
  for(date in c(T,F)){ #date true means to set date to 14/08/1992, will only be done for noise==noNA==F
    for(noise in c(T,F)){
      for(noNA in c(T,F)){
        if(date & (noise | noNA)){
          next
        }
        data = read.csv(dataPath)
        data = select(data, -STATION)


        data["NAME"] <- replace(data["NAME"], data["NAME"]=="MIAMI INTERNATIONAL AIRPORT, FL US","Miami")
        data["NAME"] <- replace(data["NAME"], data["NAME"]=="EDINBURGH ROYAL BOTANIC GARDE, UK","Edinburgh")
        data["NAME"] <- replace(data["NAME"], data["NAME"]=="DUBLIN PHOENIX PARK, EI","Dublin")

        startDate = "1992-08-08"
        if(date){
          startDate = "19992-08-14"
          datesToDelete = c("1992-08-08","1992-08-09","1992-08-10","1992-08-11","1992-08-12","1992-08-13")
          data <- data[-which(data$DATE %in% datesToDelete),]
        }
        #remove na values
        if(noNA){
          data <- na.omit(data)
        }

        location = c("Miami","Edinburgh","Dublin")


        #for loop for counting number of patterns containing identical values for d=3
        #na.omit() needs to be activated above
        # count = c()
        # for (i in location){
        #   tmax = as.numeric(unlist(data[data["NAME"]==i,]["TMAX"]))
        #   n = length(tmax)
        #   counter=0
        #   for (j in 1:(n-2)){
        #     if(tmax[j]==tmax[j+1] | tmax[j]==tmax[j+2] | tmax[j+1]==tmax[j+2])
        #       counter = counter + 1
        #   }
        #   count = append(count,counter)
        # }
        # print(count)

        #add noise that only effect patterns with identical values. Both statcomp and pdc will always assign a observation containing
        #multiple identical values to the same pattern instead of a more evenly 50/50 divide between the two possible patterns.

        #commenting this line changes this entropy vastly especially for miami
        if(noise){
          data[,"TMAX"]=data[,"TMAX"]+runif(nrow(data),0,0.5)
        }


        d = 3
        #plot raw data
        #p <- ggplot(data=data,aes(DATE,TMAX,group=1))+geom_line()
        #p + facet_grid(rows=vars(NAME))


        #By using both na.omit() and adding tiny noise to all data points statcomp and pdc behaves in the same way, however if either
        #of the two preproccesing is lacking they start giving different entropies. The tiny noise is very important since 40% for miami and around 20-25%
        #Of the patterns for dublin and edinburgh have identical values. Statcomp and pdc assign these sets of observation to different patterns, which causes an error.
        #NA values are also handled differently between them


        #generate entropy tables
        d = 3
        entropyDf = data.frame()
        switch = c(1,4,2,5,3,6) #statcomp has a different ordinal pattern order then article implementation
        for (i in location){
          tmax = as.numeric(unlist((data[data["NAME"]==i,"TMAX"])))
          opd = ordinal_pattern_distribution(tmax,d) #statcomp
          n = sum(opd)
          opd = opd[switch]*(1/n)
          e = HShannon(opd)

          opd1 = codebook(tmax,m=d) #pdc
          e1 = HShannon(opd1)


          hist = histD3(tmax) #hist implementation of odp
          opd2 = c(0,0,0,0,0,0)
          for (j in unlist(hist)){
            opd2[j] = opd2[j]+1
          }
          n = sum(opd2)
          opd2 = opd2*(1/n) #opd article implementation
          e2= HShannon(opd2)

          e3= NA
          if(noNA){
            opd3 = OPprob(tmax,d) #StatOrdPattHxC newest implementation from article authors
            e3 = HShannon(opd3)
          }

          opd4 = formationPatternMagnus(tmax,D=d,tau=1) #my implementation with even split on ties
          n = sum(opd4)
          opd4 = opd4/n
          e4 = HShannon(opd4)



          entropyDf = rbind(entropyDf,c(e,e1,e2,e3,e4))
        }
        colnames(entropyDf) = c("statcomp","pdc","article implementation","StatOrdPattHxC","my implementation")
        entropyDf = round(entropyDf,digits=5)
        entropyDf$location <- location
        entropyDf$noiseAdded <- c(noise,"","")
        entropyDf$naOmitted <-c(noNA,"","")
        entropyDf$StartDate <- c(startDate,"","")
        entropyDf = relocate(entropyDf,noiseAdded,naOmitted,StartDate,location)
        df = rbind(df,entropyDf)
      }
    }
  }
  as.pdf(df,stem="entropyTable",dir="./../../Figures/PDFjpg/Weather/")
}
#####################################################################################

#Entropy plots where both noise and na.omit is activated. Confidence intervals around entropy
#####################################################################################
if(runConfidenceIntervalPlot){
  set.seed(123)

  data = read.csv(dataPath)
  data = select(data, -STATION)


  data["NAME"] <- replace(data["NAME"], data["NAME"]=="MIAMI INTERNATIONAL AIRPORT, FL US","Miami")
  data["NAME"] <- replace(data["NAME"], data["NAME"]=="EDINBURGH ROYAL BOTANIC GARDE, UK","Edinburgh")
  data["NAME"] <- replace(data["NAME"], data["NAME"]=="DUBLIN PHOENIX PARK, EI","Dublin")

  data <- na.omit(data)

  data[,"TMAX"]=data[,"TMAX"]+runif(nrow(data),0,0.5)

  location = c("Miami","Edinburgh","Dublin")

  df = data.frame()
  d=3
  switch = c(1,4,2,5,3,6) #statcomp has a different ordinal pattern order then article implementation
  for(i in location){
    tmax = as.numeric(unlist((data[data["NAME"]==i,"TMAX"])))
    opd = ordinal_pattern_distribution(tmax,ndemb=d)
    n = sum(opd)
    opd = opd[switch]*(1/n)
    entropyComplexity = global_complexity(opd=opd)[1:2]
    newVariance = sigma2q(tmax,d,"S")/n
    newSd = newVariance**(1/2)#/(n**(1/2)) #standard deviation
    z = qnorm(1-0.01/2) #z score
    newCI = z*newSd #One side of the confidence interval
    df = rbind(df,c(entropyComplexity,oldVariance,sigmanq,newVariance,2*oldCI,2*newCI))
  }
  df$location = location
  colnames(df) = c("entropy","complexity","oldVariance","sigmanq","newVariance","oldCILength","newCILength","location")

  print(df)
}


#comparison of noise vs theoretical evenly split of ties
#####################################################################################
if(runNoise){
  set.seed(123)

  iteration = 1000
  data = read.csv(dataPath)
  data = select(data, -STATION)


  data["NAME"] <- replace(data["NAME"], data["NAME"]=="MIAMI INTERNATIONAL AIRPORT, FL US","Miami")
  data["NAME"] <- replace(data["NAME"], data["NAME"]=="EDINBURGH ROYAL BOTANIC GARDE, UK","Edinburgh")
  data["NAME"] <- replace(data["NAME"], data["NAME"]=="DUBLIN PHOENIX PARK, EI","Dublin")

  data <- na.omit(data)

  df = data.frame()
  rows = nrow(data)
  location = c("Miami","Edinburgh","Dublin")
  d=3

  for (i in 1:iteration){
    rand = runif(rows,0,0.5)

    data[,"TMAX"]=data[,"TMAX"]+rand


    entropyList = c()
    for(i in location){
      tmax = as.numeric(unlist((data[data["NAME"]==i,"TMAX"])))
      e = global_complexity(tmax,ndemb=d)[1]
      entropyList = append(entropyList,e)
    }
    df = rbind(df,entropyList)
    data[,"TMAX"]=data[,"TMAX"]-rand
  }
  colnames(df) = location

  te = c() #Theoretiacal Entropies
  for(i in location){
    tmax = as.numeric(unlist((data[data["NAME"]==i,"TMAX"])))
    opd4 = formationPatternMagnus(tmax,D=d,tau=1) #my implementation with even split on ties
    n = sum(opd4)
    opd4 = opd4/n
    te = append(te,HShannon(opd4))
  }


  meanM = round(mean(df[,1]),digits=7)
  meanE = round(mean(df[,2]),digits=7)
  meanD = round(mean(df[,3]),digits=7)

  teM = round(te[1],digits=7) #miami
  teE = round(te[2],digits=7) #edinburgh
  teD= round(te[3],digits=7) # dublin

  diffM = round(abs(mean(df[,1])-te[1]),digits=7)
  diffE = round(abs(mean(df[,2])-te[2]),digits=7)
  diffD = round(abs(mean(df[,3])-te[3]),digits=7)

  printDf = data.frame(matrix(c(meanM,meanE,meanD,teM,teE,teD,diffM,diffE,diffD),nrow=3,ncol=3,byrow=F))
  colnames(printDf) = c("meanRandom","theoreticalSplit","diff")
  printDf$location = location
  printDf = relocate(printDf,location)
  as.pdf(printDf,stem="random_vs_theoreticalSplit",dir="./../../Figures/PDFjpg/Weather/")

  pdf("./../../Figures/PDFjpg/Weather/noiseStochasticTheoretical.pdf")
  par(mfrow=c(3,1))
  plot(x=(1:iteration),y=sort(df[,"Miami"]),type="l",xlab="iteration",ylab="Entropy",main="Miami")
  abline(h=te[1],v=iteration/2)

  plot(x=(1:iteration),y=sort(df[,"Edinburgh"]),type="l",xlab="iteration",ylab="Entropy",main="Edinburgh")
  abline(h=te[2],v=iteration/2)

  plot(x=(1:iteration),y=sort(df[,"Dublin"]),type="l",xlab="iteration",ylab="Entropy",main="Dublin")
  abline(h=te[3],v=iteration/2)
  dev.off()

}


#####################################################################################
if(runWhiteNoise){
  d=3
  set.seed(123)

  data = read.csv(dataPath)
  data = select(data, -STATION)


  data["NAME"] <- replace(data["NAME"], data["NAME"]=="MIAMI INTERNATIONAL AIRPORT, FL US","Miami")
  data["NAME"] <- replace(data["NAME"], data["NAME"]=="EDINBURGH ROYAL BOTANIC GARDE, UK","Edinburgh")
  data["NAME"] <- replace(data["NAME"], data["NAME"]=="DUBLIN PHOENIX PARK, EI","Dublin")

  data <- na.omit(data)

  iteration = 1000

  n=1000
  constant = rep(c(0),n)

  opd = formationPatternMagnus(constant,d,1,output=1)
  opd = opd/sum(opd)
  te = HShannon(opd)

  entropyList = c()
  for (i in 1:iteration){
    rand = runif(n,0,1)
    constant=constant+rand
    entropyList = append(entropyList,global_complexity(constant,ndemb=d)[1])
    constant=constant-rand
  }

  meanEntropy = mean(entropyList)
  diff = abs((te-meanEntropy))

  printDf = data.frame(matrix(c(meanEntropy,te,diff),nrow=1,ncol=3,byrow=F))
  colnames(printDf) = c("meanRandom","theoreticalSplit","diff")
  as.pdf(printDf,stem="random_vs_theoreticalSplitWhiteNoise",dir="./../../Figures/PDFjpg/Weather/")

  pdf("./../../Figures/PDFjpg/Weather/constantWithWhiteNoiseStochasticTheoretical.pdf")
  plot(x=(1:iteration),y=sort(entropyList),type="l",xlab="iteration",ylab="Entropy",main="ConstantValueAddedWhiteNoise")
  abline(h=te,v=iteration/2)
  dev.off()

}


################################################################################

if(runPvalTheoretical){
  d=3
  set.seed(123)
  
  data = read.csv(dataPath)
  data = select(data, -STATION)


  data["NAME"] <- replace(data["NAME"], data["NAME"]=="MIAMI INTERNATIONAL AIRPORT, FL US","Miami")
  data["NAME"] <- replace(data["NAME"], data["NAME"]=="EDINBURGH ROYAL BOTANIC GARDE, UK","Edinburgh")
  data["NAME"] <- replace(data["NAME"], data["NAME"]=="DUBLIN PHOENIX PARK, EI","Dublin")


  data <- na.omit(data)
  
  pvalDf = data.frame()

  miami = as.numeric(unlist((data[data["NAME"]=="Miami","TMAX"])))
  edinburgh = as.numeric(unlist((data[data["NAME"]=="Edinburgh","TMAX"])))
  dublin = as.numeric(unlist((data[data["NAME"]=="Dublin","TMAX"])))
  
  
  me = pval(miami,edinburgh,d,"S") #miami edinburgh
  md = pval(miami,dublin,d,"S") #miami dublin
  ed = pval(edinburgh,dublin,d,"S") #dublin edinburgh
  pvalDf = rbind(pvalDf,c(i,me,md,ed))
  data[,"TMAX"]=data[,"TMAX"]

  colnames(pvalDf) = c("Iteration","Miami-Edinburgh","Miami-Dublin","Edinburgh-Dublin")
  as.pdf(pvalDf,stem="pValuesTheoreticalTest",dir="./../../Figures/PDFjpg/Weather/")
}

if(runPvalRandom){
  d=3
  set.seed(123)
  iterations = 3
  
  data = read.csv(dataPath)
  data = select(data, -STATION)
  
  
  data["NAME"] <- replace(data["NAME"], data["NAME"]=="MIAMI INTERNATIONAL AIRPORT, FL US","Miami")
  data["NAME"] <- replace(data["NAME"], data["NAME"]=="EDINBURGH ROYAL BOTANIC GARDE, UK","Edinburgh")
  data["NAME"] <- replace(data["NAME"], data["NAME"]=="DUBLIN PHOENIX PARK, EI","Dublin")
  
  
  data <- na.omit(data)
  
  pvalDf = data.frame()
  for (i in 1:iterations){
    rand = runif(nrow(data),0,0.5)
    data[,"TMAX"]=data[,"TMAX"]+rand
    
    miami = as.numeric(unlist((data[data["NAME"]=="Miami","TMAX"])))
    edinburgh = as.numeric(unlist((data[data["NAME"]=="Edinburgh","TMAX"])))
    dublin = as.numeric(unlist((data[data["NAME"]=="Dublin","TMAX"])))
    
    
    me = pval(miami,edinburgh,d,"S") #miami edinburgh
    md = pval(miami,dublin,d,"S") #miami dublin
    ed = pval(edinburgh,dublin,d,"S") #dublin edinburgh
    pvalDf = rbind(pvalDf,c(i,me,md,ed))
    data[,"TMAX"]=data[,"TMAX"]-rand
  }
  colnames(pvalDf) = c("Iteration","Miami-Edinburgh","Miami-Dublin","Edinburgh-Dublin")
  sorted = pvalDf[order(pvalDf$`Edinburgh-Dublin`),]
  as.pdf(sorted,stem=paste("pValuesTheoretical,",iterations,"=Iterations,Sorted",sep=""),dir="./../../Figures/PDFjpg/Weather/")
}
