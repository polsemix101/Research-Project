library(statcomp)
library(ggplot2)
library(dplyr)
library(adpe)
library(pdc)
library(this.path)

setwd(this.path::here())

source("ADPE/AsymptoticEntropyVariancesComparison.R")
dataPath1 = "./../../Data/CSV/weatherBig.csv"


data1 = read.csv(dataPath1)
data1 = select(data1, -STATION)


data1["NAME"] <- replace(data1["NAME"], data1["NAME"]=="MIAMI INTERNATIONAL AIRPORT, FL US","Miami")
data1["NAME"] <- replace(data1["NAME"], data1["NAME"]=="EDINBURGH ROYAL BOTANIC GARDE, UK","Edinburgh")
data1["NAME"] <- replace(data1["NAME"], data1["NAME"]=="DUBLIN PHOENIX PARK, EI","Dublin")
data1 <- data1[data1["NAME"]!="SAO PAULO AEROPORT, BR",]

#remove na values
data1 <- na.omit(data1)



location = list("Miami","Edinburgh","Dublin")

#for loop for counting number of patterns containing identical values for d=3
# count = c()
# for (i in location){
#   tmax = as.numeric(unlist(data1[data1["NAME"]==i,]["TMAX"]))
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

for (i in 1:nrow(data1)){
  data1[i,"TMAX"]=data1[i,"TMAX"]+runif(1,0,0.9)
}

d = 3 
#plot raw data
#p <- ggplot(data1=data1,aes(DATE,TMAX,group=1))+geom_line()
#p + facet_grid(rows=vars(NAME))




#By using both na.omit() and adding tiny noise to all data points statcomp and pdc behaves in the same way, however if either
#of the two preproccesing is lacking they start giving different entropies. The tiny noise is very important since 40% for miami and around 20-25%
#Of the patterns for dublin and edinburgh have identical values. Statcomp and pdc assign these sets of observation to different patterns, which causes an error. 
#NA values are also handled differently between them

df = data.frame()

for (i in location){
  dataTemp = data1[data1["NAME"]==i,]
  tmax = as.numeric(unlist(dataTemp["TMAX"]))
  fji = c() #firstJanurayIndex,point of following loops is to find the first valid date of each year, which will make it easier to window
  endFound=F
  i = 1
  while(endFound==F){
    startDate = dataTemp[i,"DATE"]
    if (substr(startDate,6,10)!="01-01"){
      startDate = paste(as.integer(substr(startDate,1,4))+1,"-01-01",sep="")
      i = which(dataTemp$DATE >= startDate)[1]
    }
    if (startDate=="2019-01-01"){
      endFound = T
    }
    fji = append(fji,i)
    i = i+1
  }
  nYears = length(fji) #number of valid years
  years = 1 #window size
  step = 1 #step size
  s = 1 #startIndex
  e = s+years #endIndex
  entropies = c()
  while (e<=nYears){
    opd = ordinal_pattern_distribution(tmax[fji[s]:fji[e]]-1,d) #statcomp
    n=sum(opd)
    opd = opd*(1/n)
    entropies = append(entropies,unname(global_complexity(opd=opd)[1]))
    s = s + step
    e = e + step 
  }
  x = 1:(ceiling((nYears-years)/step))
  df = data.frame(x,entropies)
  names(df) <- c("x","y")
  p <- ggplot(df,aes(x=x,y=y))+ geom_point()
  print(p)
  entropyOfEntropy = global_complexity(entropies,ndemb=3)[1]
  print(entropyOfEntropy)
  
  # opd = ordinal_pattern_distribution(tmax,d) #statcomp
  # n=sum(opd)
  # opd = opd*(1/n)
  # entropyComplexity = global_complexity(opd=opd)[1:2]
  # oldVariance = v.n.q(n,d,opd) #function from ADPE/AsymptoticEntropyVariancesComparison.R
  # newVariance = sigma.n.q(n,d,opd) #function from ADPE/AsymptoticEntropyVariancesComparison.R
  # oldSd = oldVariance**(1/2) #standard deviation
  # newSd = newVariance**(1/2)#/(n**(1/2)) #standard deviation
  # z = qnorm(1-0.01/2) #z score
  # oldCI = z*oldSd #One side of the confidence interval
  # newCI = z*newSd #One side of the confidence interval
  # df = rbind(df,c(entropyComplexity,oldVariance,newVariance,2*oldCI,2*newCI))
}
# df$location = location
# colnames(df) = c("entropy","complexity","oldVariance","newVariance","oldCILength","newCILength","location")
# 
# print(df)
