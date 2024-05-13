library(statcomp)
library(ggplot2)
library(dplyr)
library(adpe)
library(pdc)
library(this.path)
library(tidyverse)
library(prodlim)

setwd(this.path::here())


path = "/StatOrdPattHxC/StatOrdPattHxC/R/"

files = list.files(paste(getwd(),path,sep=""),"\\.R$")

for (i in files){
  source(paste(getwd(),path,i,sep=""))
}
source("C:/Users/magnu/Desktop/skole/Vuw/BA/Research-Project/Code/R/ADPE/AsymptoticEntropyVariancesComparison.R")


dataPath = "./../../Data/CSV/weather.csv"
data = read.csv(dataPath)
data = select(data, -STATION)

data["NAME"] <- replace(data["NAME"], data["NAME"]=="MIAMI INTERNATIONAL AIRPORT, FL US","Miami")
data["NAME"] <- replace(data["NAME"], data["NAME"]=="EDINBURGH ROYAL BOTANIC GARDE, UK","Edinburgh")
data["NAME"] <- replace(data["NAME"], data["NAME"]=="DUBLIN PHOENIX PARK, EI","Dublin")

#remove na values
data <- na.omit(data)


location = list("Miami","Edinburgh","Dublin")


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
set.seed(123)
data[,"TMAX"]=data[,"TMAX"]+runif(nrow(data),0,0.5)


# n=nrow(data[data["NAME"]=="Dublin",])
# rand = runif(n,0,0.5)
# data[data["NAME"]=="Dublin","TMAX"]=data[data["NAME"]=="Dublin","TMAX"]+rand
# 
# 
# n1=nrow(data[data["NAME"]=="Edinburgh",])
# rand = runif(n1,0,0.5)
# data[data["NAME"]=="Edinburgh","TMAX"]=data[data["NAME"]=="Edinburgh","TMAX"]+rand
# 
# 
# n2=nrow(data[data["NAME"]=="Miami",])
# rand = runif(n2,0,0.5)
# data[data["NAME"]=="Miami","TMAX"]=data[data["NAME"]=="Miami","TMAX"]+rand



d = 3 
#plot raw data
#p <- ggplot(data=data,aes(DATE,TMAX,group=1))+geom_line()
#p + facet_grid(rows=vars(NAME))




#By using both na.omit() and adding tiny noise to all data points statcomp and pdc behaves in the same way, however if either
#of the two preproccesing is lacking they start giving different entropies. The tiny noise is very important since 40% for miami and around 20-25%
#Of the patterns for dublin and edinburgh have identical values. Statcomp and pdc assign these sets of observation to different patterns, which causes an error. 
#NA values are also handled differently between them

df = data.frame()
d = 3

switch = c(1,4,2,5,3,6)
for (i in location){
  tmax = as.numeric(unlist((data[data["NAME"]==i,"TMAX"])))
  opd = ordinal_pattern_distribution(tmax,d)
  n = sum(opd)
  opd = opd[switch]*(1/n)

  opd1 = OPprob(tmax,d)

  entropyComplexity = global_complexity(opd=opd)[1:2]
  oldVariance = v.n.q(n,d,opd) #function from ADPE/AsymptoticEntropyVariancesComparison.R
  newVariance = sigma.n.q(n,d,opd) #function from ADPE/AsymptoticEntropyVariancesComparison.R
  oldSd = oldVariance**(1/2) #standard deviation
  newSd = newVariance**(1/2)#/(n**(1/2)) #standard deviation
  z = qnorm(1-0.01/2) #z score
  oldCI = z*oldSd #One side of the confidence interval
  newCI = z*newSd #One side of the confidence interval
  df = rbind(df,c(entropyComplexity,oldVariance,newVariance,2*oldCI,2*newCI))
}
df$location = location
colnames(df) = c("entropy","complexity","oldVariance","newVariance","oldCILength","newCILength","location")

print(df)

