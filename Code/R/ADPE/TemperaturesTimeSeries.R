library(statcomp)
library(ggplot2)
library(ggthemes)
library(scales)
library(reshape2)
library(rstudioapi)
library(cowplot)
library(gtools)
library(prodlim)

setwd(dirname(getActiveDocumentContext()$path))
print(getwd())
source("Bandt-Pompe.R")

path = "StatOrdPattHxC/StatOrdPattHxC/R/"

files = list.files(paste(getwd(),path,sep=""),"\\.R$")

for (i in files){
  source(paste("C:/Users/magnu/Desktop/skole/Vuw/BA/Research-Project/Code/R/",path,i,sep=""))
}
source("C:/Users/magnu/Desktop/skole/Vuw/BA/Research-Project/Code/R/ADPE/AsymptoticEntropyVariancesComparison.R")

#####################################################################################

## Read, organize, and select meteo data

meteo = read.csv("./../../../Data/CSV/weatherBig.csv")



set.seed(123)

# Data max
#index 16292 is 08/08/1992
meteo_dublin_max = meteo[meteo$NAME == "DUBLIN PHOENIX PARK, EI",][16292:26297,"TMAX"]
meteo_miami_max = meteo[meteo$NAME == "MIAMI INTERNATIONAL AIRPORT, FL US",][16292:26297,"TMAX"]
meteo_edinburgh_max = meteo[meteo$NAME == "EDINBURGH ROYAL BOTANIC GARDE, UK",][11910:21915,"TMAX"]



meteo_dublin_max = na.omit(meteo_dublin_max)
meteo_edinburgh_max = na.omit(meteo_edinburgh_max)
meteo_miami_max = na.omit(meteo_miami_max)

path = "/StatOrdPattHxC/StatOrdPattHxC/R/"

files = list.files(paste(getwd(),path,sep=""),"\\.R$")

for (i in files){
  source(paste(getwd(),path,i,sep=""))
}
source("C:/Users/magnu/Desktop/skole/Vuw/BA/Research-Project/Code/R/ADPE/AsymptoticEntropyVariancesComparison.R")


dataPath = "./../../../Data/CSV/weather.csv"
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
# set.seed(123)
# data[,"TMAX"]=data[,"TMAX"]+runif(nrow(data),0,0.5)

statcompOrdinalPattern = c()

n2 = length(meteo_miami_max)
rand = runif(n2,0,0.5)
meteo_miami_max=meteo_miami_max+rand
data[data["NAME"]=="Miami","TMAX"]=data[data["NAME"]=="Miami","TMAX"]+rand

statcompOrdinalPattern = append(statcompOrdinalPattern,global_complexity(meteo_miami_max,ndemb=3)[1])


n1 = length(meteo_edinburgh_max)
rand = runif(n1,0,0.5)

meteo_edinburgh_max=meteo_edinburgh_max+rand
data[data["NAME"]=="Edinburgh","TMAX"]=data[data["NAME"]=="Edinburgh","TMAX"]+rand

statcompOrdinalPattern = append(statcompOrdinalPattern,global_complexity(meteo_edinburgh_max,ndemb=3)[1])

n = length(meteo_dublin_max)
rand = runif(n,0,0.5)
meteo_dublin_max=meteo_dublin_max+rand
data[data["NAME"]=="Dublin","TMAX"]=data[data["NAME"]=="Dublin","TMAX"]+rand

statcompOrdinalPattern = append(statcompOrdinalPattern,global_complexity(meteo_dublin_max,ndemb=3)[1])




names(statcompOrdinalPattern) = c("miami","edinburgh","dublin")
print(statcompOrdinalPattern)


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
  a = sigma2q(tmax,d,"S")
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


# temp = formationPatternMagnus(meteo_edinburgh_max,D=3,tau=1)
# n = sum(temp)
# opd = temp/n
# print(opd)
# entropy = global_complexity(opd=opd)[1]
# print(entropy)






histD3 <- function(series){

  p.patterns = formationPattern(series, D = 3, tau = 1, 0)
  n.symbols = dim(p.patterns)[1]
  symbol = matrix(c(0,1,2,
                    0,2,1,
                    1,0,2,
                    1,2,0,
                    2,0,1,
                    2,1,0), ncol = 3, byrow = TRUE)
  index.rep = array(0, n.symbols)

  for(i in 1:n.symbols){
    for(j in 1:6){
      if(all(p.patterns[i,] == symbol[j, ])){
        index.rep[i]=j
        break
      }
    }
  }

  index.rep = index.rep[1:n.symbols]
  index.rep = data.frame(i = index.rep)
  return(index.rep)
}

hist_dublin_max = histD3(meteo_dublin_max)
hist_edinburgh_max = histD3(meteo_edinburgh_max)
hist_miami_max = histD3(meteo_miami_max)


hist_entropy = c()
for (j in c(hist_dublin_max,hist_edinburgh_max,hist_miami_max)){
  a = c(0,0,0,0,0,0)
  for (i in unlist(j)){
    a[i] = a[i]+1
  }
  n = sum(a)
  a = a*(1/n)
  entropy = global_complexity(opd=a)[1]
  hist_entropy = append(hist_entropy,entropy)
}
names(hist_entropy) = c("dublin","edinburgh","miami")
print(hist_entropy)




