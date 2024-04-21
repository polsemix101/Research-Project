library(statcomp)
library(ggplot2)
library(dplyr)
library(adpe)
        
dataPath = "./../../Data/CSV/weather.csv"

data = read.csv(dataPath)
data = select(data, -STATION)

data["NAME"] <- replace(data["NAME"], data["NAME"]=="MIAMI INTERNATIONAL AIRPORT, FL US","Miami")
data["NAME"] <- replace(data["NAME"], data["NAME"]=="EDINBURGH ROYAL BOTANIC GARDE, UK","Edinburgh")
data["NAME"] <- replace(data["NAME"], data["NAME"]=="DUBLIN PHOENIX PARK, EI","Dublin")



#if the loop below is commented out patterns with NA values are simply skipped
#Set NA values equal to previous value, so data[i]=data[i-1]
# for (i in 1:nrow(data)){
#    if (is.na(data[i,"TMAX"])){
#      data[i,"TMAX"]= data[i-1,"TMAX"]
#    }
#   data[i,"TMAX"]=data[i,"TMAX"]+runif(1,0,0.9)
# }

#plot raw data
#p <- ggplot(data=data,aes(DATE,TMAX,group=1))+geom_line() 
#p + facet_grid(rows=vars(NAME))


d=3

location = list("Miami","Edinburgh","Dublin")
df = data.frame()

for (i in location){
  tmax = data[data["NAME"]==i,]["TMAX"]
  entropyComplexity = global_complexity(unlist(tmax),ndemb=d)[1:2]
  odp = ordinal_pattern_distribution(unlist(tmax),ndemb=d)
  n = sum(odp)
  oldVariance = v.n.q(n,d,odp*(1/n))
  newVariance = sigma.n.q(n,d,odp*(1/n))
  oldSd = oldVariance**(1/2) #standard deviation
  newSd = newVariance**(1/2)#/(n**(1/2)) #standard deviation
  zAlpha = qnorm(1-0.01/2) #z score
  oldCI = zAlpha*oldSd #half of the confidence interval
  newCI = zAlpha*newSd #half of the confidence interval
  df = rbind(df,c(entropyComplexity,oldVariance,newVariance,2*oldCI,2*newCI))
}
df$location = location
colnames(df) = c("entropy","complexity","oldVariance","newVariance","oldCILength","newCILength","location")

print(df)






