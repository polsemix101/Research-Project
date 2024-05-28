library(pdc)
library(devtools)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path ))

path <-"./../../data/TXT/pallcacocha_red_intensity.txt"
data <- read.delim(path, sep="\t", skip=103, nrow=15918)


nrow = nrow(data) #number of rows total



sData = list() #splitted/windowed data
tData = list() #time for each datapoint

D=5



crow = 0 #current row

while (crow+1000 <= nrow) {
  sData <- append(sData,list(data[crow:(crow+1000),3]))
  tData <- append(tData, list(as.double(data[crow:(crow+1000),2])))
  crow <- crow + 100
}


f = list() #frequency of patterns for each window. 
for (i in 1:length(sData)) {
  f <- append(f,list(codebook(as.numeric(unlist(sData[[i]])),m=D, t=1)))
}

se <-function(lst){ #shannon entropy
  entropy = 0
  for (i in 1:length(lst)){
    if (lst[i]!=0) {
      entropy <- entropy - lst[i]*log(lst[i])/log(factorial(D)) 
    }
  }
  entropy
} 

avg <-function(lst){ #average time data, to get average time for each window
  avg = 0.0
  for (i in 1:length(lst)){
    avg <- avg + as.double(lst[i])/1000.0
  avg
  }
}

median <- function(lst){#take median of time series, to be used as x-axis value
  lst[500]
}


entropy <- lapply(f, se)


time1 <- lapply(tData,median)



srow <- 11000

path = "./../../Figures/PDFjpg/ElNino/"
pdf(paste(path,"Entropy.pdf"))
plot(time1,entropy,pch=16, xlab='Time(cal. yr BP)',ylab='Normalized Shannon Entropy')
dev.off()
pdf(paste(path,"RedColor.pdf"))
plot(data[1:srow,2],data[1:srow,3], type='l', xlab='Time(cal. yr BP)', ylab='Red color intensity')
