library(pdc)
library(devtools)
library(statcomp)


#install_gitlab(repo="freryal/asymptotic-distribution-of-various-types-of-entropyunder-the-multinomial-law",host='https://gitlab.ecs.vuw.ac.nz/')

path <-"C:/Users/magnu/Desktop/skole/Vuw/BA/Research-Project/Data/TXT/Nonlinear_time_series_analysis_of_palaeoclimate_proxy_records/"

#file names
fn <- list("medisect","ODP662_sst","ODP722_sst","ODP659_dust","ODP659_isotope","ODP721_dust","ODP967_dust") 

data <- list()
for (i in 1:3){
  data <- append(data,list(read.delim(paste(path,fn[i],".txt",sep=""), sep="",skip=0)))
}
for (i in 1:4){
  data <- append(data,list(read.delim(paste(path,fn[i+3],".txt",sep=""), sep="",skip=1)))
}

data[[2]] <- (data[[2]])[-c(616,617,618), ] #Delete line 616-618 in "ODP662_sst"



se <-function(lst){ #shannon entropy
  entropy = 0
  for (i in 1:length(lst)){
    if (lst[i]!=0) {
      entropy <- entropy - lst[i]*log(lst[i])#/log(factorial(D))
    }
  }
  entropy
}

  
D=3
t=2

pdf(paste(path,"plots.pdf"), width=4, height=4)
for (i in 1:7){
  
  sData <- list() #splitted/windowed data
  tData <- list() #time for each datapoint
  
  #Find the correct columns for each dataset

  age <- colnames(data[[i]])[grepl("age|Age",colnames(data[[i]]))]
  if (length(age)>1){
    age <- age[3]
  }
  yName <- colnames(data[[i]])[grepl("d18O|flux|SST",colnames(data[[i]]))]
  if (length(yName)>1){
    yName <- yName[3]
  }
  
  #split data up into windows of size 410ka, with shift 41ka
  x <- data[[i]][,age]
  y <-(data[[i]])[,yName]
  nrow <- length(x) #number of rows
  min <- x[1] #min age
  max <- x[nrow] # max age
  c <- min #current age
  j=1 #start index for a given window
  k=1 #end index for a given window
  while (c+410 < max){
    while (x[k]<c+410){
      k<-k+1
    }
    k <- k-1
    
    sData <- append(sData,list(y[j:k]))
    tData <- append(tData,x[ceiling((j+k)/2)])
    
    while (x[j]<c+41){
      j<-j+1
    }
    c <- c+41
  }
  
  #calculate entropy
  nw <- length(tData)
  entropy <- list()
  for (n in 1:nw){
    f <- codebook(as.numeric(unlist(sData[n])),m=D,t=t)
    entropy <- append(entropy, se(f))
  }
  
  plot(as.numeric(tData),as.numeric(entropy),type="l", xlab='Age (ka BP)',ylab='S_order',main=fn[i])

}
dev.off()

# 