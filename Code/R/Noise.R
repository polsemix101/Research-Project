library(stats)
library(statcomp)

path <- "C:/Users/magnu/Desktop/skole/Vuw/BA/Research-Project/"

N <- c(10^3,10^4,10^5,10^6)
D <- c(3,4,5,6)
alpha <- c(0,0.5,1,1.5,2,2.5,3)
color <- c("blue","brown","green","red","yellow","pink","cyan")


generate_data <- function(N,D,alpha,path){
  for (i in 1:4){
    n <- N[i] #length of data series
    
    for (j in 1:4){
      d <- D[j] #embedding dimension
      for(m in 1:7){ #alpha/k value
        k <- alpha[m]
        E <- c()
        C <- c()
        for (r in 1:301){ #increase r interval to 1:301, when program is done
          data <- powernoise(k,n)[["x"]]
          EC <- unname(global_complexity(x=data,ndemb=d)) #Entropy and Complexity
          E <- append(E, c(EC[1])) #entropy
          C <- append(C, c(EC[2])) #complexity
        }
        df <- cbind(data.frame(E),data.frame(C))
        write.table(df, paste(path,"Data/CSV/","n=",n,",D=",d,"k=",k,".txt",sep=""),row.names=F, sep=",")
      }
    }
  }  
}

plot_data <- function(N,D,alpha,color,path){
  par(mfrow=c(4,4))
  for (i in 1:4){
    n <- N[i]
    for (j in 1:4){
      d <- D[j]
      pdf(file=paste(path,"Figures/PDFjpg/Noise/n=",n,",D=",d,".pdf",sep=""))
      plot(0:1,0:1, type="n", main=paste("n=",n,", D=",d), xlab="Entropy",ylab="Complexity")
      grid(nx = NULL, ny = NULL,
           lty = 1, col = "gray", lwd = 2)
      for (m in 1:7){
        k <- alpha[m]
        data <- read.delim(paste(path,"Data/CSV/n=",n,",D=",d,"k=",k,".txt",sep=""),sep=",",header = T)
        points(E,C,pch=19, col=color[m])
      }
      dev.off()
    }
  }
}



load_data <- function(N,D,alpha,path){
  pdf(file=paste(path,"Figures/PDFjpg/Noise/allPlots.pdf",sep=""))
  for (i in 1:4){
    n <- N[i]
    for (j in 1:4){
      d <- D[j]
      plot(0:1,0:1, type="n", main=paste("n=",n,", D=",d), xlab="Entropy",ylab="Complexity")
      grid(nx = NULL, ny = NULL,
           lty = 1, col = "gray", lwd = 2)
      x <- c()
      y <- c()
      for (m in 1:7){
        k <- alpha[m]
        data <- read.delim(paste(path,"Data/CSV/n=",n,",D=",d,"k=",k,".txt",sep=""),sep=",",header = T)
        pca <- prcomp(data, scale=T,center=T,retx=T)
        a <- pca$x[1:301,1]
        median <- which(a==sort(a)[151])
        points(data[1:301,1],data[1:301,2],pch=19, col=color[m])
        x <- append(x,c(data[median,1]))
        y <- append(y,c(data[median,2]))
      }
      lines(x,y,type="l")
    }
  }
  dev.off()
}


load_data(N,D,alpha,path)





