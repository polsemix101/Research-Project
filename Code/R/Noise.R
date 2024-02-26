library(stats)
library(statcomp)
library(ggplot2)
library(ggthemes)

path <- "C:/Users/magnu/Desktop/skole/Vuw/BA/Research-Project/"

N <- c(10^3,10^4,10^5,10^6)
D <- c(3,4,5,6)
alpha <- c(0,0.5,1,1.5,2,2.5,3)
color <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2")

load(file=paste(path,"Data/Rmd/EdgesHxC.RData", sep=""))


#s is a boolean and should be true if raw data needs to be plotted - the function hasn't been tested for
#multiple n and d values while s is true. By setting S false the script runs faster.
generate_data <- function(N,D,alpha,path,rep, s){
  HC = NULL
  data <- NULL
  for (n in N){
    for (d in D){
      for(k in alpha){ 
        dataTemp <- NULL
        E <- NULL
        C <- NULL
        for (r in 1:rep){ 
          if (s){
            newData <- powernoise(k,n)[["x"]]
            dataTemp <- append(dataTemp,list(newData))  
          }

          EC <- unname(global_complexity(x=newData,ndemb=d)) #Entropy and Complexity
          E <- c(E,EC[1])
          C <- c(C,EC[2])
        }
        pca <- prcomp(data.frame(E,C), scale=T,center=T,retx=T)
        a <- pca$x[1:rep,1]
        median <- as.integer(rep/2+0.5)
        if (rep%%2==0){
          median <- c(as.integer(rep/2),as.integer(rep/2+1))
        }
        median <- which(a==sort(a)[median])
        bool <- rep(F,rep)
        bool[median] <- T

        HC <- rbind(HC,data.frame(k,n,d,E,C,bool))
        
        if (s){
          df <- data.frame(rep(k,n),1:n,dataTemp[median])
          names(df) <- c("k","t","Noise")
          data <- rbind(data,df)
        }

        #write to txt file not nessecary right now with small n, can be used for n=6 to generate sample data once.
        #df <- cbind(data.frame(E),data.frame(C))
        #write.table(df, paste(path,"Data/CSV/","n=",n,",D=",d,"k=",k,".txt",sep=""),row.names=F, sep=",")
      }
    }
  }
  HC <- data.frame(HC)
  names(HC) <- c("k", "Length", "Dimension", "Entropy", "Complexity","median")
  HC$k <- as.factor(HC$k)
  HC$Length <- as.factor(HC$Length)
  HC$Dimension <- as.factor(HC$Dimension)
  
  
  p1 <- ggplot(HC, aes(x=Entropy, y=Complexity, col=k)) +
    geom_point(size=1, alpha=.7) + 
    geom_line(data = subset(LinfLsup, Dimension==d & Side == "Lower"), 
      aes(x = V1, y = V2, col = NULL), alpha = 0.3) +
    geom_line(data = subset(LinfLsup, Dimension==d & Side == "Upper"), 
      aes(x = V1, y = V2, col = NULL), alpha = 0.3) +
    geom_line(data = subset(HC, median==TRUE), 
      aes(x = Entropy, y = Complexity, col = NULL), alpha = 0.3) +
    ggtitle(expression("301 points of"~italic(f)^{italic(-k)}~"noise")) +
    xlab(expression("Shannon Entropy"~italic(H))) +
    ylab(expression("Statistical Complexity  "~italic(C))) +
    scale_color_manual(values=color,
    name = expression(italic(k))) +
    guides(colour = guide_legend(override.aes = list(size=3, alpha=1))) +
    theme_tufte()
  
  pdf(paste(path,"Figures/PDFjpg/Noise/finalPlots.pdf",sep=""))
  plot(p1)
  
  if (s) {
    names(data) <- c("k","t","noiseData")
    p2 <- ggplot(data, aes(x=t,y=noiseData, col=k)) + geom_line() +
      facet_wrap(vars(k))
    plot(p2)  
  }
  dev.off()
  
    
  
 }

generate_data(10^3,6,alpha,path,11,T)




#this function was made for loading entropy-complexity data, so the generation part could be skipped.
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


#load_data(N,D,alpha,path)





