library(ggplot2)
library(ggthemes)
library(statcomp)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

source("AsymptoticEntropyVariancesComparison.R")

#################################################################################

# Time series

ts_gen_Linear <- function(a, n, eps){
  
  TS <- c(a, a-2)
  b <- mean(TS[(length(TS)-1):length(TS)])
  
  for (i in 1:n){
    TS <- c(TS, b)
    TS <- c(TS, mean(TS[(length(TS)-1):length(TS)]))
    TS <- c(TS, TS[length(TS)-1] + 0.5)
    TS <- c(TS, TS[length(TS)-1] - 0.5)
    TS <- c(TS, TS[length(TS)-1] + 0.5)
    TS <- c(TS, TS[length(TS)-1] - 0.5)
    TS <- c(TS, mean(TS[(length(TS)-1):length(TS)]))
    TS <- c(TS, TS[length(TS)-1] - 0.5)
    TS <- c(TS, TS[length(TS)] + (1:7)*eps)
    TS <- c(TS, TS[length(TS)-1] - 0.5)
    TS <- c(TS, TS[length(TS)] - (1:5)*eps)
    b <- mean(TS[(length(TS)-1):length(TS)])
  }
  
  return(TS)
}

ene <- 126
eneL <- 6

tse <- data.frame(t = 1:(eneL * 21 + 2), X = ts_gen_Linear(20, eneL, 0.01)/20)

ggplot(tse, aes(x = t, y = X)) +
  geom_point(color = "purple", size=0.25) +
  geom_line(color = "purple") +
  # ggtitle("Time series") +
  # theme_hc() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(family = "serif", size = 25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank())

## OP histogram

OPFr <- ordinal_pattern_distribution(tse$X, ndemb = 3)
OPFreq <- c(OPFr[1], OPFr[4], OPFr[2], OPFr[5], OPFr[3], OPFr[6])

barplot(OPFreq, col = "darkorchid1", border = "darkorchid4",
        names.arg = c(expression(pi[1]), expression(pi[2]), expression(pi[3]),
                      expression(pi[4]), expression(pi[5]), expression(pi[6])),
        space = 0, yaxt = "n")

## Asymptotic Distribution

k <- 6
pL <- (1:k) / (k*(k+1)/2)

load(file="../../Data/R/Output.Rdata")

Output.melt <- reshape2::melt(Output, 
                              measure.vars=5:10,
                              id.vars=1:4,
                              variable.name="Type",
                              value.name="Entropy")

dataLinear <- subset(Output.melt, 
                     subset = Output.melt$Model=="Linear" &
                       Output.melt$k=="6" & 
                       Output.melt$n=="60000" &
                       Output.melt$epsilon == "0",
                     select = 5:6)

dataLinear_HS <- subset(dataLinear,
                        subset = dataLinear$Type=="HS",
                        select = 2)

muLinear_HS <- -sum(pL*log(pL))/log(6)
sigmaLinear_HS <- sqrt(sigma.n.q(270000, 3, pL)) / (log(6)) 


ggplot(dataLinear_HS, aes(x = Entropy)) + 
  geom_histogram(aes(y =..density..),
                 bins=nclass.FD(dataLinear_HS$Entropy),
                 colour = "blue", 
                 fill = "darkorchid1") +  
  stat_function(fun = dnorm,  args = list(mean = muLinear_HS, sd = sigmaLinear_HS),
                colour = "red", size = 1.2) +
  xlab("Normalized Permutation Entropy") +
  ylab("Density") +
  scale_y_continuous(labels = scales::comma) +
  theme_hc() +
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15))
