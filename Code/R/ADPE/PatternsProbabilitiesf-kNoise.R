### This program generates and plots
### the probability function of ordinal patterns of embedding dimension 3
### for f^{-k} noise with varying k and varying time series size

library(ggplot2)
library(ggthemes)
library(fftw)
library(statcomp)
library(reshape2)
set.seed(seed = 1234567890, kind = "Mersenne-Twister")

series.generator.fk <- function(pp, y, n, k){
  Series = vector(mode="numeric")
  filtro = (1:n)^-(k/2)
  filtro = filtro / sum(filtro)
  y1 = y * filtro    
  x1 = IFFT(y1, plan = pp)  
  Series = c(Re(x1)) 
  return(Series)
}

N <- c(300, 3000, 30000, 300000, 3000000) # Time series lengths
K <- c(.5, 1., 1.5, 2., 2.5, 3) # The exponent of the f^{-k} power spectrum
Rep <- 1000 # Number of replicas

pi.vector <- NULL # 
  # Columns 1-6: mean of pi hat
  # Column 7: N
  # Column *: K

for(current.K in K){
  for(current.N in N) {
    phat.sum <- rep(0, 6) # Stores the sum of phat over the replications
    for(current.R in 1:Rep) {
      
      x <- rnorm(current.N)  # Base series  
      x <- x - mean(x) # Subtract the mean
      pp <- planFFT(current.N) # Start the FFT process
      y <- FFT(x, plan = pp) # Compute the FFT of x
      fk <- series.generator.fk(pp, y, current.N, k=current.K) # Inject correlation 
      
      phat.sum <- phat.sum +
        ordinal_pattern_distribution(fk, ndemb=3)/(current.N-2)
    }
    phat.sum <- phat.sum / Rep # The average of phat over the replications
    pi.vector <- rbind(pi.vector,
                       c(phat.sum, current.N, current.K))
  }
}

pi.vector <- data.frame(pi.vector)
names(pi.vector) <- c("pi1", "pi2", "pi3", "pi4", "pi5", "pi6", "n", "k")
pi.vector$n <- as.factor(pi.vector$n)
pi.vector$k <- as.factor(pi.vector$k)

summary(pi.vector)

pi.vector.molten <- melt(pi.vector, measure.vars = 1:6, variable.name = "Pattern", value.name = "Probability")

pis <- c(expression(pi[1]), 
         expression(pi[2]), 
         expression(pi[3]), 
         expression(pi[4]), 
         expression(pi[5]),
         expression(pi[6]))

ggplot(pi.vector.molten, aes(x=Pattern, y=Probability, col=n)) +
  geom_point(size=3) +
  xlab(expression(pi)) +
  scale_x_discrete(labels=pis) +
  geom_line(aes(group=n)) +
  facet_wrap(~k) +
  theme_pander() +
  theme(text=element_text(family="serif"))
