library(ggplot2)
library(ggthemes)
library(statcomp)
library(gridExtra)
library(latex2exp)

# Let's see the position of the permutations using ordinal_pattern_distribution()

pi1 <- c(1,2,3)
plot(pi1, type="o")
ordinal_pattern_distribution(pi1, ndemb = 3)

pi2 <- c(1,3,2)
plot(pi2, type="o")
ordinal_pattern_distribution(pi2, ndemb = 3)

pi3 <- c(2,1,3)
plot(pi3, type="o")
ordinal_pattern_distribution(pi3, ndemb = 3)

pi4 <- c(2,3,1)
plot(pi4, type="o")
ordinal_pattern_distribution(pi4, ndemb = 3)

pi5 <- c(3,1,2)
plot(pi5, type="o")
ordinal_pattern_distribution(pi5, ndemb = 3)

pi6 <- c(3,2,1)
plot(pi6, type="o")
ordinal_pattern_distribution(pi6, ndemb = 3)

# So, the order is: pi1, pi3, pi5, pi2, pi4, pi6

########################### GENERATION FUNCTIONS ####################################

# Model: One is one


ts_gen_OiO <- function(a, n, eps){

  left_tail <- c(a+0.5, a, a+0.75, a+1, a+0.25)
  
  middle <- rep(c(a+0.5, a, a+0.75, a+0.25), n/6)
  
  d <- (1:(n/3-1)) * rep(eps, n/3-1)
  right_tail <- c(a+0.375, a+0.375-d)
  
  TS <- c(left_tail, middle[-c(length(middle)-2, length(middle)-1, length(middle))], right_tail)
  return(TS)
}

# Model: Half and half


ts_gen_HaH <- function(a, n, eps){
  
  left_tail <- c(a, a + 0.5)
  for (i in 1:(n/6+eps)){
    left_tail <- c(left_tail, a + 1, a + 0.75, a + 1.25)
    a <- left_tail[length(left_tail)]
  }
  
  right_tail <- vector()
  b <- left_tail[length(left_tail)-1] - 0.25
  for (i in 1:(n/6-eps)){
    right_tail <- c(right_tail, b, b - 0.75, b - 0.5)
    b <- right_tail[length(right_tail)-1] - 0.25      
  }
  
  TS <- c(left_tail, right_tail)
  return(TS)
}

# Model: Linear

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


########################### EXAMPLES ####################################

# Time series

ene <- 126
eneL <- 6

ts1 <- data.frame(t = 1:(ene + 2), X = ts_gen_OiO(0, ene, 0.01))
ts2 <- data.frame(t = 1:(ene + 2), X = ts_gen_HaH(0, ene, 7)/30)
ts3 <- data.frame(t = 1:(eneL * 21 + 2), X = ts_gen_Linear(20, eneL, 0.01)/20)

# Plots

p1 <- ggplot(ts1, aes(x = t, y = X)) +
  geom_point(color = 22, size=0.25) +
  geom_line(color = 22) +
  xlab(expression(italic(t))) +
  ylab(expression(italic(x)[italic(t)])) +
  # ggtitle("Model OiO") +
  theme_hc() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(family = "serif", size = 25))

p2 <- ggplot(ts2, aes(x = t, y = X)) +
  geom_point(color = 22, size=0.25) +
  geom_line(color = 22) +
  xlab(expression(italic(t))) +
  ylab(expression(italic(x)[italic(t)])) +
  # ggtitle("Model HaH") +
  theme_hc() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(family = "serif", size = 25))

p3 <- ggplot(ts3, aes(x = t, y = X)) +
  geom_point(color = 22, size=0.25) +
  geom_line(color = 22) +
  xlab(expression(italic(t))) +
  ylab(expression(italic(x)[italic(t)])) +
  # ggtitle("Model Lin") +
  theme_hc() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(family = "serif", size = 25))

# Altogether

grid.arrange(p1, p2, p3, ncol=1)

# Vectos N and p

ordinal_pattern_distribution(ts1$X, ndemb = 3)
ordinal_pattern_distribution(ts1$X, ndemb = 3)/ene

ordinal_pattern_distribution(ts2$X, ndemb = 3)
ordinal_pattern_distribution(ts2$X, ndemb = 3)/ene

ordinal_pattern_distribution(ts3$X, ndemb = 3)
sort(ordinal_pattern_distribution(ts3$X, ndemb = 3)/(21 * 6))

