library(tuneR)
library(pdc)
library(ggplot2)
library(ggthemes)

# Shannon Entropy function -------------------------------------------------------------------------------
shannon.entropy <- function(prob){
  entropy = prob * log(prob)
  entropy[is.nan(entropy)] = 0
  return(-sum(entropy))
}

shannon.entropy.normalized <- function(prob){
  entropy = (shannon.entropy(prob)/log(length(prob)))
  entropy[is.nan(entropy)] = 0
  return(entropy)
}

constant <- function(prob){
  k = (0.5)/length(prob)
  a1 = (0.5 + k) * log(0.5 + k)
  a2 = (length(prob) - 1) * k * log(k)
  a3 = (1 - 0.5) * log(length(prob))
  b = -1/(a1 + a2 + a3)
  return(b)
}

# Jensen Divergence function ------------------------------------------------------------------------------
jensen.divergence <- function(prob){
  cc = rep(1/length(prob),length(prob))
  s_p = shannon.entropy(prob)
  s_q = shannon.entropy(cc)
  s_pq = shannon.entropy((prob + cc)/2)
  divergence = sum(s_pq - (s_p/2) - (s_q/2))
  return(divergence)
}

Ccomplexity<-function(prob){
  cc = jensen.divergence(prob) * constant(prob) * shannon.entropy.normalized(prob)
  return(cc)
}


# plot(noise(kind="power", alpha = 2, duration = 1000))

set.seed(1234567890, kin="Mersenne-Twister")

# Load the boundaries in the HxC plane
load("../../Data/Rmd/EdgesHxC.RData")

Rep <- 301
Length <- c(1000, 10000, 100000, 1000000)
Dimension <- c(3, 4, 5, 6)
K <- c(0, .5, 1., 1.5, 2., 2.5, 3)

HC <- NULL
for(k in K){
  for(len in Length){
    for(d in Dimension){
      for(r in 1:Rep){

        x <-  noise(kind="power", alpha = k, duration = len)
        ophist <- codebook(x@left, m=d, t=1, use.fast=TRUE, normalized=TRUE)
        HC <- rbind(HC,
                    c(k, len, d,
                      shannon.entropy.normalized(ophist),
                      Ccomplexity(ophist)))
      }
    }
  }
}

names(HC) <- c("k", "Length", "Dimension", "Entropy", "Complexity")
HC <- data.frame(HC)
HC$k <- as.factor(HC$k)
HC$Length <- as.factor(HC$Length)
HC$Dimension <- as.factor(HC$Dimension)

## Activate if needed
save(file="../../Data/Rmd/HxCf-k.Rdata", HC)
load(file="../../Data/Rmd/HxCf-k.Rdata")

# HC.melt <- melt(HC, measure.vars = c("Entropy", "Complexity"))


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2")

ggplot(HC, aes(x=Entropy, y=Complexity, col=k)) +
         geom_point(size=1, alpha=.7) +
  geom_line(data = subset(LinfLsup, Side == "Lower"), 
            aes(x = V1, y = V2, col = NULL), alpha = 0.3) +
  geom_line(data = subset(LinfLsup, Side == "Upper"), 
            aes(x = V1, y = V2, col = NULL), alpha = 0.3) +
  xlim(range(HC$Entropy)) +
  # ylim(range(HC$Complexity)) +
  facet_grid(Length~Dimension) +
  ggtitle(expression("301 points of"~italic(f)^{italic(-k)}~"noise")) +
  xlab(expression("Shannon Entropy"~italic(H))) +
  ylab(expression("Statistical Complexity"~italic(C))) +
  scale_color_manual(values=cbbPalette, 
                     name = expression(italic(k))) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha=1))) +
  theme_tufte()

## Only D=6 and L = 1000
ggplot(subset(HC, Dimension=="6" & Length=="1000"), aes(x=Entropy, y=Complexity, col=k)) +
  geom_boxplot(aes(x=Entropy, col=k), orientation="y", notch=TRUE) +
  # geom_boxplot(aes(y=Complexity, col=k), orientation="x", notch=TRUE) +
  geom_point(size=1, alpha=.2) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension=="6"),
             aes(x = V1, y = V2, col = NULL), alpha = 0.3) +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension=="6"),
             aes(x = V1, y = V2, col = NULL), alpha = 0.3) +
  xlim(range(HC$Entropy)) +
  # ylim(range(HC$Complexity)) +
  labs(title=expression("301 points of each type of"~italic(f)^{italic(-k)}~"noise"),
       subtitle=expression(italic(D==6)~"and"~italic(L==1000))) +
  # ggtitle(expression("301 points of each type of"~italic(f)^{italic(-k)}~"noise")) +
  # ggsubtitle("") +
  xlab(expression("Shannon Entropy"~italic(H))) +
  ylab(expression("Statistical Complexity"~italic(C))) +
  scale_color_manual(values=cbbPalette, 
                     name = expression(italic(k))) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha=1))) +
  theme_tufte()


## Only D=6, L=1000 and k in 0, 0.5, 1
ggplot(subset(HC, Dimension=="6" & Length=="1000" & (k=="0" | k=="0.5" | k=="1")), aes(x=Entropy, y=Complexity, col=k)) +
  geom_boxplot(aes(x=Entropy, col=k), orientation="y", notch=TRUE) +
  # geom_boxplot(aes(y=Complexity, col=k), orientation="x", notch=TRUE) +
  geom_point(size=1, alpha=.8) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension=="6"),
            aes(x = V1, y = V2, col = NULL), alpha = 0.3) +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension=="6"),
            aes(x = V1, y = V2, col = NULL), alpha = 0.3) +
  xlim(range(subset(HC, Dimension=="6" & Length=="1000" & (k=="0" | k=="0.5" | k=="1"))$Entropy)) +
  ylim(range(subset(HC, Dimension=="6" & Length=="1000" & (k=="0" | k=="0.5" | k=="1"))$Complexity)) +
  labs(title=expression("301 points of each type of"~italic(f)^{italic(-k)}~"noise"),
       subtitle=expression(italic(D==6)~","~italic(L==1000))) +
  xlab(expression("Shannon Entropy"~italic(H))) +
  ylab(expression("Statistical Complexity"~italic(C))) +
  scale_color_manual(values=cbbPalette, 
                     name = expression(italic(k))) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha=1))) +
  theme_tufte()


## Only D=3, L=1000 and k in 0, 0.5
ggplot(subset(HC, Dimension=="3" & Length=="1000" & (k=="0" | k=="0.5" )), aes(x=Entropy, y=Complexity, col=k)) +
  geom_boxplot(aes(x=Entropy, col=k), orientation="y", notch=TRUE) +
  # geom_boxplot(aes(y=Complexity, col=k), orientation="x", notch=TRUE) +
  geom_point(size=1, alpha=.8) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension=="3"),
            aes(x = V1, y = V2, col = NULL), alpha = 0.3) +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension=="3"),
            aes(x = V1, y = V2, col = NULL), alpha = 0.3) +
  xlim(range(subset(HC, Dimension=="3" & Length=="1000" & (k=="0" | k=="0.5" ))$Entropy)) +
  ylim(range(subset(HC, Dimension=="3" & Length=="1000" & (k=="0" | k=="0.5" ))$Complexity)) +
  labs(title=expression("301 points of each type of"~italic(f)^{italic(-k)}~"noise"),
       subtitle=expression(italic(D==3)~","~italic(L==1000))) +
  xlab(expression("Shannon Entropy"~italic(H))) +
  ylab(expression("Statistical Complexity"~italic(C))) +
  scale_color_manual(values=cbbPalette, 
                     name = expression(italic(k))) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha=1))) +
  theme_tufte()

