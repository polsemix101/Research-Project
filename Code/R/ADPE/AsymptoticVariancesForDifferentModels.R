library(rstudioapi)
library(tidyverse)
library(ggthemes)
library(latex2exp)

setwd(dirname(getActiveDocumentContext()$path))

source("AsymptoticEntropyVariancesComparison.R")


k <- 6
nfactor <- c(100, seq(500, 10000, 500))

Output <- NULL

for(nf in nfactor){
  
  # Compute p
  
  pOiO <- c(1, nf/3-1, rep(nf/6,4)) / nf
  pHaH <- c(rep(nf/6 + nf/12, 3), rep(nf/6 - nf/12, 3)) / nf
  pLin <- (1:k) / 21 

  # Compute asymptotic variances
  
  Output <- rbind(Output,
                    cbind(
                      rep(c("OiO", "HaH", "Lin"), 2),
                      rep(nf, 6),
                      c(rep("sigma", 3), rep("nu", 3)),
                      c(sigma.n.q(nf, 3, pOiO), sigma.n.q(nf, 3, pHaH), sigma.n.q(nf, 3, pLin),
                      v.n.q(nf, 3, pOiO), v.n.q(nf, 3, pHaH), v.n.q(nf, 3, pLin))
                      
    ))
}

Output <- data.frame(Output)
names(Output) <- c("Model", "n", "OPdist", "Var")
Output$Model <- as.factor(Output$Model)
Output$n <- as.numeric(Output$n)
Output$OPdist <- as.factor(Output$OPdist)
Output$Var <- as.numeric(Output$Var)

# Plot altogether

neworder <- c("OiO","HaH","Lin")
Output <- arrange(transform(Output,
                           Model=factor(Model,levels=neworder)),Model)

ggplot(Output, aes(x = n, y = sqrt(Var), color = OPdist)) +
  facet_wrap(~Model, ncol = 1) +
  geom_point() +
  geom_line() +
  ylab("") +
  scale_color_manual(name = "Asymptotic standard deviation", values = c(21, 22), 
                     labels = c(expression(nu), expression(sigma))) +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 10))

# Differences

OutOiO <- Output %>% subset(Output$Model == "OiO")
OutOiOs <- OutOiO %>% subset(OutOiO$OPdist == "sigma")
OutOiOn <- OutOiO %>% subset(OutOiO$OPdist == "nu")
pos <- which(OutOiOs$Var - OutOiOn$Var < 10^(-3))
OutOiOs$n[pos]

OutHaH <- Output %>% subset(Output$Model == "HaH")
OutHaHs <- OutHaH %>% subset(OutHaH$OPdist == "sigma")
OutHaHn <- OutHaH %>% subset(OutHaH$OPdist == "nu")
pos <- which(OutHaHs$Var - OutHaHn$Var < 10^(-3))
OutHaHs$n[pos]

OutLin <- Output %>% subset(Output$Model == "Lin")
OutLins <- OutLin %>% subset(OutLin$OPdist == "sigma")
OutLinn <- OutLin %>% subset(OutLin$OPdist == "nu")
pos <- which(OutLins$Var - OutLinn$Var < 10^(-3))
OutLins$n[pos]

# Separated plots

ggplot(OutOiO, aes(x = n, y = sqrt(Var), color = OPdist)) +
  geom_point() +
  geom_line() +
  xlab(expression(italic(n))) +
  ylab("") +
  scale_color_manual(name = "Asymptotic standard deviation", values = c(21, 22),
                     breaks = c("sigma", "nu"),
                     labels = c(expression(sigma), expression(nu))) +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 25))

ggplot(OutHaH, aes(x = n, y = sqrt(Var), color = OPdist)) +
  geom_point() +
  geom_line() +
  xlab(expression(italic(n))) +
  ylab("") +
  scale_color_manual(name = "Asymptotic standard deviation", values = c(21, 22),
                     breaks = c("sigma", "nu"),
                     labels = c(expression(sigma), expression(nu))) +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 25))

ggplot(OutLin, aes(x = n, y = sqrt(Var), color = OPdist)) +
  geom_point() +
  geom_line() +
  xlab(expression(italic(n))) +
  ylab("") +
  scale_color_manual(name = "Asymptotic standard deviation", values = c(21, 22),
                     breaks = c("sigma", "nu"),
                     labels = c(expression(sigma), expression(nu))) +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 25))

####################################################################

# Comparing length of time series

new_n_OiO <- round((OutOiOs$n * OutOiOn$Var) / OutOiOs$Var)
new_n_HaH <- round((OutHaHs$n * OutHaHn$Var) / OutHaHs$Var)
new_n_Lin <- round((OutLins$n * OutLinn$Var) / OutLins$Var)

df_enes <- data.frame(n = OutOiOs$n,
                      nOiO = new_n_OiO,
                      nHaH = new_n_HaH,
                      nLin = new_n_Lin)

# Separated plots

ggplot(df_enes, aes(x = n, y = nOiO)) +
  geom_point(color = 22) +
  geom_line(color = 22) +
  geom_line(aes(x=n, y=n), col="gray") +
  geom_point(aes(x=n, y=n), col="gray") +
  coord_fixed() +
  xlab(expression(italic(n))) +
  ylab(TeX(r'($italic(n)\prime$)')) +
  ggtitle("One is One (OiO) Model") +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 20))
ggsave("../../Figures/PDF/new_n_OiO.pdf", width = 14, height=12, units="cm")

ggplot(df_enes, aes(x = n, y = nHaH)) +
  geom_point(color = 22) +
  geom_line(color = 22) +
  geom_line(aes(x=n, y=n), col="gray") +
  geom_point(aes(x=n, y=n), col="gray") +
  coord_fixed() +
  xlab(expression(italic(n))) +
  ylab(TeX(r'($italic(n)\prime$)')) +
  ggtitle("Half and Half (HaH) Model") +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 20))
ggsave("../../Figures/PDF/new_n_HaH.pdf", width = 14, height=12, units="cm")

ggplot(df_enes, aes(x = n, y = nLin)) +
  geom_point(color = 22) +
  geom_line(color = 22) +
  geom_line(aes(x=n, y=n), col="gray") +
  geom_point(aes(x=n, y=n), col="gray") +
  #scale_y_log10() +
  coord_fixed() +
  xlab(expression(italic(n))) +
  ylab(TeX(r'($italic(n)\prime$)')) +
  ggtitle("Linear (Lin) Model") +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 20))
ggsave("../../Figures/PDF/new_n_Lin.pdf", width = 14, height=12, units="cm")

####################################################################

# Colored noises

load("../../Data/R/pivector.Rdata")

View(pi.vector)

result <- NULL

for(ka in levels(pi.vector$k)){
  
  df <- pi.vector %>% subset(pi.vector$k == ka)

  for(ene in levels(pi.vector$n)){
  
    # Fix probability vector
  
    qu <- df %>% subset(df$n == ene, 1:6)
  
    # Compute asymptotic variances
  
    result <- rbind(result,
                  cbind(
                    rep(ka, 2),
                    rep(ene, 2),
                    c("sigma", "nu"),
                    c(sigma.n.q(as.numeric(ene), 3, as.vector(unlist(qu))),
                      v.n.q(as.numeric(ene), 3, as.vector(unlist(qu))))
                  ))
  }
}

result <- data.frame(result)
names(result) <- c("k", "n", "OPdist", "Var")
result$k <- as.factor(result$k)
result$n <- as.numeric(result$n)
result$OPdist <- as.factor(result$OPdist)
result$Var <- as.numeric(result$Var)

res05 <- result %>% subset(result$k == 0.5)
res1 <- result %>% subset(result$k == 1)
res15 <- result %>% subset(result$k == 1.5)
res2 <- result %>% subset(result$k == 2)
res25 <- result %>% subset(result$k == 2.5)
res3 <- result %>% subset(result$k == 3)

options(scipen = 999)

ggplot(res05, aes(x = as.factor(n), y = sqrt(Var), color = OPdist, group = OPdist)) +
  geom_point() +
  geom_line() +
  xlab(expression(italic(n))) +
  ylab("") +
  
  scale_color_manual(name = "Asymptotic standard deviation", values = c(21, 22),
                     breaks = c("sigma", "nu"),
                     labels = c(expression(sigma), expression(nu))) +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 25))

ggplot(res1, aes(x = as.factor(n), y = sqrt(Var), color = OPdist, group = OPdist)) +
  geom_point() +
  geom_line() +
  xlab(expression(italic(n))) +
  ylab("") +
  
  scale_color_manual(name = "Asymptotic standard deviation", values = c(21, 22),
                     breaks = c("sigma", "nu"),
                     labels = c(expression(sigma), expression(nu))) +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 25))

ggplot(res15, aes(x = as.factor(n), y = sqrt(Var), color = OPdist, group = OPdist)) +
  geom_point() +
  geom_line() +
  xlab(expression(italic(n))) +
  ylab("") +
  
  scale_color_manual(name = "Asymptotic standard deviation", values = c(21, 22),
                     breaks = c("sigma", "nu"),
                     labels = c(expression(sigma), expression(nu))) +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 25))

ggplot(res2, aes(x = as.factor(n), y = sqrt(Var), color = OPdist, group = OPdist)) +
  geom_point() +
  geom_line() +
  xlab(expression(italic(n))) +
  ylab("") +
  
  scale_color_manual(name = "Asymptotic standard deviation", values = c(21, 22),
                     breaks = c("sigma", "nu"),
                     labels = c(expression(sigma), expression(nu))) +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 25))

ggplot(res25, aes(x = as.factor(n), y = sqrt(Var), color = OPdist, group = OPdist)) +
  geom_point() +
  geom_line() +
  xlab(expression(italic(n))) +
  ylab("") +
  scale_color_manual(name = "Asymptotic standard deviation", values = c(21, 22),
                     breaks = c("sigma", "nu"),
                     labels = c(expression(sigma), expression(nu))) +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 25))

ggplot(res3, aes(x = as.factor(n), y = sqrt(Var), color = OPdist, group = OPdist)) +
  geom_point() +
  geom_line() +
  xlab(expression(italic(n))) +
  ylab("") +
  scale_color_manual(name = "Asymptotic standard deviation", values = c(21, 22),
                     breaks = c("sigma", "nu"),
                     labels = c(expression(sigma), expression(nu))) +
  theme_hc() +
  theme(text = element_text(family = "serif", size = 25))

