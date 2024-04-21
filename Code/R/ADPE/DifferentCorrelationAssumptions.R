library(ggplot2)
library(ggthemes)
library(ggbreak)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("InformationTheory.R")
source("Bandt-Pompe.R")
source("AsymptoticEntropyVariancesComparison.R")
# require(reshape2)

### FUNCTIONS

# Function to find the normalized Shannon entropy

HShannon <- function(p){
  
  prob <- p[p > 0]
  N <- length(p)
  H <- -sum(prob * log(prob)) / log(N)
  
  return(H)
}


# Function to find the p-value

pval_ind <- function(p1, p2, n1, n2){
  
  H1 <- HShannon(p1)
  H2 <- HShannon(p2)
  
  V1 <- v.n.q(n1,3,p1)
  V2 <- v.n.q(n2,3,p2)
  
  epsilon <- abs(H1 - H2)
  sigma <- sqrt(V1 + V2)
  
  pv <- 2 - 2 * pnorm(epsilon / sigma)
  
  return(pv)
}

pval_corr <- function(p1, p2, n1, n2){
  
  H1 <- HShannon(p1)
  H2 <- HShannon(p2)
  
  V1 <- sigma.n.q(n1,3,p1)
  V2 <- sigma.n.q(n2,3,p2)
  
  epsilon <- abs(H1 - H2)
  sigma <- sqrt(V1 + V2)
  
  pv <- 2 - 2 * pnorm(epsilon / sigma)
  
  return(pv)
}


# Load Data --------------------------------------------------------------------
meteo = read.csv("../../Data/meteo/3010488.csv")

meteo_dublin = meteo[meteo$NAME == "DUBLIN PHOENIX PARK, EI",][16298:26297, "TMAX"]
meteo_miami = meteo[meteo$NAME == "MIAMI INTERNATIONAL AIRPORT, FL US",][16298:26297, "TMAX"]
meteo_edinburgh = meteo[meteo$NAME == "EDINBURGH ROYAL BOTANIC GARDE, UK",][11916:21915, "TMAX"]

# D = 3

df_data = data.frame(
  "Dublin" = meteo_dublin,
  "Edinburgh" = meteo_edinburgh,
  "Miami" = meteo_miami
  #date
)

    
probs_dublin <- bandt.pompe(unlist(df_data["Dublin"]), dimension=3, delay=1)
probs_edinburgh <- bandt.pompe(unlist(df_data["Edinburgh"]), dimension=3, delay=1)
probs_miami <- bandt.pompe(unlist(df_data["Miami"]), dimension=3, delay=1)

n_d<- length(unlist(df_data["Dublin"])) -2
n_e<- length(unlist(df_data["Edinburgh"])) -2
n_m<- length(unlist(df_data["Miami"])) -2

## Hypothesis Testing

# Dublin vs Edinburgh
p_DE_i<-pval_ind(probs_dublin,probs_edinburgh,n_d,n_e)
p_DE_c<-pval_corr(probs_dublin,probs_edinburgh,n_d,n_e)

# Dublin vs Miami
p_DM_i<-pval_ind(probs_dublin,probs_miami,n_d,n_m)
p_DM_c<-pval_corr(probs_dublin,probs_miami,n_d,n_m)

# Edinburgh vs Miami
p_EM_i<-pval_ind(probs_edinburgh,probs_miami,n_e,n_m)
p_EM_c<-pval_corr(probs_edinburgh,probs_miami,n_e,n_m)

#####################################################################

# Minumun temperature

meteo_dublin_min = meteo[meteo$NAME == "DUBLIN PHOENIX PARK, EI",][16298:26297, "TMIN"]
meteo_miami_min = meteo[meteo$NAME == "MIAMI INTERNATIONAL AIRPORT, FL US",][16298:26297, "TMIN"]
meteo_edinburgh_min = meteo[meteo$NAME == "EDINBURGH ROYAL BOTANIC GARDE, UK",][11916:21915, "TMIN"]

# D = 3

df_data_min = data.frame(
  "Dublin" = meteo_dublin_min,
  "Edinburgh" = meteo_edinburgh_min,
  "Miami" = meteo_miami_min
  #date
)


probs_dublin_min <- bandt.pompe(unlist(df_data_min["Dublin"]), dimension = 3, delay = 1)
probs_edinburgh_min <- bandt.pompe(unlist(df_data_min["Edinburgh"]), dimension = 3, delay = 1)
probs_miami_min <- bandt.pompe(unlist(df_data_min["Miami"]), dimension = 3, delay = 1)

n_d_min <- length(unlist(df_data_min["Dublin"])) - 2
n_e_min <- length(unlist(df_data_min["Edinburgh"])) - 2
n_m_min <- length(unlist(df_data_min["Miami"])) - 2

## Hypothesis Testing

# Dublin vs Edinburgh
p_DE_i_min <- pval_ind(probs_dublin_min, probs_edinburgh_min, n_d_min, n_e_min)
p_DE_c_min <- pval_corr(probs_dublin_min, probs_edinburgh_min, n_d_min, n_e_min)

# Dublin vs Miami
p_DM_i_min <- pval_ind(probs_dublin_min, probs_miami_min, n_d_min, n_m_min)
p_DM_c_min <- pval_corr(probs_dublin_min, probs_miami_min, n_d_min, n_m_min)

# Edinburgh vs Miami
p_EM_i_min <- pval_ind(probs_edinburgh_min, probs_miami_min, n_e_min, n_m_min)
p_EM_c_min <- pval_corr(probs_edinburgh_min, probs_miami_min, n_e_min, n_m_min)


#####################################################################

# Precipitations

meteo_dublin_pre = meteo[meteo$NAME == "DUBLIN PHOENIX PARK, EI",][16298:26297, "PRCP"]
meteo_miami_pre = meteo[meteo$NAME == "MIAMI INTERNATIONAL AIRPORT, FL US",][16298:26297, "PRCP"]
meteo_edinburgh_pre = meteo[meteo$NAME == "EDINBURGH ROYAL BOTANIC GARDE, UK",][11916:21915, "PRCP"]

# D = 3

df_data_pre = data.frame(
  "Dublin" = meteo_dublin_pre,
  "Edinburgh" = meteo_edinburgh_pre,
  "Miami" = meteo_miami_pre
  #date
)


probs_dublin_pre <- bandt.pompe(unlist(df_data_pre["Dublin"]), dimension = 3, delay = 1)
probs_edinburgh_pre <- bandt.pompe(unlist(df_data_pre["Edinburgh"]), dimension = 3, delay = 1)
probs_miami_pre <- bandt.pompe(unlist(df_data_pre["Miami"]), dimension = 3, delay = 1)

n_d_pre <- length(unlist(df_data_pre["Dublin"])) - 2
n_e_pre <- length(unlist(df_data_pre["Edinburgh"])) - 2
n_m_pre <- length(unlist(df_data_pre["Miami"])) - 2

## Hypothesis Testing

# Dublin vs Edinburgh
p_DE_i_pre <- pval_ind(probs_dublin_pre, probs_edinburgh_pre, n_d_pre, n_e_pre)
p_DE_c_pre <- pval_corr(probs_dublin_pre, probs_edinburgh_pre, n_d_pre, n_e_pre)

# Dublin vs Miami
p_DM_i_pre <- pval_ind(probs_dublin_pre, probs_miami_pre, n_d_pre, n_m_pre)
p_DM_c_pre <- pval_corr(probs_dublin_pre, probs_miami_pre, n_d_pre, n_m_pre)

# Edinburgh vs Miami
p_EM_i_pre <- pval_ind(probs_edinburgh_pre, probs_miami_pre, n_e_pre, n_m_pre)
p_EM_c_pre <- pval_corr(probs_edinburgh_pre, probs_miami_pre, n_e_pre, n_m_pre)


# df <- data.frame(Day = c(1:(n_d_pre+2), 1:(n_e_pre+2), 1:(n_m_pre+2)),
#                  PRCP = c(unlist(df_data_pre["Dublin"]), 
#                           unlist(df_data_pre["Edinburgh"]),
#                           unlist(df_data_pre["Miami"])),
#                  City = c(rep("Dublin", n_d_pre+2), 
#                           rep("Edinburgh", n_e_pre+2),
#                           rep("Miami", n_m_pre+2)))
# 
# ggplot(df, aes(x = Day, y = PRCP, col = City)) +
#   geom_line()
# 
# df_red <- data.frame(Day = rep(1:500, 3),
#                  PRCP = c(unlist(df_data_pre["Dublin"])[1:500], 
#                           unlist(df_data_pre["Edinburgh"])[1:500],
#                           unlist(df_data_pre["Miami"])[1:500]),
#                  City = c(rep("Dublin", 500), 
#                           rep("Edinburgh", 500),
#                           rep("Miami", 500)))
# 
# ggplot(df_red, aes(x = Day, y = PRCP, col = City)) +
#   geom_line()
# 
# df_DE <- data.frame(Day = rep(1:500, 2),
#                      PRCP = c(unlist(df_data_pre["Dublin"])[1:500], 
#                               unlist(df_data_pre["Edinburgh"])[1:500]),
#                      City = c(rep("Dublin", 500), 
#                               rep("Edinburgh", 500)))
# 
# ggplot(df_DE, aes(x = Day, y = PRCP, col = City)) +
#   geom_line()

#######################################################################

# Table

pvals <- rbind(c(p_DE_c_min, p_DM_c_min, p_EM_c_min),
               c(p_DE_i_min, p_DM_i_min, p_EM_i_min),
               c(p_DE_c, p_DM_c, p_EM_c),
               c(p_DE_i, p_DM_i, p_EM_i),
               c(p_DE_c_pre, p_DM_c_pre, p_EM_c_pre),
               c(p_DE_i_pre, p_DM_i_pre, p_EM_i_pre)
               )

#######################################################################

### Confidence intervals ###

n2 <- 10000

# Means computation for the whole series

mean_dublin_min <- HShannon(probs_dublin_min)
mean_edinburgh_min <- HShannon(probs_edinburgh_min)
mean_miami_min <- HShannon(probs_miami_min)

mean_dublin_max <- HShannon(probs_dublin)
mean_edinburgh_max <- HShannon(probs_edinburgh)
mean_miami_max <- HShannon(probs_miami)

mean_dublin_prcp <- HShannon(probs_dublin_pre)
mean_edinburgh_prcp <- HShannon(probs_edinburgh_pre)
mean_miami_prcp <- HShannon(probs_miami_pre)


# Standard deviation computation for the whole series with the true model

sd_dublin_min_c <- sqrt(sigma.n.q(n2, 3, probs_dublin_min))
sd_edinburgh_min_c <- sqrt(sigma.n.q(n2, 3, probs_edinburgh_min))
sd_miami_min_c <- sqrt(sigma.n.q(n2, 3, probs_miami_min))

sd_dublin_max_c <- sqrt(sigma.n.q(n2, 3, probs_dublin))
sd_edinburgh_max_c <- sqrt(sigma.n.q(n2, 3, probs_edinburgh))
sd_miami_max_c <- sqrt(sigma.n.q(n2, 3, probs_miami))

sd_dublin_prcp_c <- sqrt(sigma.n.q(n2, 3, probs_dublin_pre))
sd_edinburgh_prcp_c <- sqrt(sigma.n.q(n2, 3, probs_edinburgh_pre))
sd_miami_prcp_c <- sqrt(sigma.n.q(n2, 3, probs_miami_pre))

# Standard deviation computation for the whole series with the multinomial model

sd_dublin_min_i <- sqrt(v.n.q(n2, 3, probs_dublin_min))
sd_edinburgh_min_i <- sqrt(v.n.q(n2, 3, probs_edinburgh_min))
sd_miami_min_i <- sqrt(v.n.q(n2, 3, probs_miami_min))

sd_dublin_max_i <- sqrt(v.n.q(n2, 3, probs_dublin))
sd_edinburgh_max_i <- sqrt(v.n.q(n2, 3, probs_edinburgh))
sd_miami_max_i <- sqrt(v.n.q(n2, 3, probs_miami))

sd_dublin_prcp_i <- sqrt(v.n.q(n2, 3, probs_dublin_pre))
sd_edinburgh_prcp_i <- sqrt(v.n.q(n2, 3, probs_edinburgh_pre))
sd_miami_prcp_i <- sqrt(v.n.q(n2, 3, probs_miami_pre))

# CI as data frames

df_min <- data.frame(
  Limits = c(mean_dublin_min - 1.96 * sd_dublin_min_c / sqrt(n2),
             mean_dublin_min + 1.96 * sd_dublin_min_c / sqrt(n2),
             mean_dublin_min - 1.96 * sd_dublin_min_i / sqrt(n2),
             mean_dublin_min + 1.96 * sd_dublin_min_i / sqrt(n2),
             mean_edinburgh_min - 1.96 * sd_edinburgh_min_c / sqrt(n2),
             mean_edinburgh_min + 1.96 * sd_edinburgh_min_c / sqrt(n2),
             mean_edinburgh_min - 1.96 * sd_edinburgh_min_i / sqrt(n2),
             mean_edinburgh_min + 1.96 * sd_edinburgh_min_i / sqrt(n2),
             mean_miami_min - 1.96 * sd_miami_min_c / sqrt(n2),
             mean_miami_min + 1.96 * sd_miami_min_c / sqrt(n2),
             mean_miami_min - 1.96 * sd_miami_min_i / sqrt(n2),
             mean_miami_min + 1.96 * sd_miami_min_i / sqrt(n2)),
  Level = c(rep(0,2), rep(1,2), rep(0,2), rep(1,2), rep(0,2), rep(1,2)), 
  Case = c(rep("Dublin (true)",2), rep("Dublin (mult.)",2),
           rep("Edinburgh (true)",2), rep("Edinburgh (mult.)",2),
           rep("Miami (true)",2), rep("Miami (mult.)",2))) 

df_max <- data.frame(
  Limits = c(mean_dublin_max - 1.96 * sd_dublin_max_c / sqrt(n2),
             mean_dublin_max + 1.96 * sd_dublin_max_c / sqrt(n2),
             mean_dublin_max - 1.96 * sd_dublin_max_i / sqrt(n2),
             mean_dublin_max + 1.96 * sd_dublin_max_i / sqrt(n2),
             mean_edinburgh_max - 1.96 * sd_edinburgh_max_c / sqrt(n2),
             mean_edinburgh_max + 1.96 * sd_edinburgh_max_c / sqrt(n2),
             mean_edinburgh_max - 1.96 * sd_edinburgh_max_i / sqrt(n2),
             mean_edinburgh_max + 1.96 * sd_edinburgh_max_i / sqrt(n2),
             mean_miami_max - 1.96 * sd_miami_max_c / sqrt(n2),
             mean_miami_max + 1.96 * sd_miami_max_c / sqrt(n2),
             mean_miami_max - 1.96 * sd_miami_max_i / sqrt(n2),
             mean_miami_max + 1.96 * sd_miami_max_i / sqrt(n2)),
  Level = c(rep(0,2), rep(1,2), rep(0,2), rep(1,2), rep(0,2), rep(1,2)), 
  Case = c(rep("Dublin (true)",2), rep("Dublin (mult.)",2),
           rep("Edinburgh (true)",2), rep("Edinburgh (mult.)",2),
           rep("Miami (true)",2), rep("Miami (mult.)",2))) 

df_prcp <- data.frame(
  Limits = c(mean_dublin_prcp - 1.96 * sd_dublin_prcp_c / sqrt(n2),
             mean_dublin_prcp + 1.96 * sd_dublin_prcp_c / sqrt(n2),
             mean_dublin_prcp - 1.96 * sd_dublin_prcp_i / sqrt(n2),
             mean_dublin_prcp + 1.96 * sd_dublin_prcp_i / sqrt(n2),
             mean_edinburgh_prcp - 1.96 * sd_edinburgh_prcp_c / sqrt(n2),
             mean_edinburgh_prcp + 1.96 * sd_edinburgh_prcp_c / sqrt(n2),
             mean_edinburgh_prcp - 1.96 * sd_edinburgh_prcp_i / sqrt(n2),
             mean_edinburgh_prcp + 1.96 * sd_edinburgh_prcp_i / sqrt(n2),
             mean_miami_prcp - 1.96 * sd_miami_prcp_c / sqrt(n2),
             mean_miami_prcp + 1.96 * sd_miami_prcp_c / sqrt(n2),
             mean_miami_prcp - 1.96 * sd_miami_prcp_i / sqrt(n2),
             mean_miami_prcp + 1.96 * sd_miami_prcp_i / sqrt(n2)),
  Level = c(rep(0,2), rep(1,2), rep(0,2), rep(1,2), rep(0,2), rep(1,2)), 
  Case = c(rep("Dublin (true)",2), rep("Dublin (mult.)",2),
           rep("Edinburgh (true)",2), rep("Edinburgh (mult.)",2),
           rep("Miami (true)",2), rep("Miami (mult.)",2))) 

# CI plots

plot_min <- ggplot(df_min, aes(x = Limits, y = Level, col = Case)) +
  # geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  # geom_hline(yintercept=1, linetype="dashed", color = "gray") +
  geom_segment(aes(x = 0.9573, y = 0, xend = 0.9573, yend = 1), linetype="dashed", color = "white") +
  geom_segment(aes(x = mean_dublin_min, y = 0, xend = mean_dublin_min, yend = 1), linetype="dashed", color = "gray") +
  geom_segment(aes(x = mean_edinburgh_min, y = 0, xend = mean_edinburgh_min, yend = 1), linetype="dashed", color = "gray") +
  geom_segment(aes(x = mean_miami_min, y = 0, xend = mean_miami_min, yend = 1), linetype="dashed", color = "gray") +
  geom_line(linewidth=1) +
  geom_point(shape = 124, size = 2.5) +
  scale_color_manual(name = "Model",
                     breaks = c("Dublin (mult.)", "Dublin (true)",
                                "Edinburgh (mult.)", "Edinburgh (true)",
                                "Miami (mult.)", "Miami (true)"),
                     values = c("Dublin (mult.)" = 22, "Dublin (true)" = 21,
                                "Edinburgh (mult.)" = 22, "Edinburgh (true)" = 21,
                                "Miami (mult.)" = 22, "Miami (true)" = 21),
                     labels = rep(c("Multinomial", "True"), 3)) +
  geom_point(aes(x = mean_dublin_min, y = 1), color = 22, size = 1) +
  geom_point(aes(x = mean_dublin_min, y = 0), color = 21, size = 1) +
  geom_point(aes(x = mean_edinburgh_min, y = 1), color = 22, size = 1) +
  geom_point(aes(x = mean_edinburgh_min, y = 0), color = 21, size = 1) +
  geom_point(aes(x = mean_miami_min, y = 1), color = 22, size = 1) +
  geom_point(aes(x = mean_miami_min, y = 0), color = 21, size = 1) +
  xlab("") +
  ylab("") +
  ylab("Minimum temperature") +
  scale_x_cut(c(0.9577, 0.97665, 0.977, 0.9793), which = c(2,4), scales = c(0,0)) +
  xlim(c(0.9572, 0.9797)) +
  scale_x_continuous(breaks = c(0.9574, 0.9575, 0.9576, 0.9767, 
                                0.9768, 0.9769, 0.9794, 0.9795)) +
  theme_pander() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(family = "serif"),
        legend.title = element_blank(),
        legend.position = "none") +
  annotate("text", x = 0.95747, y = 1.1, label = "Dublin", size = 6, family = "serif") +
  annotate("text", x = 0.976797, y = 1.1, label = "Edinburgh", size = 6, family = "serif") +
  annotate("text", x = 0.97946, y = 1.1, label = "Miami", size = 6, family = "serif")

plot_max <-
  ggplot(df_max, aes(x = Limits, y = Level, col = Case)) +
  # geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  # geom_hline(yintercept=1, linetype="dashed", color = "gray") +
  geom_segment(aes(x = 0.91785, y = 0, xend = 0.91785, yend = 1), linetype="dashed", color = "white") +
  geom_segment(aes(x = mean_dublin_max, y = 0, xend = mean_dublin_max, yend = 1), linetype="dashed", color = "gray") +
  geom_segment(aes(x = mean_edinburgh_max, y = 0, xend = mean_edinburgh_max, yend = 1), linetype="dashed", color = "gray") +
  geom_segment(aes(x = mean_miami_max, y = 0, xend = mean_miami_max, yend = 1), linetype="dashed", color = "gray") +
  geom_line(linewidth=1) +
  geom_point(shape = 124, size = 2.5) +
  scale_color_manual(name = "Model",
                     breaks = c("Dublin (mult.)", "Dublin (true)",
                                "Edinburgh (mult.)", "Edinburgh (true)",
                                "Miami (mult.)", "Miami (true)"),
                     values = c("Dublin (mult.)" = 22, "Dublin (true)" = 21,
                                "Edinburgh (mult.)" = 22, "Edinburgh (true)" = 21,
                                "Miami (mult.)" = 22, "Miami (true)" = 21),
                     labels = rep(c("Multinomial", "True"), 3)) +
  geom_point(aes(x = mean_dublin_max, y = 1), color = 22, size = 1) +
  geom_point(aes(x = mean_dublin_max, y = 0), color = 21, size = 1) +
  geom_point(aes(x = mean_edinburgh_max, y = 1), color = 22, size = 1) +
  geom_point(aes(x = mean_edinburgh_max, y = 0), color = 21, size = 1) +
  geom_point(aes(x = mean_miami_max, y = 1), color = 22, size = 1) +
  geom_point(aes(x = mean_miami_max, y = 0), color = 21, size = 1) +
  xlab("") +
  ylab("") +
  ylab("Maximum temperature") +
  scale_x_cut(c(0.9183, 0.96735, 0.9677, 0.9714), which = c(2,4), scales = c(0,0)) +
  xlim(c(0.9178, 0.97175)) +
  scale_x_continuous(breaks = c(0.9179, 0.9180, 0.9181, 0.9182, 0.9674,
                                0.9675, 0.9676, 0.9715, 0.9716, 0.9717)) +
  theme_pander() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(family = "serif"),
        legend.title = element_blank(),
        legend.position = "none") +
  annotate("text", x = 0.918065, y = 1.1, label = "Dublin", size = 6, family = "serif") +
  annotate("text", x = 0.96755, y = 1.1, label = "Edinburgh", size = 6, family = "serif") +
  annotate("text", x = 0.97159, y = 1.1, label = "Miami", size = 6, family = "serif")


plot_prcp <-
  ggplot(df_prcp, aes(x = Limits, y = Level, col = Case)) +
  # geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  # geom_hline(yintercept=1, linetype="dashed", color = "gray") +
  geom_segment(aes(x = 0.7836, y = 0, xend = 0.7836, yend = 1), linetype="dashed", color = "white") +
  geom_segment(aes(x = mean_dublin_prcp, y = 0, xend = mean_dublin_prcp, yend = 1), linetype="dashed", color = "gray") +
  geom_segment(aes(x = mean_edinburgh_prcp, y = 0, xend = mean_edinburgh_prcp, yend = 1), linetype="dashed", color = "gray") +
  geom_segment(aes(x = mean_miami_prcp, y = 0, xend = mean_miami_prcp, yend = 1), linetype="dashed", color = "gray") +
  geom_line(linewidth=1) +
  geom_point(shape = 124, size = 2.5) +
  scale_color_manual(name = "Model",
                     breaks = c("Dublin (mult.)", "Dublin (true)",
                                "Edinburgh (mult.)", "Edinburgh (true)",
                                "Miami (mult.)", "Miami (true)"),
                     values = c("Dublin (mult.)" = 22, "Dublin (true)" = 21,
                                "Edinburgh (mult.)" = 22, "Edinburgh (true)" = 21,
                                "Miami (mult.)" = 22, "Miami (true)" = 21),
                     labels = rep(c("Multinomial", "True"), 3)) +
  geom_point(aes(x = mean_dublin_prcp, y = 1), color = 22, size = 1) +
  geom_point(aes(x = mean_dublin_prcp, y = 0), color = 21, size = 1) +
  geom_point(aes(x = mean_edinburgh_prcp, y = 1), color = 22, size = 1) +
  geom_point(aes(x = mean_edinburgh_prcp, y = 0), color = 21, size = 1) +
  geom_point(aes(x = mean_miami_prcp, y = 1), color = 22, size = 1) +
  geom_point(aes(x = mean_miami_prcp, y = 0), color = 21, size = 1) +
  xlab("") +
  ylab("") +
  ylab("Daily precipitation") +
  scale_x_cut(c(0.7843, 0.9028, 0.9034, 0.9237), which = c(2,4), scales = c(0,0)) +
  xlim(c(0.78355, 0.92425)) +
  scale_x_continuous(breaks = c(0.7838, 0.7840, 0.7842, 0.9029,
                                0.9031, 0.9033, 0.9238, 0.9240)) +
  theme_pander() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 15),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    text = element_text(family = "serif"),
    legend.title = element_blank(),
    legend.position = "none"
    ) +
  annotate("text", x = 0.78393, y = 1.1, label = "Dublin", size = 6, family = "serif") +
  annotate("text", x = 0.90311, y = 1.1, label = "Edinburgh", size = 6, family = "serif") +
  annotate("text", x = 0.92396, y = 1.1, label = "Miami", size = 6, family = "serif")


plot_min
plot_max
plot_prcp

########################################################################

## Data for table

feature <- c(rep("Minimum temperature", 3),
             rep("Maximum temperature", 3),
             rep("Daily precipitation", 3))

city <- rep(c("Dublin", "Edinburgh", "Miami"), 3)

means <- c(mean_dublin_min, mean_edinburgh_min, mean_miami_min,
           mean_dublin_max, mean_edinburgh_max, mean_miami_max,
           mean_dublin_prcp, mean_edinburgh_prcp, mean_miami_prcp)

ltrue <- c(2 * 1.96 * sd_dublin_min_c / sqrt(n2),
                2 * 1.96 * sd_edinburgh_min_c / sqrt(n2),
                2 * 1.96 * sd_miami_min_c / sqrt(n2),
                2 * 1.96 * sd_dublin_max_c / sqrt(n2),
                2 * 1.96 * sd_edinburgh_max_c / sqrt(n2),
                2 * 1.96 * sd_miami_max_c / sqrt(n2),
                2 * 1.96 * sd_dublin_prcp_c / sqrt(n2),
                2 * 1.96 * sd_edinburgh_prcp_c / sqrt(n2),
                2 * 1.96 * sd_miami_prcp_c / sqrt(n2))

lmult <- c(2 * 1.96 * sd_dublin_min_i / sqrt(n2),
                2 * 1.96 * sd_edinburgh_min_i / sqrt(n2),
                2 * 1.96 * sd_miami_min_i / sqrt(n2),
                2 * 1.96 * sd_dublin_max_i / sqrt(n2),
                2 * 1.96 * sd_edinburgh_max_i / sqrt(n2),
                2 * 1.96 * sd_miami_max_i / sqrt(n2),
                2 * 1.96 * sd_dublin_prcp_i / sqrt(n2),
                2 * 1.96 * sd_edinburgh_prcp_i / sqrt(n2),
                2 * 1.96 * sd_miami_prcp_i / sqrt(n2))



info <- data.frame(Feature = feature,
                   City = city,
                   Mean = means,
                   LengthTrue = ltrue,
                   LengthMult = lmult,
                   Difference = ltrue - lmult)

