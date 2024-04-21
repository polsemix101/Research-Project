library(statcomp)
library(ggplot2)
library(ggthemes)
library(scales)
library(reshape2)
library(rstudioapi)
library(cowplot)

setwd(dirname(getActiveDocumentContext()$path))

source("Bandt-Pompe.R")

#####################################################################################

## Read, organize, and select meteo data

meteo = read.csv("../../Data/meteo/3010488.csv")

# Data date

date <- as.Date(meteo[meteo$NAME == "DUBLIN PHOENIX PARK, EI",][16298:26297, "DATE"],
                "%Y-%m-%d")
range(date)


# Data min

meteo_dublin_min = meteo[meteo$NAME == "DUBLIN PHOENIX PARK, EI",][16298:26297, "TMIN"]
meteo_miami_min = meteo[meteo$NAME == "MIAMI INTERNATIONAL AIRPORT, FL US",][16298:26297, "TMIN"]
meteo_edinburgh_min = meteo[meteo$NAME == "EDINBURGH ROYAL BOTANIC GARDE, UK",][11916:21915, "TMIN"]

df_data_min = data.frame(
  "Dublin" = meteo_dublin_min,
  "Edinburgh" = meteo_edinburgh_min,
  "Miami" = meteo_miami_min,
  date
)

meteo.df.melt.min <- melt(df_data_min, 
                      measure.vars = 1:3, 
                      variable.name = "City", 
                      value.name = "T")

# Data max

meteo_dublin_max = meteo[meteo$NAME == "DUBLIN PHOENIX PARK, EI",][16298:26297, "TMAX"]
meteo_miami_max = meteo[meteo$NAME == "MIAMI INTERNATIONAL AIRPORT, FL US",][16298:26297, "TMAX"]
meteo_edinburgh_max = meteo[meteo$NAME == "EDINBURGH ROYAL BOTANIC GARDE, UK",][11916:21915, "TMAX"]

df_data_max = data.frame(
  "Dublin" = meteo_dublin_max,
  "Edinburgh" = meteo_edinburgh_max,
  "Miami" = meteo_miami_max,
  date
)

meteo.df.melt.max <- melt(df_data_max, 
                      measure.vars = 1:3, 
                      variable.name = "City", 
                      value.name = "T")

# Data precipitation

meteo_dublin_prcp = meteo[meteo$NAME == "DUBLIN PHOENIX PARK, EI",][16298:26297, "PRCP"]
meteo_miami_prcp = meteo[meteo$NAME == "MIAMI INTERNATIONAL AIRPORT, FL US",][16298:26297, "PRCP"]
meteo_edinburgh_prcp = meteo[meteo$NAME == "EDINBURGH ROYAL BOTANIC GARDE, UK",][11916:21915, "PRCP"]

df_data_prcp = data.frame(
  "Dublin" = meteo_dublin_prcp,
  "Edinburgh" = meteo_edinburgh_prcp,
  "Miami" = meteo_miami_prcp,
  date
)


meteo.df.melt.prcp <- melt(df_data_prcp, 
                      measure.vars = 1:3, 
                      variable.name = "City", 
                      value.name = "T")

# Plots

ggplot(meteo.df.melt.min, aes(x=date, y=T, col=City)) + 
  geom_line(linewidth=.2) +
  xlab("Date") +
  ylab("Minimum Temperature [ºF]") +
  facet_grid(City ~.) +
  theme_pander() + 
  theme(text=element_text(size=8,
                          family="serif"),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title = element_blank()
  )  +
  scale_color_manual(values = c("Dublin" = "#1b9e77", "Edinburgh" = "#7570b3", "Miami" = "#d95f02")) +
  scale_x_date(date_breaks="1 year", date_labels="%Y")

ggsave(file="../../Figures/PDF/MinTemperatureTimeSeries.pdf",
       width=30, height=9, units="cm")

ggplot(meteo.df.melt.max, aes(x=date, y=T, col=City)) + 
  geom_line(linewidth=.2) +
  xlab("Date") +
  ylab("Maximum Temperature [ºF]") +
  facet_grid(City ~.) +
  theme_pander() + 
  theme(text=element_text(size=8,
                          family="serif"),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title = element_blank()
  )  +
  scale_color_manual(values = c("Dublin" = "#1b9e77", "Edinburgh" = "#7570b3", "Miami" = "#d95f02")) +
  scale_x_date(date_breaks="1 year", date_labels="%Y")

ggsave(file="../../Figures/PDF/MaxTemperatureTimeSeries.pdf",
       width=30, height=9, units="cm")

ggplot(meteo.df.melt.prcp, aes(x=date, y=T, col=City)) + 
  geom_line(linewidth=.2) +
  xlab("Date") +
  ylab("Precipitation") +
  facet_grid(City ~., scales="free") +
  theme_pander() + 
  theme(text=element_text(size=8,
                          family="serif"),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title = element_blank()
  )  +
  scale_color_manual(values = c("Dublin" = "#1b9e77", "Edinburgh" = "#7570b3", "Miami" = "#d95f02")) +
  scale_x_date(date_breaks="1 year", date_labels="%Y")
  

ggsave(file="../../Figures/PDF/PrecipitationTimeSeries.pdf",
       width=30, height=9, units="cm")

#############################################################################

## Histograms of ordinal patterns

histD3 <- function(series){
  
  p.patterns = formationPattern(series, D = 3, tau = 1, 0)
  n.symbols = dim(p.patterns)[1]
  symbol = matrix(c(0,1,2,
                    0,2,1,
                    1,0,2,
                    1,2,0,
                    2,0,1,
                    2,1,0), ncol = 3, byrow = TRUE)
  index.rep = array(0, n.symbols)
  
  for(i in 1:n.symbols){
    for(j in 1:6){
      if(all(p.patterns[i,] == symbol[j, ])){
        index.rep[i]=j
        break
      }
    }
  }
  
  index.rep = index.rep[1:n.symbols]
  index.rep = data.frame(i = index.rep)
  return(index.rep)
}


# Min

hist_dublin_min = histD3(meteo_dublin_min)
hist_edinburgh_min = histD3(meteo_edinburgh_min)
hist_miami_min = histD3(meteo_miami_min)

seqs.symbols.min <- data.frame(
  Dublin = unlist(hist_dublin_min),
  Edinburgh = unlist(hist_edinburgh_min),
  Miami = unlist(hist_miami_min)
)

meltMin <- melt(seqs.symbols.min, variable.name = "City")

PropMin <- ggplot(meltMin, aes(x=value, col=City, fill=City)) +
  geom_histogram(aes(y=..density..), alpha=.9,
                 bins=6, binwidth=.95, 
                 position=position_dodge(.85)) +
  theme_pander() + 
  theme(text=element_text(size=24,
                          family="serif"),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title = element_blank()
  ) +
  xlab("Pattern") + ylab("Proportion") +
  scale_color_manual(values = c("Dublin" = "#1b9e77", "Edinburgh" = "#7570b3", "Miami" = "#d95f02")) +
  scale_fill_manual(values = c("Dublin" = "#1b9e77", "Edinburgh" = "#7570b3", "Miami" = "#d95f02")) +
  scale_x_continuous(breaks=1:6,
                     labels=c(expression(pi[1]),
                              expression(pi[2]),
                              expression(pi[3]),
                              expression(pi[4]),
                              expression(pi[5]),
                              expression(pi[6]))) +
  ylim(0, .57)


ggsave(file="../../Figures/PDF/histMin.pdf", width=9, height=5, units="cm")

# Max

hist_dublin_max = histD3(meteo_dublin_max)
hist_edinburgh_max = histD3(meteo_edinburgh_max)
hist_miami_max = histD3(meteo_miami_max)

seqs.symbols.max <- data.frame(
  Dublin = unlist(hist_dublin_max),
  Edinburgh = unlist(hist_edinburgh_max),
  Miami = unlist(hist_miami_max)
)

meltMax <- melt(seqs.symbols.max, variable.name = "City")

PropMax <- ggplot(meltMax, aes(x=value, col=City, fill=City)) +
  geom_histogram(aes(y=..density..), alpha=.9,
                 bins=6, binwidth=.95, 
                 position=position_dodge(.85)) +
  theme_pander() + 
  theme(text=element_text(size=24,
                          family="serif"),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title = element_blank()
  ) +
  xlab("Pattern") + ylab("Proportion") +
  scale_color_manual(values = c("Dublin" = "#1b9e77", "Edinburgh" = "#7570b3", "Miami" = "#d95f02")) +
  scale_fill_manual(values = c("Dublin" = "#1b9e77", "Edinburgh" = "#7570b3", "Miami" = "#d95f02")) +
  scale_x_continuous(breaks=1:6,
                     labels=c(expression(pi[1]),
                              expression(pi[2]),
                              expression(pi[3]),
                              expression(pi[4]),
                              expression(pi[5]),
                              expression(pi[6]))) +
  ylim(0, .57)


ggsave(file="../../Figures/PDF/histMax.pdf", width=9, height=5, units="cm")

# Prcp

hist_dublin_prcp = histD3(meteo_dublin_prcp)
hist_edinburgh_prcp = histD3(meteo_edinburgh_prcp)
hist_miami_prcp = histD3(meteo_miami_prcp)

seqs.symbols.prcp <- data.frame(
  Dublin = unlist(hist_dublin_prcp),
  Edinburgh = unlist(hist_edinburgh_prcp),
  Miami = unlist(hist_miami_prcp)
)

meltPrcp <- melt(seqs.symbols.prcp, variable.name = "City")

PropPrecipitation <- ggplot(meltPrcp, aes(x=value, col=City, fill=City)) +
  geom_histogram(aes(y=..density..), alpha=.9,
                 bins=6, binwidth=.95, 
                 position=position_dodge(.85)) +
  theme_pander() + 
  theme(text=element_text(size=24,
                          family="serif"),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title = element_blank()
  ) +
  xlab("Pattern") + ylab("Proportion") +
  scale_color_manual(values = c("Dublin" = "#1b9e77", "Edinburgh" = "#7570b3", "Miami" = "#d95f02")) +
  scale_fill_manual(values = c("Dublin" = "#1b9e77", "Edinburgh" = "#7570b3", "Miami" = "#d95f02")) +
  scale_x_continuous(breaks=1:6,
                     labels=c(expression(pi[1]),
                              expression(pi[2]),
                              expression(pi[3]),
                              expression(pi[4]),
                              expression(pi[5]),
                              expression(pi[6]))) +
  ylim(0, .57)

ggsave(file="../../Figures/PDF/histPrcp.pdf", width=9, height=5, units="cm")



