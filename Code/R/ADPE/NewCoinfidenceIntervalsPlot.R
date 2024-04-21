# Correct (correlated) model
AsymptoticQuantities <- data.frame(
  means=c(mean_dublin_min, mean_edinburgh_min, mean_miami_min, 
          mean_dublin_max, mean_edinburgh_max, mean_miami_max, 
          mean_dublin_prcp, mean_edinburgh_prcp,mean_miami_prcp
          ),
  standard.deviations=c(
    sd_dublin_min_c, sd_edinburgh_min_c, sd_miami_min_c,
    sd_dublin_max_c, sd_edinburgh_max_c, sd_miami_max_c,
    sd_dublin_prcp_c, sd_edinburgh_prcp_c, sd_miami_prcp_c)
)

AsymptoticQuantities$City <- rep(c("Dublin", "Edinburgh", "Miami"), 3)
AsymptoticQuantities$Measurement <- rep(c("Minimum", "Maximum", "Precipitation"), each=3)

# Multinomial model

ApproximateQuantities <- data.frame(
  means=c(mean_dublin_min, mean_edinburgh_min, mean_miami_min, 
          mean_dublin_max, mean_edinburgh_max, mean_miami_max, 
          mean_dublin_prcp, mean_edinburgh_prcp,mean_miami_prcp
  ),
  standard.deviations=c(sd_dublin_min_i, sd_edinburgh_min_i, sd_miami_min_i,
                        sd_dublin_max_i, sd_edinburgh_max_i, sd_miami_max_i,
                        sd_dublin_prcp_i, sd_edinburgh_prcp_i, sd_miami_prcp_i)
  )
ApproximateQuantities$City <- rep(c("Dublin", "Edinburgh", "Miami"), 3)
ApproximateQuantities$Measurement <- rep(c("Minimum", "Maximum", "Precipitation"), each=3)

# Plot
qnormal <- qnorm((1-.01/2))
ggplot(subset(AsymptoticQuantities, City!="Miami"), aes(x=City, y=means)) +
  geom_point(shape=21, size=1, fill="white") +
  geom_errorbar(aes(x=City, 
                    ymin=means-qnormal*standard.deviations/100, 
                    ymax=means+qnormal*standard.deviations/100),
                col="red", width=.4) +
  geom_errorbar(data=subset(ApproximateQuantities, City!="Miami"),
                aes(x=City,
                    ymin=means-qnormal*standard.deviations/100, 
                    ymax=means+qnormal*standard.deviations/100),
                col="darkgrey", width=.4) +
  ylab("Mean Entropy") +
  facet_grid(~Measurement) +
  theme_pander() +
  theme(text=element_text(family="serif"))
