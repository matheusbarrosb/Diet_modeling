#-------------------------------------------------------------------------------
# (ANALOGOUS OF) A BAYESIAN TWO-WAY ANOVA
#-------------------------------------------------------------------------------
require(rjags)
require(R2jags)
require(see)
require(ggmcmc)
require(ggplot2)

#-------------------------------------------------------------------------------
# FORMATTING DATA

diet.raw <- Seine_Data_DT

LAGRHO <- diet.raw %>%
  filter(diet.raw$`Species code` == "LAGRHO")


LAGRHO <- LAGRHO %>%
  group_by(`Fish ID`, Year, Treatment, Length, Site, `Wet Weight`,`Gut weight`) %>%
  summarize(count = n())

# MODEL SPECIFICATION

AOV_model <- "
model {

# PRIORS -----------------------------------------------------------------------

sigma_mean ~ dunif(0,100) 

for (i in 1:N_sites) {
  for (j in 1:N_treat) {
  
    group_mean[i,j] ~ dnorm(mu_mean, 1/(sigma_mean*sigma_mean))
  
  }
}

# hyperparameters for lognormal distribution for GROUP MEANS
mu_mean ~ dnorm(log_meanY, tau_mu)
sigma_mu ~ dunif(0,100)
tau_mu <- pow(sigma_mu, -2)

#-------------------------------------------------------------------------------
# now let's set hyperparameters for the likelihood function

# hyperparameters for tau
scale ~ dgamma(0.001, 0.001)
rate ~ dgamma(0.001, 0.001)

for (i in 1:N_sites) {
  for (j in 1:N_treat) {

    sigma_[i,j] ~ dgamma(scale, rate) # set up unequal variances
    tau[i,j] <- 1/(sigma_[i,j]*sigma_[i,j])

  }
}

# continuous covariate
beta_TL ~ dnorm(0, 0.001)

# LIKELIHOOD -------------------------------------------------------------------

for (i in 1:N_obs) {

  mean[i] <- group_mean[site[i], treat[i]] + beta_TL*TL[i]
  y[i] ~ dlnorm(mean[i], tau[site[i], treat[i]])

}

}"

AOV_string <- textConnection(AOV_model)


# DATA WRANGLING FOR JAGS ------------------------------------------------------
TL <-LAGRHO$Length
y <- (as.numeric(LAGRHO$`Gut weight`)/LAGRHO$`Wet Weight`)
site <- as.numeric(as.factor(LAGRHO$Site))
treat <- as.numeric(as.factor(LAGRHO$Treatment))
N_sites <- length(unique(site))
N_treat = length(unique(treat))

data.df <- data.frame(y, site, treat, TL)
data.df <- na.exclude(data.df)

data_AOV <- list(y = data.df$y, site = data.df$site,
                 treat = data.df$treat, TL = data.df$TL,
                 N_sites = N_sites, N_treat = N_treat,
                 N_obs = length(data.df$y),
                 K = N_sites*N_treat, log_meanY = log(mean(data.df$y)))

ni = 10000; nb = 1000; nt = 1; nc = 3

fit.aov <- jags(data = data_AOV,
                inits = NULL,
                parameters = c("group_mean", "mu_mean"),
                model.file = AOV_string,
                n.thin = nt,
                n.chains = nc,
                n.burnin = nb,
                n.iter = ni)

plot(fit.aov)

# PLOTTING

require(ggmcmc)
require(ggplot2)

mcmc <- ggs(as.mcmc(fit.aov))

MCMC.df <- filter(mcmc, Parameter %in% c("group_mean[1,1]", "group_mean[2,1]", "group_mean[3,1]",
                                         "group_mean[5,1]", "group_mean[6,1]",
                                         "group_mean[4,2]", "group_mean[6,2]", "group_mean[7,2]"))

ggplot(MCMC.df) +
  geom_violinhalf(aes(x = Parameter, y = exp(value)),
                  position = position_nudge(x = .2, y = 0),
                  color = NA, fill = "grey72") + ylim(0,0.25)+
  geom_jitter(data = data.df, aes(x = site, y = y), alpha = 0.35, width = 0.15) +
  theme_bw() 


