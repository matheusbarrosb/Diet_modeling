#-------------------------------------------------------------------------------
# (ANALOGOUS OF) A BAYESIAN TWO-WAY ANOVA
#-------------------------------------------------------------------------------

# MODEL SPECIFICATION

AOV_model <- "
model {

# PRIORS -----------------------------------------------------------------------

for (i in 1:N_sites) {
  for (j in 1:N_treat) {
  
    group_mean[i,j] ~ dnorm(0, 0.0001) T(0,)
  
  }
}

# hyperparameters for tau
scale ~ dunif(0,1)
rate ~ dunif(0,1)

for (i in 1:N_sites) {
  for (j in 1:N_treat) {

    tau[i,j] ~ dgamma(scale, rate) # set up unequal variances
    sigma_sq[i,j] <- 1/tau[i,j]

  }
}

# LIKELIHOOD -------------------------------------------------------------------

for (i in 1:N_obs) {

  mean[i] <- group_mean[site[i], treat[i]]
  y[i] ~ dnorm(mean[i], tau[site[i], treat[i]])

}

}"

AOV_string <- textConnection(AOV_model)


# DATA WRANGLING
TL <-LAGRHO$Length
y <- (as.numeric(LAGRHO$`Gut weight`)/TL)*100
site <- as.numeric(as.factor(LAGRHO$Site))
treat <- as.numeric(as.factor(LAGRHO$Treatment))
N_sites <- length(unique(data.df$site))
N_treat = length(unique(data.df$treat))

data.df <- data.frame(y, site, treat)
data.df <- na.exclude(data.df)

data_AOV <- list(y = data.df$y, site = data.df$site, treat = data.df$treat, 
                 N_sites = N_sites, N_treat = N_treat,
                 N_obs = length(data.df$y),
                 K = N_sites*N_treat) 

ni = 5000; nb = 1000; nt = 1; nc = 3

fit.aov <- jags(data = data_AOV,
            inits = NULL,
            parameters = c("group_mean"),
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

levels(mcmc$Parameter)

mcmc.df <- mcmc %>% 
  filter(Parameter ) %>%
  droplevels()

MCMC.df <- filter(mcmc, Parameter %in% c("group_mean[1,1]", "group_mean[2,1]", "group_mean[3,1]",
                              "group_mean[5,1]", "group_mean[6,1]", "group_mean[1,2]",
                              "group_mean[4,2]", "group_mean[6,2]", "group_mean[7,2]"))

ggs_caterpillar(MCMC.df)











