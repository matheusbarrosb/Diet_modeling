#-------------------------------------------------------------------------------
# RANDOM BAYESIAN MULTIVARIATE LOGISTIC MODEL FOR FISH DIET USING JAGS
#-------------------------------------------------------------------------------

# this is a bunch of code to fit a multivariate logistic model to fish diet 
# composition data. I use a shared random intercept structure to allow to
# account for the correlation structure between the multiple outcome variables

# OBS.:
# this version uses a matrix as a single response variable instead of fitting
# multiple univariate regressions and is more computationally efficient
require(dplyr)
require(tidyverse)
require(rjags)
require(R2jags)
require(ggmcmc)
#-------------------------------------------------------------------------------
# formatting data

diet.raw <- Seine_Data_DT

LAGRHO <- diet.raw %>%
  filter(diet.raw$`Species code` == "LAGRHO")


LAGRHO <- LAGRHO %>%
  group_by(`Fish ID`, Year, Treatment, Group, Length, Site, `Gut weight`) %>%
  summarize(count = n())

LAGRHO <- LAGRHO %>%
  pivot_wider(names_from = Group, values_from = count, values_fill = 0)


LAGRHO.df <- data.frame(LAGRHO$Treatment, LAGRHO$Site, LAGRHO$Length,
                        LAGRHO$SAV, LAGRHO$Amphipod, LAGRHO$Crustacean,
                        LAGRHO$Polychaete, LAGRHO$Tanaidacea)

names(LAGRHO.df) <- c('treat', 'site', 'TL', 'SAV', 'amphipod', 'crustacean',
                      'polychaetae', 'tanaidacea')

LAGRHO.df <- LAGRHO.df %>%
  pivot_longer(names_to = "prey", values_to = "PA", cols = c('SAV', 'amphipod', 'crustacean', 'polychaetae', 'tanaidacea'))

# -----------------------------------------------------------------------------
# transforming variables for JAGS
site <- as.factor(LAGRHO$Site)
site <- as.numeric(site)

treat <- as.factor(LAGRHO$Treatment)
treat <- as.numeric(treat)

TL <- LAGRHO$Length

amph <- as.numeric(LAGRHO$Amphipod)
pol <- as.numeric(LAGRHO$Polychaete)
SAV <- as.numeric(LAGRHO$SAV)
crust <- as.numeric(LAGRHO$Crustacean)
isop <- as.numeric(LAGRHO$Isopod)
tanad <- as.numeric(LAGRHO$Tanaidacea)

amph <- ifelse(amph > 1, 1, amph)
pol <- ifelse(pol > 1, 1, pol)
SAV <- ifelse(SAV > 1, 1, SAV)
crust <- ifelse(crust > 1, 1, crust)
isop <- ifelse(isop > 1, 1, isop)
tanad <- ifelse(tanad > 1, 1, tanad)

prey.df <- data.frame(amph, pol, SAV, crust, isop, tanad)
prey.mat <- as.matrix(prey.df)

X <- model.matrix(prey.mat~ site + treat + TL)
K <- ncol(X)
re <- as.numeric(site)
N_sites <- length(unique(site))
prey <- as.factor(seq(1,ncol(prey.mat)))
Nre <- length(unique(prey))

model.data <- list(
  y = prey.mat,
  z = prey.mat,
  TL = TL,
  site = as.factor(LAGRHO$Site),
  treat = as.factor(treat),
  N_obs = length(amph),
  N_prey = ncol(prey.mat),
  N_sites = N_sites,
  re = prey,
  a0 = rep(0,Nre),
  A0 = diag(Nre))


modelString <- "model

{
  
  # LIKELIHOOD ----------------------------------------------------------------- 
  
  for(k in 1:N_prey){

    for(i in 1:N_obs) {

      logit(psi[i, k]) <- alpha[re[k]] + beta_TL[k]*TL[i] + beta_site[k]*site[i] + beta_treat[k]*treat[i]
  
      z[i, k] ~ dbern(psi[i, k])
  
    }  
  
    # PRIORS (PREY LEVEL) ------------------------------------------------------
    
    # FIXED INTERCEPT - AVERAGE LOG ODDS OF ENCOUNTER FOR PREY K
    b0[k] ~ dnorm(mu_b0, tau_b0)
  
    mu_eta[k] <- mu_lp + rho*sigma_lp/sigma_b0*(b0[k] - mu_b0)
    lp[k] ~ dnorm(mu_eta[k], tau_eta)
    p[k] <- ilogit(lp[k]) # back to probability scale
  
    beta_TL[k] ~ dnorm(0, 0.001)
    beta_site[k] ~ dnorm(0, 0.001)
    beta_treat[k] ~ dnorm(0, 0.001)
    
  }
  
  # HYPERPRIORS ----------------------------------------------------------------
  
  alpha ~ dmnorm(a0, A0[,]) # random intercept follows a MVN to account for correlations between prey
  
  #FIXED INTERCEPT CONFIG.
  b0_mean ~ dbeta(1, 1) # mean fixed intercept on probability scale
  mu_b0 <- logit(b0_mean) # mean fixed intercept on logit scale
  sigma_b0 ~ dgamma(3, 1) # sd of fixed intercept
  tau_b0 <- 1/(sigma_b0^2) # precision
  
  #PRIORS FOR PREY-SPECIFIC PROBABILITIES
  p_mean ~ dbeta(1, 1) # probability scale
  mu_lp <- logit(p_mean) # logit scale
  sigma_lp ~ dunif(0, 5) # sd
  tau_lp <- 1/(sigma_lp^2) # precision
  
  rho ~ dunif(-1, 1)
  tau_eta <- tau_lp/(1 - rho^2)
  
  # DERIVED QUANTITIES
  
  for (k in 1:N_prey) {
    PoE[k] <- ilogit(b0[k])
  }
  
}"

model.spec<-textConnection(modelString)

#parameters to be monitored
params <- c("beta_TL", "beta_site", "beta_treat", "alpha", "PoE")

nt = 1; nc = 3; nb = 1000; ni = 5000

fit <- jags(data = model.data,
            inits = NULL,
            parameters = params,
            model.file = model.spec,
            n.thin = nt,
            n.chains = nc,
            n.burnin = nb,
            n.iter = ni) 

print(fit, intervals=c(0.2, 0.8), digits=2)
plot(fit)



