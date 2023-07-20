#-------------------------------------------------------------------------------
# RANDOM BAYESIAN MULTIVARIATE LOGISTIC MODEL FOR FISH DIET USING JAGS
#-------------------------------------------------------------------------------
require(rstanarm)
require(dplyr)
require(tidyverse)
require(rjags)
require(R2jags)
require(rstan)

#-------------------------------------------------------------------------------
# formatting data

diet.raw <- Seine_Data_DT

LAGRHO <- diet.raw %>%
  filter(diet.raw$`Species code` == "LAGRHO")


LAGRHO <- LAGRHO %>%
  group_by(`Fish ID`, Year, Treatment, Group, Length, Site) %>%
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


X <- model.matrix(~ site + treat + TL)
K <- ncol(X)
re <- as.numeric(site)
Nre <- length(unique(site))

model.data <- list(
  amph = amph,  
  pol = pol,
  SAV = SAV,
  crust = crust,
  isop = isop,
  tanad = tanad,
  X = X,                                                               # Covariates
  K = K,                                                               # Num. betas
  N = length(amph),                                              # Sample size
  re = site,                                           # Random effects
  b0 = rep(0,K),
  B0 = diag(0.0001, K),
  a0 = rep(0,Nre),
  A0 = diag(Nre))


modelString <- "
    model {
        #-----------------------------------------------------------------------
        #                            FIXED EFFECTS
        #-----------------------------------------------------------------------
        betas_amph ~ dmnorm(b0[], B0[,])
        betas_pol ~ dmnorm(b0[], B0[,])
        betas_SAV ~ dmnorm(b0[], B0[,])
        betas_crust ~ dmnorm(b0[], B0[,])
        betas_isop ~ dmnorm(b0[], B0[,])
        betas_tanad ~ dmnorm(b0[], B0[,])

        #-----------------------------------------------------------------------
        #                            RANDOM EFFECTS
        #-----------------------------------------------------------------------
        # Priors for random effect group
        a ~ dmnorm(a0, tau.re * A0[,])
        num ~ dnorm(0, 0.0016)
        denom ~ dnorm(0, 1)

        sigma.re ~ dunif(0,100)
        tau.re <- 1 / (sigma.re * sigma.re)

        #-----------------------------------------------------------------------
        #                              LIKELIHOOD
        #-----------------------------------------------------------------------
        for (i in 1:N) {
            amph[i] ~ dbern(p_amph[i])
            pol[i] ~ dbern(p_pol[i])
            SAV[i] ~ dbern(p_SAV[i])
            crust[i] ~ dbern(p_crust[i])
            isop[i] ~ dbern(p_isop[i])
            tanad[i] ~ dbern(p_tanad[i])

            logit(p_amph[i]) <- max(-20, min(20, eta_amph[i]))
            logit(p_pol[i]) <- max(-20, min(20, eta_pol[i]))
            logit(p_SAV[i]) <- max(-20, min(20, eta_SAV[i]))
            logit(p_crust[i]) <- max(-20, min(20, eta_crust[i]))
            logit(p_isop[i]) <- max(-20, min(20, eta_isop[i]))
            logit(p_tanad[i]) <- max(-20, min(20, eta_tanad[i]))

            eta_amph[i] <- inprod(betas_amph[], X[i,]) + a[re[i]]
            eta_pol[i] <- inprod(betas_pol[], X[i,]) + a[re[i]]
            eta_SAV[i] <- inprod(betas_SAV[], X[i,]) + a[re[i]]
            eta_crust[i] <- inprod(betas_crust[], X[i,]) + a[re[i]]
            eta_isop[i] <- inprod(betas_isop[], X[i,]) + a[re[i]]
            eta_tanad[i] <- inprod(betas_tanad[], X[i,]) + a[re[i]]
        }
    }"

model.spec<-textConnection(modelString)

#parameters to be monitored
params <- c("betas_amph", "betas_pol", "betas_SAV", "a", "sigma.re", "tau.re", "p_amph", "p_pol")


fit <- jags(data = model.data,
             inits = NULL,
             parameters = params,
             model.file = model.spec,
             n.thin = 1,
             n.chains = 3,
             n.burnin = 5000,
             n.iter = 15000) # working as of 20230720. Need to increase iterations

print(fit, intervals=c(0.025, 0.975), digits=3)
plot(fit)








