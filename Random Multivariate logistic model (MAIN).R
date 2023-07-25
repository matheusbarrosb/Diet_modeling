#-------------------------------------------------------------------------------
# RANDOM BAYESIAN MULTIVARIATE LOGISTIC MODEL FOR FISH DIET USING JAGS
#-------------------------------------------------------------------------------

# this is a bunch of code to fit a multivariate logistic model to fish diet 
# composition data. I use a shared random intercept structure to allow to
# account for the correlation structure between the multiple outcome variables

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
  X = X,                                                               
  K = K,                                                               
  N = length(amph),                                           
  re = site,                                       
  b0 = rep(0,K),
  B0 = diag(0.001, K),
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
        alpha ~ dmnorm(a0, tau.re * A0[,])
        num ~ dnorm(0, 0.0016)
        denom ~ dnorm(0, 1)

        sigma.re ~ dunif(0,10)
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

            # constraining the link function to log odds between 20 and -20
            # to avoid extreme values
            logit(p_amph[i]) <- max(-20, min(20, eta_amph[i]))
            logit(p_pol[i]) <- max(-20, min(20, eta_pol[i]))
            logit(p_SAV[i]) <- max(-20, min(20, eta_SAV[i]))
            logit(p_crust[i]) <- max(-20, min(20, eta_crust[i]))
            logit(p_isop[i]) <- max(-20, min(20, eta_isop[i]))
            logit(p_tanad[i]) <- max(-20, min(20, eta_tanad[i]))

            eta_amph[i] <- inprod(betas_amph[], X[i,]) + alpha[re[i]]
            eta_pol[i] <- inprod(betas_pol[], X[i,]) + alpha[re[i]]
            eta_SAV[i] <- inprod(betas_SAV[], X[i,]) + alpha[re[i]]
            eta_crust[i] <- inprod(betas_crust[], X[i,]) + alpha[re[i]]
            eta_isop[i] <- inprod(betas_isop[], X[i,]) + alpha[re[i]]
            eta_tanad[i] <- inprod(betas_tanad[], X[i,]) + alpha[re[i]]
        }
    }"

model.spec<-textConnection(modelString)

#parameters to be monitored
params <- c("betas_amph","betas_pol","betas_SAV", "betas_crust", "betas_isop", "betas_isop")


nt = 1; nc = 3; nb = 5000; ni = 50000

fit <- jags(data = model.data,
             inits = NULL,
             parameters = params,
             model.file = model.spec,
             n.thin = nt,
             n.chains = nc,
             n.burnin = nb,
             n.iter = ni) # working as of 20230720. Need to increase iterations

print(fit, intervals=c(0.2, 0.8), digits=2)
plot(fit)


MCMC <- as.mcmc(fit)

MCMC.df <- ggs(MCMC)
head(MCMC.df)

# PLOTTING

# 1. Amphipods

ggs_caterpillar(MCMC.df, greek = TRUE,
                family = "betas_amph") + 
  theme_custom() +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  ggtitle("Amphipods")


# 2. Polychaetes

ggs_caterpillar(MCMC.df, greek = TRUE,
                family = "betas_pol") + 
  theme_custom() +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  ggtitle("Polychaetes")


# 3. Other crustaceans

ggs_caterpillar(MCMC.df, greek = TRUE,
                family = "betas_crust") + 
  theme_custom() +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  ggtitle("Other crustaceans")

ggmcmc(MCMC.df)



# 4. SAV

ggs_caterpillar(MCMC.df, greek = TRUE,
                family = "betas_SAV") + 
  theme_custom() +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  ggtitle("SAV")



# 5. Isopods

ggs_caterpillar(MCMC.df, greek = TRUE,
                family = "betas_isop") + 
  theme_custom() +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  ggtitle("Isopods")


# 6. Tanaidacea

ggs_caterpillar(MCMC.df, greek = TRUE,
                family = "betas_tanad") + 
  theme_custom() +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  ggtitle("Tanaidacea")


# PLOTTING POSTERIOR PROBABILITIES OF ECOUNTER FOR A PREY

probs <- fit$BUGSoutput$median$p_SAV
sites <- LAGRHO$Site
treat <- LAGRHO$Treatment

plot(probs ~ as.factor(sites))

prob.df <- data.frame(probs, sites, TL)

ggplot(data = prob.df, aes(x = sites, y = probs, fill = treat)) +
  geom_jitter(width = 0.2) + 
  theme_custom() + geom_flat_violin(position = position_nudge(x = .25),
                                    alpha = 0.4) +
  xlab("") + scale_fill_manual(values = c("grey", "seagreen")) +
  ylab("Posterior probability of occurrence") +
  ggtitle("Pinfish eat more SAV at certain sites")


TL2 <- runif(1000, 30, 300)

prob_SAV_mean <- 1/(1 + exp(-(-8.16 + 0.04*TL2)))
prob_SAV_low <- 1/(1 + exp(-(-5.16 + 0.03*TL2)))
prob_SAV_high <- 1/(1 + exp(-(-12.16 + 0.05*TL2)))

log_odds_mean <- log(prob_SAV_mean/(1-prob_SAV_mean))
log_odds_low <- log(prob_SAV_low/(1-prob_SAV_low))
log_odds_high <- log(prob_SAV_high/(1-prob_SAV_high))

probs.df <- data.frame(TL2, prob_SAV_mean, prob_SAV_low, prob_SAV_high,
                      log_odds_high, log_odds_low, log_odds_mean)

ggplot(data = probs.df, aes(x = TL2, y = prob_SAV_mean)) +
  geom_line() + geom_ribbon(aes(ymin = prob_SAV_low,
                                ymax = prob_SAV_high),
                            alpha = 0.3, fill = "seagreen") +
  theme_custom() + ylab("Posterior probability of ocurrence") +
  xlab("TL (mm)") + ggtitle("SAV ingestion is more likely for bigger fish")





