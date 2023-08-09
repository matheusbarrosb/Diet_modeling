#-------------------------------------------------------------------------------
# MULTIVARIATE BETA-BINOMIAL MODEL FOR FISH DIET USING JAGS
#-------------------------------------------------------------------------------

# this is a bunch of code to fit a multivariate beta-binomial model to estimate
# probabilities of encounter of prey for different sites
require(dplyr)
require(tidyverse)
require(rjags)
require(R2jags)
require(ggmcmc)
#-------------------------------------------------------------------------------
# FORMATING RAW DATA  ----------------------------------------------------------

diet.raw <- Seine_Data_DT

LAGRHO <- diet.raw %>%
  filter(diet.raw$`Species code` == "LAGRHO")


LAGRHO <- LAGRHO %>%
  group_by(`Fish ID`,Year, Treatment, Group, Length, Site) %>%
  summarize(count = n())

LAGRHO <- LAGRHO %>%
  filter(Site %in% c("AM", "DR", "HWP", "LB", "NEPaP", "SA"))

LAGRHO <- LAGRHO %>%
  pivot_wider(names_from = Group, values_from = count, values_fill = 0)


LAGRHO.df <- data.frame(LAGRHO$Treatment, LAGRHO$Site, LAGRHO$Length,
                        LAGRHO$SAV, LAGRHO$Amphipod, LAGRHO$Crustacean,
                        LAGRHO$Polychaete, LAGRHO$Tanaidacea)

names(LAGRHO.df) <- c('treat', 'site', 'TL', 'SAV', 'amphipod', 'crustacean',
                      'polychaetae', 'tanaidacea')

LAGRHO.df <- LAGRHO.df %>% group_by(site, treat) %>% 
  summarise(TL = mean(TL), SAV = sum(SAV), amphipod = sum(amphipod), crustacean = sum(crustacean),
            polychaete = sum(polychaetae), tanaidacea = sum(tanaidacea))


# FORMATTING DATA FOR MODEL INPUTS ---------------------------------------------

y = as.matrix(LAGRHO.df[,4:8]) # response variables
site = as.numeric(as.factor(LAGRHO.df$site)) # categorical covariate 01
treat = as.numeric(as.factor(LAGRHO.df$treat)) 

# number of observations ('visits') for each site for binomial distribution
N_obs_site1 <- LAGRHO %>% 
  group_by(Site, Treatment) %>%
  summarise(no_rows = length(Site))
N_obs_site <- N_obs_site1$no_rows

# constrain prey observations to <= No of observations per site for binomial distribution
y <- ifelse(y > N_obs_site, (N_obs_site-1), y)

sitelevels = levels(as.factor(LAGRHO.df$site)) # levels for categorical covariate
treatlevels = levels(as.factor(LAGRHO.df$treat)) 
Ntotal = length(LAGRHO.df$site) # number of sites
N_lvls_site = length(unique(site)) # number of levels for categorical covariate
N_lvls_treat = length(unique(treat))

# Specify the data in a list for sending to JAGS:
dataList = list(
  y = y,
  TL = LAGRHO.df$TL,
  N_obs_site = N_obs_site,
  site = site,
  treat = treat,
  N_prey = ncol(y),
  Ntotal = Ntotal,
  N_lvls_site = N_lvls_site,
  N_lvls_treat = N_lvls_treat
)


modelstring = "
  model {
  
    for (k in 1:N_prey) {
      for (i in 1:Ntotal) {
    
        # OBSERVATION MODEL
        y[i,k] ~ dbin(p[i,k], N_obs_site[i])
  
        # PROCESS MODEL
        p[i,k] ~ dbeta(omega[site[i], treat[i]]*(kappa[k]-2)+1, (1-omega[site[i], treat[i]])*(kappa[k]-2)+1)

      }
    }

    for (k in 1:N_lvls_treat) {
      for (j in 1:N_lvls_site) {
    
        omega[j,k] <- ilogit(a0 + a[j,k] + beta_TL * TL[j])
      
        a[j,k] ~ dnorm(0, 1/(aSigma*aSigma))
      
      }

    }

    a0 ~ dnorm(0.0 , 1/2^2) 
    aSigma ~ dgamma(0.01, 0.01)  
    beta_TL ~ dnorm(0, 0.001)

    for (k in 1:N_prey) {

    kappaMinusTwo[k] ~ dgamma(0.01 , 0.01)  # varying K per prey group
    kappa[k] <- kappaMinusTwo[k] + 2

    }    

    # Convert a0,a[] to sum-to-zero b0,b[] :
    
    for (k in 1:N_lvls_treat) {
      for (j in 1:N_lvls_site) { 

        m[j,k] <- a0 + a[j,k] 
      
      }  
    }
    
    b0 <- mean(m)

    for (k in 1:N_lvls_treat) {
      for (j in 1:N_lvls_site) {

        b[j,k] <- m[j,k] - b0

      }
    }

  }
  " # end of model specification

model.spec <- textConnection(modelstring)

#parameters to be monitored
params <- c("p")

nt = 1; nc = 3; nb = 1000; ni = 5000

fit <- jags(data = dataList,
            inits = NULL,
            parameters = params,
            model.file = model.spec,
            n.thin = nt,
            n.chains = nc,
            n.burnin = nb,
            n.iter = ni)

print(fit, intervals=c(0.2, 0.8), digits=2)

mcmc.df <- ggs(as.mcmc(fit))


ggs_caterpillar(mcmc.df) + theme_nice() + xlim(0,1)



site_treat<- factor(interaction(LAGRHO.df$site, LAGRHO.df$treat, sep = "-"))

p_SAV <- NULL
for (i in 1:Ntotal) {
  p_SAV[i] <- fit$BUGSoutput$mean$p[i,1]
}


plot(p_SAV ~ site_treat)

# POSTERIOR PREDICTIVE CHECKS

ob_p <- as.data.frame(dataList$y/dataList$N_obs_site) 

pred_p <- matrix(NA, nrow = dataList$Ntotal, ncol = dataList$N_prey)
for(i in 1:(dataList$N_prey)) {
  for (j in 1:(dataList$Ntotal)) {pred_p[j,i] <- fit$BUGSoutput$mean$p[j,i]}
}

pred_p <- as.data.frame(pred_p)
names(pred_p) <- c("SAV", "amphipod", "crustacean", "polychaete", "tanaidacea")

SAV_pred <- pred_p$SAV
amph_pred <- pred_p$amphipod
crust_pred <- pred_p$crustacean
pol_pred <- pred_p$polychaete
tanad_pred <- pred_p$tanaidacea

predicted <- c(SAV_pred, amph_pred, crust_pred, pol_pred, tanad_pred)

SAV_obs <- ob_p$SAV
amph_obs <- ob_p$amphipod
crust_obs <- ob_p$crustacean
pol_obs <- ob_p$polychaete
tanad_obs <- ob_p$tanaidacea

observed <- c(SAV_obs, amph_obs, crust_obs, pol_obs, tanad_obs)

plot(observed ~ predicted)



