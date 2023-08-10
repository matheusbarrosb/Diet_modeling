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
require(see)
#-------------------------------------------------------------------------------
# FORMATING RAW DATA  ----------------------------------------------------------

diet.raw <- Seine_Data_DT_1_

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
                        LAGRHO$Polychaete, LAGRHO$Isopod, LAGRHO$Tanaidacea,
                        LAGRHO$Fish)

names(LAGRHO.df) <- c('treat', 'site', 'TL', 'SAV', 'amphipod', 'crustacean',
                      'polychaetae', 'isopod', 'tanaidacea', 'fish')

LAGRHO.df <- LAGRHO.df %>% group_by(site, treat) %>% 
  summarise(TL = mean(TL), SAV = sum(SAV), amphipod = sum(amphipod), crustacean = sum(crustacean),
            polychaete = sum(polychaetae), tanaidacea = sum(tanaidacea), isopod = sum(isopod),
            fish = sum(fish))


# FORMATTING DATA FOR MODEL INPUTS ---------------------------------------------

y = as.matrix(LAGRHO.df[,4:10]) # response variables
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
  N_lvls_treat = N_lvls_treat,
  total_prey_obs = colSums(y),
  total_obs = sum(N_obs_site)
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

  # OVERALL PoE FOR EACH PREY K

    for (k in 1:N_prey) {

      total_prey_obs[k] ~ dbin(ov_POE[k], total_obs)

      ov_POE[k] ~ dbeta(omega_ov[k]*(kappa_ov[k]-2)+1, (1-omega_ov[k])*(kappa_ov[k]-2)+1)

      omega_ov[k] ~ dbeta(1,1)

      kappa_ov[k] ~ dgamma(0.01,0.01)

    }

  }
  " # end of model specification

model.spec <- textConnection(modelstring)

#parameters to be monitored
params <- c("p")

nt = 1; nc = 3; nb = 1000; ni = 10000

fit <- jags(data = dataList,
            inits = NULL,
            parameters = params,
            model.file = model.spec,
            n.thin = nt,
            n.chains = nc,
            n.burnin = nb,
            n.iter = ni)

plot(fit)
print(fit, intervals=c(0.2, 0.8), digits=2)

site_treat<- factor(interaction(LAGRHO.df$site, LAGRHO.df$treat, sep = "-"))

# get PoEs for each prey at each site
p_SAV <- NULL
p_amph <- NULL
p_crust <- NULL
p_pol <- NULL
p_tanad <- NULL
p_isop <- NULL
p_fish <- NULL

for (i in 1:Ntotal) {
  p_SAV[i] <- fit$BUGSoutput$mean$p[i,1]
  p_amph[i] <- fit$BUGSoutput$mean$p[i,2]
  p_crust[i] <- fit$BUGSoutput$mean$p[i,3]
  p_pol[i] <- fit$BUGSoutput$mean$p[i,4]
  p_tanad[i] <- fit$BUGSoutput$mean$p[i,5]
  p_isop[i] <- fit$BUGSoutput$mean$p[i,6]
  p_fish[i] <- fit$BUGSoutput$mean$p[i,7]
}

plot(p_SAV ~ site_treat)
plot(p_amph ~ site_treat)
plot(p_crust ~ site_treat)
plot(p_pol ~ site_treat)
plot(p_tanad ~ site_treat)
plot(p_isop ~ site_treat)
plot(p_fish ~ site_treat)

# POSTERIOR PREDICTIVE CHECKS

ob_p <- as.data.frame(dataList$y/dataList$N_obs_site) # observed PoEs

pred_p <- matrix(NA, nrow = dataList$Ntotal, ncol = dataList$N_prey)
for(i in 1:(dataList$N_prey)) {
  for (j in 1:(dataList$Ntotal)) {pred_p[j,i] <- fit$BUGSoutput$mean$p[j,i]}
} # get predicted PoEs

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
abline(lm(observed ~ predicted), col = "blue")
abline(0,1, col = "red")
legend(x = "topleft",          
       legend = c("Expected", "Predicted"),
       col = c("blue", "red"), 
       lwd = 2)  


# FORMATTING GGS DATAFRAME 

mcmc.df <- ggs(as.mcmc(fit))

mcmc_df <- mcmc.df %>%
  mutate(Prey = case_when(Parameter == 'p[1,1]' | Parameter == 'p[2,1]' | Parameter == 'p[3,1]' |
                          Parameter == 'p[4,1]' | Parameter == 'p[5,1]' | Parameter == 'p[6,1]' |
                          Parameter == 'p[7,1]' ~ 'SAV',
                          Parameter == 'p[1,2]' | Parameter == 'p[2,2]' | Parameter == 'p[3,2]' |
                          Parameter == 'p[4,2]' | Parameter == 'p[5,2]' | Parameter == 'p[6,2]' |
                          Parameter == 'p[7,2]' ~ 'Amphipod',
                          Parameter == 'p[1,3]' | Parameter == 'p[2,3]' | Parameter == 'p[3,3]' |
                          Parameter == 'p[4,3]' | Parameter == 'p[5,3]' | Parameter == 'p[6,3]' |
                          Parameter == 'p[7,3]' ~ 'Crustacean',
                          Parameter == 'p[1,4]' | Parameter == 'p[2,4]' | Parameter == 'p[3,4]' |
                          Parameter == 'p[4,4]' | Parameter == 'p[5,4]' | Parameter == 'p[6,4]' |
                          Parameter == 'p[7,4]' ~ 'Polychaete',
                          Parameter == 'p[1,5]' | Parameter == 'p[2,5]' | Parameter == 'p[3,5]' |
                          Parameter == 'p[4,5]' | Parameter == 'p[5,5]' | Parameter == 'p[6,5]' |
                          Parameter == 'p[7,5]' ~ 'Tanaidacea',
                          Parameter == 'p[1,6]' | Parameter == 'p[2,6]' | Parameter == 'p[3,6]' |
                          Parameter == 'p[4,6]' | Parameter == 'p[5,6]' | Parameter == 'p[6,6]' |
                          Parameter == 'p[7,6]' ~ 'Isopod',
                          Parameter == 'p[1,7]' | Parameter == 'p[2,7]' | Parameter == 'p[3,7]' |
                          Parameter == 'p[4,7]' | Parameter == 'p[5,7]' | Parameter == 'p[6,7]' |
                          Parameter == 'p[7,7]' ~ 'Fish'))

mcmc_df <- mcmc_df %>%
  mutate(Site = case_when(Parameter == 'p[1,1]' | Parameter == 'p[1,2]' | Parameter == 'p[1,3]' |
                          Parameter == 'p[1,4]' | Parameter == 'p[1,5]' | Parameter == 'p[1,6]' |
                          Parameter == 'p[1,7]' ~ 'AM-C',
                          Parameter == 'p[2,1]' | Parameter == 'p[2,2]' | Parameter == 'p[2,3]' |
                          Parameter == 'p[2,4]' | Parameter == 'p[2,5]' | Parameter == 'p[2,6]' |
                          Parameter == 'p[2,7]' ~ 'DR-C',
                          Parameter == 'p[3,1]' | Parameter == 'p[3,2]' | Parameter == 'p[3,3]' |
                          Parameter == 'p[3,4]' | Parameter == 'p[3,5]' | Parameter == 'p[3,6]' |
                          Parameter == 'p[3,7]' ~ 'HWP-T',
                          Parameter == 'p[4,1]' | Parameter == 'p[4,2]' | Parameter == 'p[4,3]' |
                          Parameter == 'p[4,4]' | Parameter == 'p[4,5]' | Parameter == 'p[4,6]' |
                          Parameter == 'p[4,7]' ~ 'PaP-C',
                          Parameter == 'p[5,1]' | Parameter == 'p[5,2]' | Parameter == 'p[5,3]' |
                          Parameter == 'p[5,4]' | Parameter == 'p[5,5]' | Parameter == 'p[5,6]' |
                          Parameter == 'p[5,7]' ~ 'LB-C',
                          Parameter == 'p[6,1]' | Parameter == 'p[6,2]' | Parameter == 'p[6,3]' |
                          Parameter == 'p[6,4]' | Parameter == 'p[6,5]' | Parameter == 'p[6,6]' |
                          Parameter == 'p[6,7]' ~ 'PaP-T',
                          Parameter == 'p[7,1]' | Parameter == 'p[7,2]' | Parameter == 'p[7,3]' |
                          Parameter == 'p[7,4]' | Parameter == 'p[7,5]' | Parameter == 'p[7,6]' |
                          Parameter == 'p[7,7]' ~ 'SA-T'))

mcmc_df <- na.exclude(mcmc_df)


ggplot(data = mcmc_df, aes(y = value, x = Site)) +
  geom_violinhalf(trim = TRUE, scale = "width",
                  fill = "seagreen", alpha = 0.3) +
  facet_wrap(~Prey) + 
  theme_custom()



