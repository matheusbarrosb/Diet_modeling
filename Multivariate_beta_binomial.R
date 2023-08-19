#.-------------------------------------------------------------------------------
# MULTIVARIATE BETA-BINOMIAL MODEL FOR FISH DIET USING JAGS
#.-------------------------------------------------------------------------------

# this is a bunch of code to fit a multivariate beta-binomial model to estimate
# probabilities of encounter of prey for different sites
require(dplyr)
require(tidyverse)
require(R2jags)
require(ggmcmc)
require(see)
require(flextable)

# PINFISH ----------------------------------------------------------------------
## FORMATING RAW DATA ----------------------------------------------------------
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


## FORMATTING DATA FOR MODEL INPUTS --------------------------------------------

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

## WRITING MODEL ---------------------------------------------------------------

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
params <- c("ov_POE")

nt = 1; nc = 3; nb = 1000; ni = 10000

fit_MBB_lAGRHO <- jags(data = dataList,
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
mean_p_SAV <- NULL;mean_p_amph <- NULL;mean_p_crust <- NULL;mean_p_pol <- NULL
mean_p_tanad <- NULL;mean_p_isop <- NULL;mean_p_fish <- NULL

up_p_SAV <- NULL;up_p_amph <- NULL;up_p_crust <- NULL;up_p_pol <- NULL
up_p_tanad <- NULL;up_p_isop <- NULL;up_p_fish <- NULL

low_p_SAV <- NULL;low_p_amph <- NULL;low_p_crust <- NULL;low_p_pol <- NULL
low_p_tanad <- NULL;low_p_isop <- NULL;low_p_fish <- NULL

for (i in 1:Ntotal) {
  mean_p_SAV[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1]
  mean_p_amph[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,2]
  mean_p_crust[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,3]
  mean_p_pol[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,4]
  mean_p_tanad[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,5]
  mean_p_isop[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,6]
  mean_p_fish[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,7]
  
  up_p_SAV[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] + 2*(fit_MBB_lAGRHO$BUGSoutput$sd$p[i,1]) 
  up_p_amph[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] + 2*(fit_MBB_lAGRHO$BUGSoutput$sd$p[i,2])
  up_p_crust[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] + 2*(fit_MBB_lAGRHO$BUGSoutput$sd$p[i,3])
  up_p_pol[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] + 2*(fit_MBB_lAGRHO$BUGSoutput$sd$p[i,3])
  up_p_tanad[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] + 2*(fit_MBB_lAGRHO$BUGSoutput$sd$p[i,4])
  up_p_isop[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] + 2*(fit_MBB_lAGRHO$BUGSoutput$sd$p[i,5])
  up_p_fish[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] + 2*(fit_MBB_lAGRHO$BUGSoutput$sd $p[i,6])
  
  low_p_SAV[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] - 2*(fit_MBB_lAGRHO$BUGSoutput$sd$p[i,1]) 
  low_p_amph[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] - 2*(fit_MBB_lAGRHO$BUGSoutput$sd$p[i,2])
  low_p_crust[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] - 2*(fit_MBB_lAGRHO$BUGSoutput$sd$p[i,3])
  low_p_pol[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] - 2*(fit_MBB_lAGRHO$BUGSoutput$sd$p[i,3])
  low_p_tanad[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] - 2*(fit_MBB_lAGRHO$BUGSoutput$sd$p[i,4])
  low_p_isop[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] - 2*(fit_MBB_lAGRHO$BUGSoutput$sd$p[i,5])
  low_p_fish[i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[i,1] - 2*(fit_MBB_lAGRHO$BUGSoutput$sd $p[i,6])
}

plot(p_SAV ~ site_treat)
plot(p_amph ~ site_treat)
plot(p_crust ~ site_treat)
plot(p_pol ~ site_treat)
plot(p_tanad ~ site_treat)
plot(p_isop ~ site_treat)
plot(p_fish ~ site_treat)

## POSTERIOR PREDICTIVE CHECKS -------------------------------------------------

ob_p <- as.data.frame(dataList$y/dataList$N_obs_site) # observed PoEs

# get predicted PoEs
pred_p <- matrix(NA, nrow = dataList$Ntotal, ncol = dataList$N_prey)
pred_p_up <- matrix(NA, nrow = dataList$Ntotal, ncol = dataList$N_prey)
pred_p_low <- matrix(NA, nrow = dataList$Ntotal, ncol = dataList$N_prey)
for(i in 1:(dataList$N_prey)) {
  for (j in 1:(dataList$Ntotal)) {
    pred_p[j,i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[j,i] # means
    pred_p_up[j,i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[j,i] + (2*fit_MBB_lAGRHO$BUGSoutput$sd$p[j,i])
    pred_p_low[j,i] <- fit_MBB_lAGRHO$BUGSoutput$mean$p[j,i] - (2*fit_MBB_lAGRHO$BUGSoutput$sd$p[j,i])
    }
}

pred_p <- as.data.frame(pred_p)
pred_p_up <- as.data.frame(pred_p_up)
pred_p_low <- as.data.frame(pred_p_low)

names(pred_p) <- c("SAV", "amphipod", "crustacean", "polychaete", "tanaidacea", "isopod", "fish")
names(pred_p_up) <- c("SAV", "amphipod", "crustacean", "polychaete", "tanaidacea", "isopod", "fish")
names(pred_p_low) <- c("SAV", "amphipod", "crustacean", "polychaete", "tanaidacea", "isopod", "fish")


#means
SAV_pred_mean <- pred_p$SAV;amph_pred_mean <- pred_p$amphipod;crust_pred_mean <- pred_p$crustacean
pol_pred_mean <- pred_p$polychaete;tanad_pred_mean <- pred_p$tanaidacea; isop_pred_mean <- pred_p$isopod; fish_pred_mean <- pred_p$fish

#ups
SAV_pred_up <- pred_p_up$SAV; amph_pred_up <- pred_p_up$amphipod; crust_pred_up <- pred_p_up$crustacean
pol_pred_up <- pred_p_up$polychaete; tanad_pred_up <- pred_p_up$tanaidacea; isop_pred_up <- pred_p_up$isopod; fish_pred_up <- pred_p_up$fish

#lows
SAV_pred_low <- pred_p_low$SAV; amph_pred_low <- pred_p_low$amphipod; crust_pred_low <- pred_p_low$crustacean
pol_pred_low <- pred_p_low$polychaete; tanad_pred_low <- pred_p_low$tanaidacea; isop_pred_low <- pred_p_low$isopod; fish_pred_low <- pred_p_low$fish

predicted_mean <- c(SAV_pred_mean, amph_pred_mean, crust_pred_mean, pol_pred_mean, tanad_pred_mean, isop_pred_mean, fish_pred_mean)
predicted_up <- c(SAV_pred_up, amph_pred_up, crust_pred_up, pol_pred_up, tanad_pred_up, isop_pred_up, fish_pred_up)
predicted_low <- c(SAV_pred_low, amph_pred_low, crust_pred_low, pol_pred_low, tanad_pred_low, isop_pred_low, fish_pred_low)

SAV_obs <- ob_p$SAV
amph_obs <- ob_p$amphipod
crust_obs <- ob_p$crustacean
pol_obs <- ob_p$polychaete
tanad_obs <- ob_p$tanaidacea
isop_obs <- ob_p$isopod
fish_obs <- ob_p$fish

observed <- c(SAV_obs, amph_obs, crust_obs, pol_obs, tanad_obs, isop_obs, fish_obs)

plot(observed ~ predicted_mean)
abline(lm(observed ~ predicted_mean), col = "blue")
abline(lm(observed ~ predicted_low), col = "blue", lty = "dashed")
abline(lm(observed ~ predicted_up), col = "blue", lty = "dashed")
abline(0,1, col = "red")
legend(x = "topleft",          
       legend = c("Expected", "Predicted"),
       col = c("blue", "red"), 
       lwd = 2)  

ppcheck_LAGRHO <- data.frame(predicted_mean, predicted_low, predicted_up, observed)

## FORMATTING GGS DATAFRAME ----------------------------------------------------

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

mcmc_sum <- mcmc_df %>% group_by(Prey, Site) %>%
  summarise(median = median(value), 
            low = quantile(value, probs = c(0.05)),
            upp = quantile(value, probs = c(0.95)))


ggplot() +
  geom_violinhalf(data = mcmc_df, aes(y = value, x = Site),
                  trim = TRUE, scale = "width",
                  fill = "lightgreen", alpha = 0.3,
                  position = "dodge") +
  facet_wrap(~Prey) + 
  theme_custom() + geom_errorbar(data = mcmc_sum, 
                aes(ymin = mcmc_sum$low,
                    ymax = mcmc_sum$upp, x = Site),
                width = 0, linewidth = 0.8) +
  geom_point(data = mcmc_sum, aes(x = Site, y = median)) +
  ylab("Posterior PoE") + xlab("")


## PLOTTING RAW FREQUENCY TABLES -----------------------------------------------
raw_freq_LAGRHO <- matrix(NA, nrow = dataList$Ntotal, ncol = dataList$N_prey)


for(i in 1:(nrow(raw_freq_LAGRHO))) {
  for (j in 1:(ncol(raw_freq_LAGRHO))) {
    raw_freq_LAGRHO[j,i] <- (dataList$y[j,i]/dataList$N_obs_site[j])*100
    }}

colnames(raw_freq_LAGRHO) <- c("SAV", "Amphipods", "Other Crustaceans", "Polychaete",
                               "Tanaidacea", "Isopods", "Fish")

raw_freq_LAGRHO <- round(raw_freq_LAGRHO, 0)
raw_freq_LAGRHO <- data.frame(cbind(site_treat, as.data.frame(raw_freq_LAGRHO)))

names(raw_freq_LAGRHO) <- c("Site", "SAV", "Amphipods", "Other Crustaceans", "Polychaete",
                               "Tanaidacea", "Isopods", "Fish")

LAGRHO_rawtable <- flextable(
  data = raw_freq_LAGRHO, 
  col_keys = c("Site", "SAV", "Amphipods", "Other Crustaceans", "Polychaete",
               "Tanaidacea", "Isopods", "Fish")) %>%
  colformat_num(suffix = "%") %>%
  bold(part = "header", bold = TRUE) %>%
  align(align = "center", part = "all")%>%
  set_caption(caption = "Pinfish",
              align_with_table = TRUE
  )


## PLOTTING OVERALL PREY PoEs---------------------------------------------------



# CROAKER ----------------------------------------------------------------------
## FORMATING RAW DATA ----------------------------------------------------------
diet.raw <- Seine_Data_DT_1_

MICUND <- diet.raw %>%
  filter(diet.raw$`Species code` == "MICUND")


MICUND <- MICUND %>%
  group_by(`Fish ID`,Year, Treatment, Group, Length, Site) %>%
  summarize(count = n())

MICUND <- MICUND %>%
  filter(Site %in% c("CI", "DR", "HB", "HWP", "NEPaP", "SA"))

MICUND <- MICUND %>%
  pivot_wider(names_from = Group, values_from = count, values_fill = 0)


MICUND.df <- data.frame(MICUND$Treatment, MICUND$Site, MICUND$Length, 
                        MICUND$Amphipod, MICUND$Crustacean,
                        MICUND$Polychaete, MICUND$Isopod, MICUND$Tanaidacea,
                        MICUND$Fish)

names(MICUND.df) <- c('treat', 'site', 'TL', 'amphipod', 'crustacean',
                      'polychaetae', 'isopod', 'tanaidacea', 'fish')

MICUND.df <- MICUND.df %>% group_by(site, treat) %>% 
  summarise(TL = mean(TL, na.rm = TRUE), amphipod = sum(amphipod), crustacean = sum(crustacean),
            polychaete = sum(polychaetae), tanaidacea = sum(tanaidacea), isopod = sum(isopod),
            fish = sum(fish))


## FORMATTING DATA FOR MODEL INPUTS --------------------------------------------

y = as.matrix(MICUND.df[,4:9]) # response variables
site = as.numeric(as.factor(MICUND.df$site)) # categorical covariate 01
treat = as.numeric(as.factor(MICUND.df$treat)) 

# number of observations ('visits') for each site for binomial distribution
N_obs_site1 <- MICUND %>% 
  group_by(Site, Treatment) %>%
  summarise(no_rows = length(Site))
N_obs_site <- N_obs_site1$no_rows

# constrain prey observations to <= No of observations per site for binomial distribution
y <- ifelse(y > N_obs_site, (N_obs_site-1), y)

sitelevels = levels(as.factor(MICUND.df$site)) # levels for categorical covariate
treatlevels = levels(as.factor(MICUND.df$treat)) 
Ntotal = length(MICUND.df$site) # number of sites
N_lvls_site = length(unique(site)) # number of levels for categorical covariate
N_lvls_treat = length(unique(treat))

# Specify the data in a list for sending to JAGS:
dataList_MICUND = list(
  y = y,
  TL = MICUND.df$TL,
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

## WRITING MODEL ---------------------------------------------------------------

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
params <- c("p", "ov_POE")

nt = 1; nc = 3; nb = 1000; ni = 10000

fit_MBB_MICUND <- jags(data = dataList_MICUND,
            inits = NULL,
            parameters = params,
            model.file = model.spec,
            n.thin = nt,
            n.chains = nc,
            n.burnin = nb,
            n.iter = ni)

plot(fit_MBB_MICUND)
print(fit, intervals=c(0.2, 0.8), digits=2)

## POSTERIOR PREDICTIVE CHECKS -------------------------------------------------

ob_p <- as.data.frame(dataList_MICUND$y/dataList_MICUND$N_obs_site) # observed PoEs

# get predicted PoEs
pred_p <- matrix(NA, nrow = dataList_MICUND$Ntotal, ncol = dataList_MICUND$N_prey)
pred_p_up <- matrix(NA, nrow = dataList_MICUND$Ntotal, ncol = dataList_MICUND$N_prey)
pred_p_low <- matrix(NA, nrow = dataList_MICUND$Ntotal, ncol = dataList_MICUND$N_prey)
for(i in 1:(dataList_MICUND$N_prey)) {
  for (j in 1:(dataList_MICUND$Ntotal)) {
    pred_p[j,i] <- fit_MBB_MICUND$BUGSoutput$mean$p[j,i] # means
    pred_p_up[j,i] <- fit_MBB_MICUND$BUGSoutput$mean$p[j,i] + (2*fit_MBB_MICUND$BUGSoutput$sd$p[j,i])
    pred_p_low[j,i] <- fit_MBB_MICUND$BUGSoutput$mean$p[j,i] - (2*fit_MBB_MICUND$BUGSoutput$sd$p[j,i])
  }
}

pred_p <- as.data.frame(pred_p)
pred_p_up <- as.data.frame(pred_p_up)
pred_p_low <- as.data.frame(pred_p_low)

names(pred_p) <- c("amphipod", "crustacean", "polychaete", "tanaidacea", "isopod", "fish")
names(pred_p_up) <- c("amphipod", "crustacean", "polychaete", "tanaidacea", "isopod", "fish")
names(pred_p_low) <- c("amphipod", "crustacean", "polychaete", "tanaidacea", "isopod", "fish")

#means
amph_pred_mean <- pred_p$amphipod;crust_pred_mean <- pred_p$crustacean
pol_pred_mean <- pred_p$polychaete;tanad_pred_mean <- pred_p$tanaidacea; isop_pred_mean <- pred_p$isopod; fish_pred_mean <- pred_p$fish

#ups
amph_pred_up <- pred_p_up$amphipod; crust_pred_up <- pred_p_up$crustacean
pol_pred_up <- pred_p_up$polychaete; tanad_pred_up <- pred_p_up$tanaidacea; isop_pred_up <- pred_p_up$isopod; fish_pred_up <- pred_p_up$fish

#lows
amph_pred_low <- pred_p_low$amphipod; crust_pred_low <- pred_p_low$crustacean
pol_pred_low <- pred_p_low$polychaete; tanad_pred_low <- pred_p_low$tanaidacea; isop_pred_low <- pred_p_low$isopod; fish_pred_low <- pred_p_low$fish

predicted_mean <- c(amph_pred_mean, crust_pred_mean, pol_pred_mean, tanad_pred_mean, isop_pred_mean, fish_pred_mean)
predicted_up <- c(amph_pred_up, crust_pred_up, pol_pred_up, tanad_pred_up, isop_pred_up, fish_pred_up)
predicted_low <- c(amph_pred_low, crust_pred_low, pol_pred_low, tanad_pred_low, isop_pred_low, fish_pred_low)

amph_obs <- ob_p$amphipod
crust_obs <- ob_p$crustacean
pol_obs <- ob_p$polychaete
tanad_obs <- ob_p$tanaidacea
isop_obs <- ob_p$isopod
fish_obs <- ob_p$fish

observed <- c(amph_obs, crust_obs, pol_obs, tanad_obs, isop_obs, fish_obs)

plot(observed ~ predicted_mean)
abline(lm(observed ~ predicted_mean), col = "blue")
abline(lm(observed ~ predicted_low), col = "blue", lty = "dashed")
abline(lm(observed ~ predicted_up), col = "blue", lty = "dashed")
abline(0,1, col = "red")
legend(x = "topleft",          
       legend = c("Expected", "Predicted"),
       col = c("blue", "red"), 
       lwd = 2)  

ppcheck_MICUND <- data.frame(predicted_mean, predicted_low, predicted_up, observed)

## FORMATTING GGS DATAFRAME ----------------------------------------------------

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

mcmc_sum <- mcmc_df %>% group_by(Prey, Site) %>%
  summarise(median = median(value), 
            low = quantile(value, probs = c(0.05)),
            upp = quantile(value, probs = c(0.95)))


ggplot() +
  geom_violinhalf(data = mcmc_df, aes(y = value, x = Site),
                  trim = TRUE, scale = "width",
                  fill = "lightgreen", alpha = 0.3,
                  position = "dodge") +
  facet_wrap(~Prey) + 
  theme_custom() + geom_errorbar(data = mcmc_sum, 
                                 aes(ymin = mcmc_sum$low,
                                     ymax = mcmc_sum$upp, x = Site),
                                 width = 0, linewidth = 0.8) +
  geom_point(data = mcmc_sum, aes(x = Site, y = median)) +
  ylab("Posterior PoE") + xlab("")


## PLOTTING RAW FREQUENCY TABLES -----------------------------------------------
raw_freq_MICUND <- matrix(NA, nrow = dataList_MICUND$Ntotal, ncol = dataList_MICUND$N_prey)

for(i in 1:(nrow(raw_freq_MICUND))) {
  for (j in 1:(ncol(raw_freq_MICUND))) {
    raw_freq_MICUND[i,j] <- (dataList_MICUND$y[i,j]/dataList_MICUND$N_obs_site[i])*100
  }}

colnames(raw_freq_MICUND) <- c("Amphipods", "Other Crustaceans", "Polychaete",
                               "Tanaidacea", "Isopods", "Fish")

raw_freq_MICUND <- round(raw_freq_MICUND, 0)
raw_freq_MICUND <- data.frame(cbind(site_treat, as.data.frame(raw_freq_MICUND)))

names(raw_freq_MICUND) <- c("Site", "Amphipods", "Other Crustaceans", "Polychaete",
                            "Tanaidacea", "Isopods", "Fish")

MICUND_rawtable <- flextable(
  data = raw_freq_MICUND, 
  col_keys = c("Site", "Amphipods", "Other Crustaceans", "Polychaete",
               "Tanaidacea", "Isopods", "Fish")) %>%
  colformat_num(suffix = "%") %>%
  bold(part = "header", bold = TRUE) %>%
  align(align = "center", part = "all")%>%
  set_caption(caption = "Croaker",
              align_with_table = TRUE
  )

# SILVER PERCH -----------------------------------------------------------------
## FORMATING RAW DATA ----------------------------------------------------------
diet.raw <- Seine_Data_DT_1_

BAICHR <- diet.raw %>%
  filter(diet.raw$`Species code` == "BAICHR")


BAICHR <- BAICHR %>%
  group_by(`Fish ID`,Year, Treatment, Group, Length, Site) %>%
  summarize(count = n())

BAICHR <- BAICHR %>%
  filter(Site %in% c("CI", "NEPaP", "SA"))

BAICHR <- BAICHR %>%
  pivot_wider(names_from = Group, values_from = count, values_fill = 0)


BAICHR.df <- data.frame(BAICHR$Treatment, BAICHR$Site, BAICHR$Length, 
                        BAICHR$Amphipod, BAICHR$Crustacean,
                        BAICHR$Polychaete, BAICHR$Isopod, BAICHR$Tanaidacea,
                        BAICHR$Fish)

names(BAICHR.df) <- c('treat', 'site', 'TL', 'amphipod', 'crustacean',
                      'polychaetae', 'isopod', 'tanaidacea', 'fish')

BAICHR.df <- BAICHR.df %>% group_by(site, treat) %>% 
  summarise(TL = mean(TL, na.rm = TRUE), amphipod = sum(amphipod), crustacean = sum(crustacean),
            polychaete = sum(polychaetae), tanaidacea = sum(tanaidacea), isopod = sum(isopod),
            fish = sum(fish))


BAICHR.df <- BAICHR.df[-c(2,5),]

## FORMATTING DATA FOR MODEL INPUTS --------------------------------------------

y = as.matrix(BAICHR.df[,4:9]) # response variables
site = as.numeric(as.factor(BAICHR.df$site)) # categorical covariate 01
treat = as.numeric(as.factor(BAICHR.df$treat)) 

# number of observations ('visits') for each site for binomial distribution
N_obs_site1 <- BAICHR %>% 
  group_by(Site, Treatment) %>%
  summarise(no_rows = length(Site))
N_obs_site <- N_obs_site1$no_rows

N_obs_site <- N_obs_site[-c(2,5)]

# constrain prey observations to <= No of observations per site for binomial distribution
y <- ifelse(y > N_obs_site, (N_obs_site-1), y)

sitelevels = levels(as.factor(BAICHR.df$site)) # levels for categorical covariate
treatlevels = levels(as.factor(BAICHR.df$treat)) 
Ntotal = length(BAICHR.df$site) # number of sites
N_lvls_site = length(unique(site)) # number of levels for categorical covariate
N_lvls_treat = length(unique(treat))

# Specify the data in a list for sending to JAGS:
dataList_BAICHR = list(
  y = y,
  TL = BAICHR.df$TL,
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

## WRITING MODEL ---------------------------------------------------------------

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
params <- c("p", "ov_POE")

nt = 1; nc = 3; nb = 1000; ni = 10000

fit_MBB_BAICHR <- jags(data = dataList_BAICHR,
                   inits = NULL,
                   parameters = params,
                   model.file = model.spec,
                   n.thin = nt,
                   n.chains = nc,
                   n.burnin = nb,
                   n.iter = ni)

print(fit_MBB_BAICHR, intervals=c(0.2, 0.8), digits=2)

## POSTERIOR PREDICTIVE CHECKS -------------------------------------------------

ob_p <- as.data.frame(dataList_BAICHR$y/dataList_BAICHR$N_obs_site) # observed PoEs

# get predicted PoEs
pred_p <- matrix(NA, nrow = dataList_BAICHR$Ntotal, ncol = dataList_BAICHR$N_prey)
pred_p_up <- matrix(NA, nrow = dataList_BAICHR$Ntotal, ncol = dataList_BAICHR$N_prey)
pred_p_low <- matrix(NA, nrow = dataList_BAICHR$Ntotal, ncol = dataList_BAICHR$N_prey)
for(i in 1:(dataList_BAICHR$N_prey)) {
  for (j in 1:(dataList_BAICHR$Ntotal)) {
    pred_p[j,i] <- fit_MBB_BAICHR$BUGSoutput$mean$p[j,i] # means
    pred_p_up[j,i] <- fit_MBB_BAICHR$BUGSoutput$mean$p[j,i] + (2*fit_MBB_BAICHR$BUGSoutput$sd$p[j,i])
    pred_p_low[j,i] <- fit_MBB_BAICHR$BUGSoutput$mean$p[j,i] - (2*fit_MBB_BAICHR$BUGSoutput$sd$p[j,i])
  }
}

pred_p <- as.data.frame(pred_p)
pred_p_up <- as.data.frame(pred_p_up)
pred_p_low <- as.data.frame(pred_p_low)

names(pred_p) <- c("amphipod", "crustacean", "polychaete", "tanaidacea", "isopod", "fish")
names(pred_p_up) <- c("amphipod", "crustacean", "polychaete", "tanaidacea", "isopod", "fish")
names(pred_p_low) <- c("amphipod", "crustacean", "polychaete", "tanaidacea", "isopod", "fish")

#means
amph_pred_mean <- pred_p$amphipod;crust_pred_mean <- pred_p$crustacean
pol_pred_mean <- pred_p$polychaete;tanad_pred_mean <- pred_p$tanaidacea; isop_pred_mean <- pred_p$isopod; fish_pred_mean <- pred_p$fish

#ups
amph_pred_up <- pred_p_up$amphipod; crust_pred_up <- pred_p_up$crustacean
pol_pred_up <- pred_p_up$polychaete; tanad_pred_up <- pred_p_up$tanaidacea; isop_pred_up <- pred_p_up$isopod; fish_pred_up <- pred_p_up$fish

#lows
amph_pred_low <- pred_p_low$amphipod; crust_pred_low <- pred_p_low$crustacean
pol_pred_low <- pred_p_low$polychaete; tanad_pred_low <- pred_p_low$tanaidacea; isop_pred_low <- pred_p_low$isopod; fish_pred_low <- pred_p_low$fish

# concatenate into one vector
predicted_mean <- c(amph_pred_mean, crust_pred_mean, pol_pred_mean, tanad_pred_mean, isop_pred_mean, fish_pred_mean)
predicted_up <- c(amph_pred_up, crust_pred_up, pol_pred_up, tanad_pred_up, isop_pred_up, fish_pred_up)
predicted_low <- c(amph_pred_low, crust_pred_low, pol_pred_low, tanad_pred_low, isop_pred_low, fish_pred_low)

#observed
amph_obs <- ob_p$amphipod; crust_obs <- ob_p$crustacean ;pol_obs <- ob_p$polychaete
tanad_obs <- ob_p$tanaidacea; isop_obs <- ob_p$isopod; fish_obs <- ob_p$fish

# concatenate again
observed <- c(amph_obs, crust_obs, pol_obs, tanad_obs, isop_obs, fish_obs)

plot(observed ~ predicted_mean)
abline(lm(observed ~ predicted_mean), col = "blue")
abline(lm(observed ~ predicted_low), col = "blue", lty = "dashed")
abline(lm(observed ~ predicted_up), col = "blue", lty = "dashed")
abline(0,1, col = "red")
legend(x = "topleft",          
       legend = c("Expected", "Predicted"),
       col = c("blue", "red"), 
       lwd = 2)  

ppcheck_BAICHR <- data.frame(predicted_mean, predicted_low, predicted_up, observed)

## FORMATTING GGS DATAFRAME ----------------------------------------------------

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

mcmc_sum <- mcmc_df %>% group_by(Prey, Site) %>%
  summarise(median = median(value), 
            low = quantile(value, probs = c(0.05)),
            upp = quantile(value, probs = c(0.95)))


ggplot() +
  geom_violinhalf(data = mcmc_df, aes(y = value, x = Site),
                  trim = TRUE, scale = "width",
                  fill = "lightgreen", alpha = 0.3,
                  position = "dodge") +
  facet_wrap(~Prey) + 
  theme_custom() + geom_errorbar(data = mcmc_sum, 
                                 aes(ymin = mcmc_sum$low,
                                     ymax = mcmc_sum$upp, x = Site),
                                 width = 0, linewidth = 0.8) +
  geom_point(data = mcmc_sum, aes(x = Site, y = median)) +
  ylab("Posterior PoE") + xlab("")


## PLOTTING RAW FREQUENCY TABLES -----------------------------------------------
raw_freq_BAICHR <- matrix(NA, nrow = dataList_BAICHR$Ntotal,
                          ncol = dataList_BAICHR$N_prey)

for(i in 1:(nrow(raw_freq_BAICHR))) {
  for (j in 1:(ncol(raw_freq_BAICHR))) {
    raw_freq_BAICHR[i,j] <- (dataList_BAICHR$y[i,j]/dataList_BAICHR$N_obs_site[i])*100
  }}

colnames(raw_freq_BAICHR) <- c("Amphipods", "Other Crustaceans", "Polychaete",
                               "Tanaidacea", "Isopods", "Fish")

raw_freq_BAICHR <- round(raw_freq_BAICHR, 0)
raw_freq_BAICHR <- data.frame(cbind(site_treat,
                                    as.data.frame(raw_freq_BAICHR)))

names(raw_freq_BAICHR) <- c("Site", "Amphipods", "Other Crustaceans", "Polychaete",
                            "Tanaidacea", "Isopods", "Fish")

BAICHR_rawtable <- flextable(
  data = raw_freq_BAICHR, 
  col_keys = c("Site", "Amphipods", "Other Crustaceans", "Polychaete",
               "Tanaidacea", "Isopods", "Fish")) %>%
  colformat_num(suffix = "%") %>%
  bold(part = "header", bold = TRUE) %>%
  align(align = "center", part = "all")%>%
  set_caption(caption = "Silver perch",
              align_with_table = TRUE
  )


# PP CHECK FOR ALL -------------------------------------------------------------
ppcheck_LAGRHO <- ppcheck_LAGRHO %>% mutate(Species = "Pinfish")
ppcheck_MICUND <- ppcheck_MICUND %>% mutate(Species = "Croaker")
ppcheck_BAICHR <- ppcheck_BAICHR %>% mutate(Species = "Silver perch")

ppcheck.df <- rbind(ppcheck_BAICHR, ppcheck_LAGRHO, ppcheck_MICUND)

ppcheck.df$Species <- ordered(ppcheck.df$Species, levels = c("Pinfish", "Croaker", "Silver perch"))

cols <- c("Data"="black","Mean predicted"="black","95% HDI"="black","Observed"="red")
line_types <- c("Data"=NA, "Mean predicted"=2, "95% HDI"=1, "Observed"=2)

ggplot(data = ppcheck.df) +
  geom_point(aes(x = predict(lm(observed ~ predicted_mean)),
                 y = observed, color = "Data"),
             alpha = 0.5, size = 1.7) +
  facet_wrap(~Species) +
  stat_smooth(aes(x = predicted_low, y = observed,
                  color = "Mean predicted"), method = "lm", se = FALSE,
              linetype = "dashed", linewidth = 0.6)+
  stat_smooth(aes(x = predicted_up, y = observed), method = "lm", se = FALSE,
              color = "black", linetype = "dashed", linewidth = 0.6)+
  stat_smooth(aes(x = predicted_mean, y = observed,
                  color = "95% HDI"), method = "lm", se = FALSE)+
  stat_smooth(aes(x = observed, y = observed, color = "Observed"),
              method = "lm", se = FALSE,
              linetype = "dashed", linewidth = 0.6)+
  theme_custom()+
  scale_color_manual(values = cols)+
  labs(y = "Predicted", x = "Observed", color = "") +
  guides(color = guide_legend(
    override.aes = list(shape = c(NA, 16, NA, NA), 
                        linetype = c(2, NA, 1, 2))))+
  theme(legend.position="top")

# OVERALL PoE for all species---------------------------------------------------

## Pinfish-----
OV_PoE_mcmc <- ggs(as.mcmc(fit_MBB_lAGRHO))


OV_PoE_mcmc <- filter(OV_PoE_mcmc, Parameter %in% c("ov_POE[1]", "ov_POE[2]", "ov_POE[3]",
                                         "ov_POE[4]", "ov_POE[5]","ov_POE[6]", "ov_POE[7]"))

OV_PoE_mcmc <- OV_PoE_mcmc %>%
  mutate(Prey = case_when(Parameter == 'ov_POE[1]' ~ 'SAV',
                          Parameter == 'ov_POE[2]' ~ 'Amphipod',
                          Parameter == 'ov_POE[3]' ~ 'Crustacean',
                          Parameter == 'ov_POE[4]' ~ 'Polychaete',
                          Parameter == 'ov_POE[5]' ~ 'Tanaidacea',
                          Parameter == 'ov_POE[6]' ~ 'Isopod',
                          Parameter == 'ov_POE[7]' ~ 'Fish'))

OV_PoE_mcmc <- OV_PoE_mcmc %>%
  group_by(Prey) %>%
  summarise(median = median(value),
            sd = sd(value))

OV_PoE_mcmc$Prey <- ordered(OV_PoE_mcmc$Prey, levels = c("Amphipod",
                                                         "Crustacean", "SAV",
                                                         "Polychaete", "Fish",
                                                         "Isopod", "Tanaidacea"))


OV_PoE_plot_LAGRHO <- ggplot(OV_PoE_mcmc) +
  geom_errorbar(aes(x = Prey, ymin = median - (2*sd),
                    ymax = median + (2*sd)), width = 0.1) +
  geom_point(aes(x = Prey, y = median), size = 1.8, color = "lightgreen") + ylim(0,.6) +
  theme_custom() + ylab("PoE") + xlab('') + ggtitle("Pinfish")

OV_PoE_plot_LAGRHO

## Croaker-----
OV_PoE_mcmc <- ggs(as.mcmc(fit_MBB_MICUND))


OV_PoE_mcmc <- filter(OV_PoE_mcmc, Parameter %in% c("ov_POE[1]", "ov_POE[2]", "ov_POE[3]",
                                                    "ov_POE[4]", "ov_POE[5]","ov_POE[6]"))

OV_PoE_mcmc <- OV_PoE_mcmc %>%
  mutate(Prey = case_when(Parameter == 'ov_POE[1]' ~ 'Amphipod',
                          Parameter == 'ov_POE[2]' ~ 'Crustacean',
                          Parameter == 'ov_POE[3]' ~ 'Polychaete',
                          Parameter == 'ov_POE[4]' ~ 'Tanaidacea',
                          Parameter == 'ov_POE[5]' ~ 'Isopod',
                          Parameter == 'ov_POE[6]' ~ 'Fish'))

OV_PoE_mcmc <- OV_PoE_mcmc %>%
  group_by(Prey) %>%
  summarise(median = median(value),
            sd = sd(value))

OV_PoE_mcmc$Prey <- ordered(OV_PoE_mcmc$Prey, levels = c("Crustacean", "Amphipod",
                                                         "Fish","Polychaete",
                                                         "Isopod", "Tanaidacea"))


OV_PoE_plot_MICUND <- ggplot(OV_PoE_mcmc) +
  geom_errorbar(aes(x = Prey, ymin = median - (1.5*sd),
                    ymax = median + (2*sd)), width = 0.1) +
  geom_point(aes(x = Prey, y = median), size = 1.8, color = "lightgreen") +
  theme_custom() + ylab("PoE") + xlab('') + ggtitle("Croaker")

OV_PoE_plot_MICUND 

## Silver perch -----
OV_PoE_mcmc <- ggs(as.mcmc(fit_MBB_BAICHR))


OV_PoE_mcmc <- filter(OV_PoE_mcmc, Parameter %in% c("ov_POE[1]", "ov_POE[2]", "ov_POE[3]",
                                                    "ov_POE[4]", "ov_POE[5]","ov_POE[6]"))

OV_PoE_mcmc <- OV_PoE_mcmc %>%
  mutate(Prey = case_when(Parameter == 'ov_POE[1]' ~ 'Amphipod',
                          Parameter == 'ov_POE[2]' ~ 'Crustacean',
                          Parameter == 'ov_POE[3]' ~ 'Polychaete',
                          Parameter == 'ov_POE[4]' ~ 'Tanaidacea',
                          Parameter == 'ov_POE[5]' ~ 'Isopod',
                          Parameter == 'ov_POE[6]' ~ 'Fish'))

OV_PoE_mcmc <- OV_PoE_mcmc %>%
  group_by(Prey) %>%
  summarise(median = median(value),
            sd = sd(value))

OV_PoE_mcmc$Prey <- ordered(OV_PoE_mcmc$Prey, levels = c("Crustacean", "Amphipod",
                                                         "Polychaete", "Fish",
                                                         "Isopod", "Tanaidacea"))


OV_PoE_plot_BAICHR <- ggplot(OV_PoE_mcmc) +
  geom_errorbar(aes(x = Prey, ymin = median - (1.5*sd),
                    ymax = median + (2*sd)), width = 0.1) +
  geom_point(aes(x = Prey, y = median), size = 1.8, color = "lightgreen") +
  theme_custom() + ylab("PoE") + xlab('') + ggtitle("Silver perch")

OV_PoE_plot_BAICHR

# arrange multiple panels
require(ggpubr)

ggarrange(OV_PoE_plot_LAGRHO, OV_PoE_plot_MICUND, OV_PoE_plot_BAICHR, ncol = 1)

# cool ggplot theme ----
theme_custom <- function(){ 
  font <- "sans"   #assign font family up front
  
  theme_minimal() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.major = element_line(color = "lightgrey",
                                      linewidth = 0.2,
                                      linetype = 11),
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_blank(),          #strip axis ticks
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 20,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 14),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 10),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 9),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10))
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}







