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
  pivot_wider(names_from = Group, values_from = count, values_fill = 0)


LAGRHO.df <- data.frame(LAGRHO$Treatment, LAGRHO$Site, LAGRHO$Length,
                        LAGRHO$SAV, LAGRHO$Amphipod, LAGRHO$Crustacean,
                        LAGRHO$Polychaete, LAGRHO$Tanaidacea)

names(LAGRHO.df) <- c('treat', 'site', 'TL', 'SAV', 'amphipod', 'crustacean',
                      'polychaetae', 'tanaidacea')

LAGRHO.df <- LAGRHO.df %>% group_by(site) %>% 
  summarise(TL = mean(TL), SAV = sum(SAV), amphipod = sum(amphipod), crustacean = sum(crustacean),
            polychaete = sum(polychaetae), tanaidacea = sum(tanaidacea))


# FORMATTING DATA FOR MODEL INPUTS ---------------------------------------------

y = as.matrix(LAGRHO.df[,3:7]) # response variables
x = as.numeric(as.factor(LAGRHO.df$site)) # categorical covariate

# number of observations ('visits') for each site for binomial distribution
N_obs_site1 <- LAGRHO %>% 
  group_by(Site) %>%
  summarise(no_rows = length(Site))
N_obs_site <- N_obs_site1$no_rows


# constrain prey observations to <= No of observations per site for binomial distribution
y <- ifelse(y > N_obs_site, (N_obs_site-1), y)

xlevels = levels(as.factor(LAGRHO.df$site)) # levels for categorical covariate
Ntotal = length(LAGRHO.df$site) # number of sites
N_lvls = length(unique(x)) # number of levels for categorical covariate

# Specify the data in a list for sending to JAGS:
dataList = list(
  y = y,
  N_obs_site = N_obs_site,
  x = x,
  N_prey = ncol(y),
  Ntotal = Ntotal,
  N_lvls = N_lvls 
)


modelstring = "
  model {
  
    for (k in 1:N_prey) {
      for (i in 1:Ntotal) {
    
        y[i,k] ~ dbin(p[i,k], N_obs_site[i])
  
        p[i,k] ~ dbeta(omega[x[i]]*(kappa[k]-2)+1, (1-omega[x[i]])*(kappa[k]-2)+1)

      }
    }


      for (j in 1:N_lvls) {
    
        omega[j] <- ilogit(a0 + a[j])
      
        a[j] ~ dnorm(0, 1/aSigma^2)
      
      }

    
    a0 ~ dnorm(0.0 , 1/2^2) 
    aSigma ~ dgamma(0.01, 0.01)  # mode=2, sd=4

    for (k in 1:N_prey) {
    kappaMinusTwo[k] ~ dgamma(0.01 , 0.01)  # varying K
    kappa[k] <- kappaMinusTwo[k] + 2
    }    

    # Convert a0,a[] to sum-to-zero b0,b[] :
    

      for (j in 1:N_lvls) { 

        m[j] <- a0 + a[j] 
      
      } # cell means 

    
    b0 <- mean(m[1:N_lvls])
    
    for (j in 1:N_lvls) {b[j] <- m[j] - b0}
    
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






