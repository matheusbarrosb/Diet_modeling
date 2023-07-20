install.packages("rstanarm")
require(rstanarm)
require(dplyr)
require(tidyverse)
require(rjags)
require(R2jags)
require(rstan)

rstan_options(auto_write = TRUE) # avoids recompilation of unchanged model files

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

#Dim, all have the same length
N <- length(amph)
N_sites <- length(unique(site))

data = data.frame(site = site, pol=pol, amph = amph,SAV = SAV,crust = crust,
            isop = isop, tanad = tanad, treat = treat,
            TL=TL,N=N, N_sites = N_sites)


f2 <- stan_mvmer(
  formula = list(ybern ~ year + (1 | id), albumin ~ sex + year + (year | id)),
  data = pbcLong,
  family = list(binomial, gaussian),
  chains = 1, cores = 1, seed = 12345, iter = 1000)

f2 <- stan_mvmer(
  formula = list(amph ~ treat + TL + (1 | site), pol ~ treat + TL + (1 | site)),
  data = data, family = list(binomial, binomial),
  chains = 1, cores = 1, seed = 12345, iter = 1000)


#-------------------------------------------------------------------------------
library("tidyverse")
library("rstan")


require(dplyr)
require(tidyverse)
require(rjags)
require(R2jags)
require(rstan)

rstan_options(auto_write = TRUE) # avoids recompilation of unchanged model files

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

#Dim, all have the same length
N <- length(amph)
N_sites <- length(unique(site))


(data("UKHouseOfCommons", package = "pscl"))
#> [1] "UKHouseOfCommons"
glimpse(UKHouseOfCommons)


LAGRHO_diet_data <- within(list(), {
  y <- as.matrix(data.frame(amph, pol, SAV))
  X <- model.matrix(~ 0 + TL + site + treat, data = LAGRHO) %>% scale()
  N <- nrow(y)
  K <- ncol(y)
  P <- ncol(X)
  alpha_loc <- rep(0, K)
  alpha_scale <- rep(10, K)
  beta_loc <- matrix(0, K, P)
  beta_scale <- matrix(2.5, K, P)
  Sigma_corr_shape <- 2
  Sigma_scale_scale <- 5
})


write("data {
  // multivariate outcome
  int N;
  int K;
  vector[K] y[N];
  // covariates
  int P;
  vector[P] X[N];
  // prior
  vector[K] alpha_loc;
  vector[K] alpha_scale;
  vector[P] beta_loc[K];
  vector[P] beta_scale[K];
  real Sigma_corr_shape;
  real Sigma_scale_scale;
}
parameters {
  // regression intercept
  vector[K] alpha;
  // regression coefficients
  vector[P] beta[K];
  // Cholesky factor of the correlation matrix
  cholesky_factor_corr[K] Sigma_corr_L;
  vector[K] Sigma_scale;
  // student-T degrees of freedom
  real nu;
}
transformed parameters {
  vector[K] mu[N]; // before it was vector[K]
  matrix[K, K] Sigma;
  // covariance matrix
  Sigma = crossprod(diag_pre_multiply(Sigma_scale, Sigma_corr_L));
  for (i in 1:N) {
    for (k in 1:K) {
      mu[i, k] = alpha[k] + dot_product(X[i], beta[k]);
    }
  }
}
model {
  for (k in 1:K) {
    alpha[k] ~ normal(alpha_loc[k], alpha_scale[k]);
    beta[k] ~ normal(beta_loc[k], beta_scale[k]);
  }
  nu ~ gamma(2, 0.1);
  Sigma_scale ~ cauchy(0., Sigma_scale_scale);
  Sigma_corr_L ~ lkj_corr_cholesky(Sigma_corr_shape);
  y ~ multinomial_logit(mu); // before: y ~ multi_student_t(nu, mu, Sigma);
}",
"modelString.stan")


model <- "modelString.stan"

stanc(model)

uk92_mod <- stan_model("modelString.stan")

uk92_fit <- sampling(uk92_mod, data = LAGRHO_diet_data, chains = 1)

fit <- stan(file = "modelString.stan", data = LAGRHO_diet_data, warmup = 200,
            iter = 600, chains = 3, cores = 2, thin = 1)


  
# SIMULATING DATA FOR STACKOVERFLOW QUESTION
 
# outcomes 
prey01 <- rbinom(100, 1,.5)
prey02 <- rbinom(100, 1,.35)  
prey03 <- rbinom(100, 1,.77)  

# explanatory variables, site and treat are code as dummy variables
TL <- round(runif(100, 20, 100),0)
site <- rbinom(100, 5,.5)
treat <- rbinom(100, 1,.5) 

#-------------------------------------------------------------------------------
require(dplyr)
require(tidyverse)
require(rjags)
require(R2jags)
require(rstan)

rstan_options(auto_write = TRUE) # avoids recompilation of unchanged model files

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

#Dim, all have the same length
N <- length(amph)
N_sites <- length(unique(site))

model.string <- "

model{

## PRIORS ----------------------------------------------------------------------

# random effects

mu_rd ~ dnorm(0, 0.001)

for (j in 1:N_sites) {

alpha[j] ~ dnorm(mu_rd[], 0.001)

}

# fixed effects
beta_site_amph ~ dnorm(0, 0.001)
beta_site_pol ~ dnorm(0, 0.001)
beta_site_SAV ~ dnorm(0, 0.001)
beta_site_crust ~ dnorm(0, 0.001)
beta_site_isop ~ dnorm(0, 0.001)
beta_site_tanad ~ dnorm(0, 0.001)
beta_TL_amph ~ dnorm(0, 0.001)
beta_TL_pol ~ dnorm(0, 0.001)
beta_TL_SAV ~ dnorm(0, 0.001)
beta_TL_crust ~ dnorm(0, 0.001)
beta_TL_isop ~ dnorm(0, 0.001)
beta_TL_tanad ~ dnorm(0, 0.001)
beta_treat_amph ~ dnorm(0, 0.001)
beta_treat_pol ~ dnorm(0,0.001)
beta_treat_SAV ~ dnorm(0,0.001)
beta_treat_crust ~ dnorm(0,0.001)
beta_treat_isop ~ dnorm(0,0.001)
beta_treat_tanad ~ dnorm(0,0.001)

## LIKELIHOOD ------------------------------------------------------------------

for (i in 1:N){
  for (j in 1:N_sites) {

## outcomes

amph[i] ~ dbern(p_amph[i]) 
pol[i] ~ dbern(p_pol[i])
SAV[i] ~ dbern(p_SAV[i])
crust[i] ~ dbern(p_crust[i])
isop[i] ~ dbern(p_isop[i])
tanad[i] ~ dbern(p_tanad[i])

## predictors
# I had to create different betas for each separarate model, given they are fixed
logit(p_amph[i]) <- alpha[i,j] + beta_TL_amph*TL[i] + beta_site_amph*site[i] + beta_treat_amph*treat[i]
logit(p_pol[i]) <- alpha[i,j] + beta_TL_pol*TL[i] + beta_site_pol*site[i] + beta_treat_pol*treat[i]
logit(p_SAV[i]) <- alpha[i,j] + beta_TL_SAV*TL[i] + beta_site_SAV*site[i] + beta_treat_SAV*treat[i]
logit(p_crust[i]) <- alpha[i,j] + beta_TL_crust*TL[i] + beta_site_crust*site[i] + beta_treat_crust*treat[i]
logit(p_isop[i]) <- alpha[i,j] + beta_TL_isop*TL[i] + beta_site_isop*site[i] + beta_treat_isop*treat[i]
logit(p_tanad[i]) <- alpha[i,j] + beta_TL_tanad*TL[i] + beta_site_tanad*site[i] + beta_treat_tanad*treat[i]
  
  }
} 

}
"
#Model
model.spec<-textConnection(model.string)

inits <- list (b1=0, b2=0, b3=0)


## fit model w JAGS
jags <- jags(model.file = model.spec,inits = NULL,
             parameters.to.save = c("alpha_amph","alpha_pol","beta_TL_amph","beta_TL_pol",
                                    "beta_TL_SAV", "beta_TL_crust", "beta_TL_isop",
                                    "p_pol", "p_amph","p_SAV","p_crust","p_isop","p_tanad", "b1"),
             data = list(site = site, pol=pol, amph = amph,SAV = SAV,crust = crust,
                         isop = isop, tanad = tanad, treat = treat,
                         TL=TL,N=N, N_sites = N_sites),
             n.chains=3,
             n.burnin = 1000, n.iter = 4000, n.thin = 3)

plot(jags)

plot(jags$BUGSoutput$mean$p_amph ~ site)

post.df <- data.frame(jags$BUGSoutput$mean$p_amph, jags$BUGSoutput$mean$p_pol,
                      jags$BUGSoutput$mean$p_SAV, jags$BUGSoutput$mean$p_crust,
                      jags$BUGSoutput$mean$p_tanad, jags$BUGSoutput$mean$p_isop,
                      site, treat)
names(post.df) <- c("amphipod", "polychaetae", "SAV", "crustacean",
                    "tanaidacea","isopod", "site", "treat")


theme_nice <- function() {
  theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(family = "Arial", face = "bold"),
          axis.title = element_text(family = "Arial"),
          strip.text = element_text(family = "Arial", face = "bold",
                                    size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA))
}


ggplot(data = post.df, aes(x = post.df$site, y = post.df$amphipod, fill = treat)) +
  geom_boxplot() + theme_nice() + scale_fill_brewer(palette = "Set1")
















  
  
  





