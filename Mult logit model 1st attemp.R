library(rstan)
library(dplyr)
library(tidyverse)

rstan_options(auto_write = TRUE) # avoids recompilation of unchanged model files

diet.raw <- Seine_Data_DT

LAGRHO <- diet.raw %>%
  filter(diet.raw$`Species code` == "LAGRHO")


LAGRHO <- LAGRHO %>%
  group_by(`Fish ID`, Year, Treatment, Group, Length, Site) %>%
  summarize(count = n())

LAGRHO <- LAGRHO %>%
  pivot_wider(names_from = Group, values_from = count, values_fill = 0)

# coding factor variables as integers
site <- as.factor(LAGRHO$Site)
site <- as.numeric(site)

treat <- as.factor(LAGRHO$Treatment)
treat <- as.numeric(treat)

TL <- LAGRHO$Length

length(site)
length(treat)
length(TL)

# setting up outcome matrix
mat <- data.frame(LAGRHO$SAV, LAGRHO$Bivalve, LAGRHO$Amphipod,
                  LAGRHO$Crustacean, LAGRHO$`Green algae`, LAGRHO$Polychaete,
                  LAGRHO$Isopod, LAGRHO$Tanaidacea, LAGRHO$Cyanobacteria)

# I'm just selecting those preys for now because the data neads a bit of a clean up 

mat <- as.matrix(mat)
colnames(mat) <- NULL # remove headers, is it necessary though?

y <- mat

N_obs <- nrow(y) # number of observations
N_prey <- ncol(y) # number of prey


stan_data <- list(N_prey = N_prey, N_obs = N_obs, y = y,
                  TL = TL, site = site, treat = treat) # data for stan

# DECLARE MODEL

write("// Multivariate logistic regression
      
      data {
      
      int site;
      int treat;
      real TL; // fish Total Length in mm
      int N_obs; // number of individual observations, rows of matrix
      int N_prey; // number of prey, columns of matrix
      matrix[N_obs, N_prey] y; // matrix with prey composition
      
      }
      
      parameters {
      
      real alpha; // intercept
      real beta01; // slope for fish size
      real beta02; // slope for site
      real beta03; // slope for treat
      real <lower = 0> sigma;
      
      }
      
      model {
      
      // Priors
      alpha ~ normal(0,0.001);
      beta01 ~ normal(0,0.001);
      beta02 ~ normal(0,0.001);
      beta03 ~ normal(0,0.001);
      
      // Likelihood
      
      for (i in 1:N_obs) {
        for (j in 1:N_prey) {
      
      y[i,j] ~ binomial_logit(alpha + beta01*TL + beta02*site + beta03*treat, sigma);
     
       }
      }
      }",
      "Mlogit_model.stan")


stanc("Mlogit_model.stan") # produces error message

stan_model <- "Mlogit_model.stan"






















