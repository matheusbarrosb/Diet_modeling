

# The factor for subjects:
n_subjects <- 5 # Number of prey
n_sample <- 100 # Number of observations for each subject
(n_obs <- n_subjects*n_sample) # Total number of data points
(subjects <- gl(n=n_subjects, k=n_sample)) # Indicator for subjects

# Continuous predictor x
original_x <- runif(n_obs, 45, 70)
summary(original_x)

# We standardize it (always good to do with continuous predictors when using JAGS/BUGS)
(mean_orig_x <- mean(original_x))
(sd_orig_x <- sd(original_x))
x <- (original_x-mean_orig_x)/sd_orig_x
round(x, 2)
summary(x)

hist(x, col="lightsteelblue1", border="white", breaks=10, freq=FALSE)
lines(density(x), col="lightsteelblue3", lwd=2)


# This is the model matrix for the means parametrization of the interaction model between the subjects and the continuous covariate x:
Xmat <- model.matrix(~subjects*x-1-x)
dim(Xmat)

# Q: where do these dimensions come from?
head(Xmat)

# -- there are 560 observations (rows) and 112 regression terms / variables (columns)
dimnames(Xmat)[[1]]
dimnames(Xmat)[[2]]

# -- there are 5 terms for subjects, the coefficients of which will provide the subject-specific intercepts (i.e., the intercept random effects, or intercept effects for short)
# -- there are 56 terms for interactions between each subject and the continuous covariate x, the coefficients of which will provide the subject-specific slopes (i.e., the slope random effects, or slope effects for short)

round(Xmat[1, ], 2) 		# Print the top row for each column
Xmat[, 1]               # Print all rows for column 1 (group 1)
round(Xmat[, 5], 2)    # Print all rows for column 57 (group 1:x)


# Parameters for the distributions of the random coefficients / random effects (note that the intercepts and slopes comes from two independent Gaussian distributions):

intercept_mean <- 0.44 # mu_alpha
intercept_sd <- 0.03 # sigma_alpha

slope_mean <- 0.05 # mu_beta
slope_sd <- 0.002 # sigma_beta


# Generate the random coefficients:

intercept_effects <- rnorm(n=n_subjects, mean=intercept_mean, sd=intercept_sd)
slope_effects <- rnorm(n=n_subjects, mean=slope_mean, sd=slope_sd)

par(mfrow=c(1, 2))
hist(intercept_effects, col="lightsteelblue1", border="white", breaks=10, freq=FALSE)
lines(density(intercept_effects), col="lightsteelblue3", lwd=2)
hist(slope_effects, col="lightsteelblue1", border="white", breaks=10, freq=FALSE)
lines(density(slope_effects), col="lightsteelblue3", lwd=2)
par(mfrow=c(1, 1))

all_effects <- c(intercept_effects, slope_effects) # Put them all together
round(all_effects, 2)

# -- thus, we have two stochastic components in our model IN ADDITION TO the usual stochastic component for the individual-level responses, to which we now turn


# Generating the continuous response variable:

# -- the deterministic part
lin_pred <- Xmat %*% all_effects # Value of lin_predictor
str(lin_pred)

# -- the stochastic part
sigma_res <- 0.01
normal_error <- rnorm(n=n_obs, mean=0, sd=sigma_res) # residuals
str(normal_error)

       # linear combination with a bias
pr = 1/(1+exp(-x))         # pass through an inv-logit function
plot(pr)
y = rbinom(n_obs,1,pr) 
plot(y~x)

# -- put the two together
y <- lin_pred+normal_error
str(y)
# or, alternatively
y <- rlogis(n=n_obs, lin_pred, sigma_res)
str(y)

# We take a look at the response variable
hist(y, col="lightsteelblue1", border="white", breaks=30, freq=FALSE)
lines(density(y), col="lightsteelblue3", lwd=2)
summary(y)

library("lattice")
xyplot(y~x|subjects)


# Analysis under a random-intercepts model

# REML analysis using R
library("lme4")
lme_fit1 <- glmer(y~x+(1|subjects),family = binomial)
print(lme_fit1, cor=FALSE)

fixef(lme_fit1)
ranef(lme_fit1)
coef(lme_fit1)

# Compare with true values:
data.frame(intercept_mean=intercept_mean, slope_mean=slope_mean, intercept_sd=intercept_sd, slope_sd=slope_sd, sigma_res=sigma_res)
print(lme_fit1, cor=TRUE)


# Bayesian analysis using JAGS

# Write model
cat("model {
    
# PRIORS
mu_int~dnorm(0, 0.0001) # Mean hyperparameter for random intercepts
sigma_int~dunif(0, 10) # SD hyperparameter for random intercepts
tau_int <- 1/(sigma_int*sigma_int)

for (i in 1:n_prey) {

    alpha[i]~dnorm(mu_int, tau_int) # Random intercepts
}

beta~dnorm(0, 0.0001) # Fixed slope
sigma_res~dunif(0, 10) # Residual standard deviation
tau_res <- 1/(sigma_res*sigma_res)

# LIKELIHOOD
for (i in 1:n_obs) {

    logit(mu[i]) <- alpha[prey[i]]+beta*TL[i] # Expectation
    y[i] ~ dbern(mu[i]) # The actual (random) responses

}
}", fill=TRUE, file="lme_model1.txt")

# Bundle data
jags_data <- list(y=LAGRHO.df2$PA, prey=as.numeric(LAGRHO.df2$prey),
                  TL=as.numeric(LAGRHO.df2$treat), n_prey=max(as.numeric(LAGRHO.df2$prey)),
                  n_obs=as.numeric(length(LAGRHO.df2$treat)))
# use as.numeric across the board for the data passed to JAGS; it might work w/o it, but this is often needed for other BUGS packages

# Inits function
inits <- function() {
  list(alpha=rnorm(n_subjects, 0, 2), beta=rnorm(1, 1, 1), mu_int=rnorm(1, 0, 1), sigma_int=rlnorm(1), sigma_res=rlnorm(1))
}

# Parameters to estimate
params <- c("alpha", "beta", "mu_int", "sigma_int", "sigma_res")

# MCMC settings
ni <- 5000; nb <- 1000; nt <- 20; nc <- 3

# Start Gibbs sampling
library("R2jags")
outj <- jags(jags_data, inits=inits, parameters.to.save=params, model.file="lme_model1.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

traceplot(outj)

# Inspect results
print(outj, dig=3)

out <- outj$BUGSoutput

# Compare with MLEs and true values
data.frame(intercept_mean=intercept_mean, slope_mean=slope_mean, intercept_sd=intercept_sd, slope_sd=slope_sd, sigma_res=sigma_res)
print(out$mean, dig=4)
print(lme_fit1, cor=FALSE)





