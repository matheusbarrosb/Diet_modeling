# BASIC LINEAR REGRESSION EXAMPLE STAN

# simulate data

x <- seq(1,100,by=1)
a <- rep(0.2,100)
b <- rep(0.54,100) 
y <- (a + rnorm(100,0,0.05)) + (b + rnorm(100,0,0.1))*x 
N <- length(x)

plot(y~x)

stan_data <- list(N = N, y = y, x = x)

# write the model

write("// Stan linear regression model
      
      data {
      int < lower = 1 > N; // Sample size
      vector[N] x; // Predictor
      vector[N] y; // Response
      }
      
      parameters {
      real alpha; // Intercept
      real beta; // Slope
      real < lower = 0 > sigma; // SD
      }
      
      model {
      y ~ normal(alpha + x * beta, sigma);
      }
      
      ",
      "stan_model.stan")


stanc("stan_model.stan")

stan_model <- "stan_model.stan"


fit <- stan(file = stan_model, data = stan_data, warmup = 500,
            iter = 1000, chains = 4, cores = 2, thin = 1)

fit

post <- extract(fit)
str(post)


# probability that beta is > a certain value:
sum(post$beta>0.5)/length(post$beta)

stan_dens(fit)






