write("data { 
  int<lower=1> K; // number of columns
  int<lower=1> J; // number of rows
  int<lower=0> N; // number of observations
  vector[J] x[N];
  vector[K] y[N]; 
}
parameters { 
  matrix[K, J] beta; 
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] L_sigma;
} 
model {
  vector[K] mu[N];
  for (n in 1:N) 
    mu[n] = beta * x[n];
  
  to_vector(beta) ~ normal(0, 2);
  L_Omega ~ lkj_corr_cholesky(1); 
  L_sigma ~ student_t(3, 0, 2);
  
  y ~ multi_normal_cholesky(mu, diag_pre_multiply(L_sigma, L_Omega));
}
generated quantities {
  matrix[K, K] Omega;
  matrix[K, K] Sigma;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  Sigma = quad_form_diag(Omega, L_sigma); 
}",
"multi.stan")

m <- "multi.stan"

stanc(m)

set.seed(42)
N <- 400
x <- runif(N, -1, 1)

Omega <- rbind( # correlation matrix
  c(1, 0.9),
  c(0.9, 1)
)
sigma <- c(0.6, 0.4) # residual SDs
Sigma <- diag(sigma) %*% Omega %*% diag(sigma) # covariance matrix
Sigma
errors <- mvtnorm::rmvnorm(N, c(0,0), Sigma)
plot(errors)
cor(errors) # realized correlation

y1 <- -0.5 + x * 1.1 + errors[,1]
y2 <- 0.8 + x * 0.4 + errors[,2]
plot(x, y1)
plot(x, y2)
plot(y1, y2)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
m <- stan(file = m, data = list(J = 2, K = 2, N = length(y1), 
                      x = matrix(c(rep(1, N), x), ncol = 2), 
                      y = matrix(c(y1, y2), ncol = 2)),
          iter = 300)
print(m)
plot(m)




















