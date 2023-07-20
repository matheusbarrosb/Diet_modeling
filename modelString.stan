data {
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
  y ~ multinomial_logit(mu); // before: y ~ student_t(nu, mu, Sigma);
}
