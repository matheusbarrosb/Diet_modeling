// Bayesian PCA
      
      data {
      int<lower=0> N; // number of observations
      int<lower=0> D; // original dimension
      int<lower=0> K; // latent dimension
      matrix[N, D] X; // data matrix
      
      }
      
      parameters {
      
      matrix[N, K] Z; // latent matrix
      matrix[D, K] W; // weight matrix
      real<lower=0> tau; //noise term
      vector<lower=0>[K] alpha; // ARD prior
      
      }
      
      transformed parameters {
      
      vector<lower=0>[K] t_alpha;
      real<lower=0> t_tau;
      t_alpha = inv(sqrt(alpha));
      t_tau = inv(sqrt(tau));
      
      }
      
      model {
      
      tau ~ gamma(1,1);
      to_vector(Z) ~ normal(0,1);
      alpha ~ gamma(1e-3,1e-3);
      for(k in 1:K) W[,k] ~ normal(0, t_alpha[k]);
      to_vector(X) ~ normal(to_vector(Z*W'), tau);
      
      }
