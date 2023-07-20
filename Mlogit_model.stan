// Multivariate logistic regression
      
      data {
      int site[232];
      int treat[232];
      real TL[232]; // fish Total Length in mm
      int N_obs; // number of individual observations, rows of matrix
      int N_prey; // number of prey, columns of matrix
      real y[232]; // matrix with prey composition
      
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
      
      y ~ logistic(alpha + beta01*TL[232] + beta02*site[232] + beta03*treat[232], sigma);

      }
      }
      
      generated quantities {
      }
      
