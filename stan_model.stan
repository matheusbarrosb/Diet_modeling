// Stan linear regression model
      
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
      
      
