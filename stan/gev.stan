// log-density function of gev-distribution
functions {
   real gen_extreme_value_lpdf(real y, real mu, real sigma, real xi) { 
     real x = (y - mu) / sigma;
     if (xi == 0) {
       return  -log(sigma) - x - exp(-x);
     } else {
       real t = 1 + xi * x;
       real inv_xi = 1 / xi;
       if (t <= 0){
         return log(0); 
         } else{ return -log(sigma) - (1 + inv_xi) * log(t) - pow(t, -inv_xi); 
       }
      }
   }
  
}

data {
  int<lower=1> n;  // number of observations
  vector[n] y;  // observations
  // hyperparameters for mu
  real mu_mu0;
  real mu_sigma0; 
  // hyperparameters for sigma
  real lsigma_mu0;
  real lsigma_sigma0; 
  // hyperparameters for xi
  real xi_mu0;
  real xi_sigma0; 
}

parameters {
  real mu;  // location
  real lsigma; // log-scale
  real xi;  // shape
}

transformed parameters {
  real sigma; // scale
  sigma = exp(lsigma);
}

model {
  for (i in 1:n) {
      target += gen_extreme_value_lpdf(y[i] | mu, sigma, xi);
  }

  // priors including constants
  target += normal_lpdf(mu | mu_mu0, mu_sigma0);
  target += normal_lpdf(lsigma | lsigma_mu0, lsigma_sigma0);
  target += normal_lpdf(xi | xi_mu0, xi_sigma0);
}


generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = gen_extreme_value_lpdf(y[i] | mu, sigma, xi);
}

