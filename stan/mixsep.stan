data {
  int<lower=1> K;    // number of components
  int<lower=1> n;    // number of observations
  vector[n] y;       // observations
  int nvec[K+1];     // cumulate number of obs. per group
  // hyperparameters for mu
  real mu_mu0;
  real<lower=0> mu_sigma0;
  // hyperparameters for sigma
  real<lower=0> sigma_alpha0;
  real<lower=0> sigma_beta0;
}

parameters {
  ordered[K] mu;
  vector<lower=0>[K] sigma;
}

model {
  // priors
  mu ~ normal(mu_mu0, mu_sigma0);
  sigma ~ gamma(sigma_alpha0, sigma_beta0);
  for (h in 1:K) {
    for (i in (nvec[h]+1):nvec[h+1]) {
      target += gumbel_lpdf(y[i] | mu[h], sigma[h]);
  }
 }
}


