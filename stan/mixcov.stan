data {
  int<lower=1> n;    // number of observations
  int<lower=0> p;    // number of covariates
  vector[n] y;       // observations 
  matrix[n, p] x;    // covariate matrix 
  // hyperparameters for mu
  real mu_mu0;
  real<lower=0> mu_sigma0;
  // hyperparameters for sigma
  real<lower=0> sigma_alpha0;
  real<lower=0> sigma_beta0;
  // hyperparameters for beta0
  real beta0_mu0;
  real<lower=0> beta0_sigma0;
  // hyperparameters for the other regression coefficients
  real beta_mu0;
  real<lower=0> beta_sigma0;
  
  
}

parameters {
  ordered[2] mu;
  real<lower=0> sigma[2];
  real beta0;
  vector[p] beta;
}

transformed parameters{
  vector<lower=0, upper=1>[n] pi;
  pi = inv_logit(x * beta + beta0);
}

model {
  // priors
  mu ~ normal(mu_mu0, mu_sigma0);
  sigma ~ gamma(sigma_alpha0, sigma_beta0);
  beta0 ~ normal(beta0_mu0, beta0_sigma0);
  beta ~ normal(beta_mu0, beta_sigma0);
  for (i in 1:n) {
    target += log_mix(pi[i],
                     gumbel_lpdf(y[i] | mu[2], sigma[2]),
                     gumbel_lpdf(y[i] | mu[1], sigma[1]));
  }
}

// log-likelihood to compute loo information criterion
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = log_mix(pi[i],
                     gumbel_lpdf(y[i] | mu[2], sigma[2]),
                     gumbel_lpdf(y[i] | mu[1], sigma[1]));
}
