library(ismev)
library(mev)

dmixgev <- function(x, K, mu, sigma, xi, pi){
  out <- 0
  for (h in 1:K) {
    out <- out + pi[h] * dgev(x, mu[h], sigma[h], xi[h])
  }
  return(out)
}


pmixgev <- function(x, K, mu, sigma, xi, pi){
  out <- 0
  for (h in 1:K) {
    out <- out + pi[h] * pgev(x, mu[h], sigma[h], xi[h])
  }
  return(out)
}



qmixgev <- function(p, K, mu, sigma, xi, pi){
  out <- NULL
  # numerical procedure to invert the CDF
  for(i in 1:length(p)){
    out <- c(out, uniroot(function(x) 
      pmixgev(x, K, mu, sigma, xi, pi) - p[i],
      lower=-1e50, upper=1e50)$root)
  }
  out
}

rmixgev <- function(n, K, mu, sigma, xi, pi){
  
  # generate the allocations 
  z <- sample(1:K, n, prob = pi, replace = TRUE)
  
  # simulate from gev distributions
  out <- rep(NA, n)
  for (h in 1:K) {
    nh <- sum(z==h)
    if(nh>0) out[z==h] <- rgev(nh, mu[h], sigma[h], xi[h]) 
  }
  
  return(list("y" = out, "z" = z))
}
