gev.est <- function(y,               # data
                    mu.init, sigma.init, xi.init,  # starting values
                    v.mu, v.phi, v.xi, # prior variances
                    nu, # covariance matrix of the proposal distribution
                    niter=10000, # number of total iterations
                    burn=1000, # number burn-in iterations
                    seed=123){
  
  # compute log posterior (apart from a constant)
  p <- function(mu, phi, xi) { 
    sum(log(dgev(y, mu, exp(phi), xi))) +
      dnorm(mu, sd = sqrt(v.mu), log = TRUE) +
      dnorm(phi, sd = sqrt(v.phi), log = TRUE) +
      dnorm(xi, sd = sqrt(v.xi), log = TRUE)
  }
  
  set.seed(seed)
  
  # initialization
  mu <- phi <- xi <- numeric()
  mu[1] <- mu.init; phi[1] <- log(sigma.init); xi[1] <- xi.init 
  
  # Metropolis--Hastings for GEV distribution
  for (i in 1:niter) { 
    # generate from proposal
    thetastar <- MASS::mvrnorm(1, mu = c(mu[i], phi[i], xi[i]), Sigma = nu)
    mustar <- thetastar[1]
    phistar <- thetastar[2]
    xistar <- thetastar[3]
    
    # evaluate log posterior of proposal vs current
    lpstar <- p(mustar, phistar, xistar)
    lp <- p(mu[i], phi[i], xi[i])
    alpha <- exp(lpstar - lp)
    
    # accept or reject
    if (alpha > runif(1)) { 
      mu[i + 1] <- mustar; phi[i + 1] <- phistar; xi[i + 1] <- xistar
    } else { 
      mu[i + 1] <- mu[i]; phi[i + 1] <- phi[i]; xi[i + 1] <- xi[i]
    }
  }
  
  return(list('mu' = mu[1:niter], 'sigma' = exp(phi[1:niter]), 'xi' = xi[1:niter]))
  
}


######################################################################################
library(LaplacesDemon)
library(mvtnorm)

mix.est <-  function(y, # data
                     K, # maximum number of components
                     niter = 1000, # number of iterations
                     a_alpha = 1, b_alpha = 1, # hyperparameters of the Gamma distribution for alpha
                     p.star = 0.234, # target acceptance probability
                     v.mu = 1e4, v.phi = 1e4, v.xi = 100,  # prior variances
                     seed = 123, ...){
  
  n <- length(y)
  set.seed(seed)
  
  # constants for the adaptive update 
  zeta0 <- -qnorm(p.star/2)
  c <- (1-1/3)*sqrt(2*pi)*exp(zeta0^2/2)/(2*zeta0) + 1/(3*p.star*(1-p.star))
  n0 <- round(5/(p.star*(1-p.star)))
  
  # initialisation of the parameters of the mixture model
  pi <- rep(1/K, K)
  mu <- runif(K, quantile(y, 0.25), quantile(y, 0.75))
  phi <- log(runif(K, 0, sd(y)))
  xi <- rnorm(K, 0, 2)
  alpha <- rgamma(1, a_alpha, b_alpha)
  
  # initialisation of the parameters for the adaptive updates
  Sigma <- array(diag(3), dim=c(3,3,K))
  tau <- rep(1, K)
  tau.start <- tau
  p.acc <- numeric()
  thetaM <- matrix(nrow=K, ncol=3)
  
  # output: a list with the values of the parameters at each iteration
  out <- list(mu = matrix(nrow = niter, ncol = K),
              phi = matrix(nrow = niter, ncol = K),
              xi = matrix(nrow = niter, ncol = K),
              pi = matrix(nrow = niter, ncol = K),
              z = matrix(nrow = niter, ncol = n), 
              alpha = rep(NA, niter))
  out$mu[1,] <- mu
  out$phi[1,] <- phi
  out$xi[1,] <- xi
  out$pi[1,] <- pi
  out$alpha[1] <- alpha
  
  for(i in 2:niter){
    ## 1) ALLOCATE EACH OBSERVATION TO A COMPONENT
    # probability of z given y
    p.z <- matrix(0, K, n)
    for (j in 1:K) {
      p.z[j,] <- pi[j]*dgev(y, mu[j], exp(phi[j]), xi[j])
    }
    p.z.given.y <- apply(p.z, 1, function(x) x/colSums(p.z))
    
    # generate z based on probability of z given y
    z <- rep(0, n)
    for(u in 1:n){
      z[u] <- sample(1:K, size = 1, prob = p.z.given.y[u,], replace=TRUE)
    }
    # save z in the output
    out$z[i,] <- z
    
    ## 2) UPDATE THE WEIGHTS
    counts <- colSums(outer(z, 1:K, FUN="=="))
    v <- rep(0, K-1)
    for(j in 1:(K - 1)){
      no1 <- FALSE
      while(!no1){ # avoid to produce values of V_h equal to 1 according to the machine
        v[j] <- rbeta(1, 1 + counts[j], alpha + sum(counts[(j + 1):K]))
        if(v[j]!=1) no1 <- TRUE 
      }
    }
    pi <- Stick(v)
    # update alpha
    alpha <- rgamma(1, a_alpha + K - 1, b_alpha - sum(log(1-v)))
    # save alpha in the output
    out$alpha[i] <- alpha
    
    ## 3) UPDATE PARAMETERS OF EACH COMPONENT VIA METROPOLIS-HASTINGS
    log.post <- function(y, mu, phi, xi){
      sum(dgev(y, mu, exp(phi), xi, log = TRUE)) +
        dnorm(mu, sd = sqrt(v.mu), log = TRUE) +
        dnorm(phi, sd = sqrt(v.phi), log = TRUE) +
        dnorm(xi, sd = sqrt(v.xi), log = TRUE)
    }
    
    
    for(j in 1:K){
      y.j <- y[z==j] # cluster data based on allocations
      mu.old <- out$mu[i-1,j]; phi.old <- out$phi[i-1,j]; xi.old <- out$xi[i-1,j] 
      # generate from proposal
      thetastar <- mvrnorm(1, mu = c(mu.old, phi.old, xi.old), Sigma = tau[j]^2*Sigma[,,j])
      mustar <- thetastar[1]
      phistar <- thetastar[2]
      xistar <- thetastar[3]
      
      # check if phistar is too big (it can happen due to huge variance in small components)
      if(phistar>500 | phistar < -500){
        mu[j] <- mu.old; phi[j] <- phi.old; xi[j] <- xi.old
        p.acc[j] <- 0
      }else{
        # check if a GEV with these parameters is defined, otherwise reject
        d <- dgev(y.j, mustar, exp(phistar), xistar)
        if(any(d<=0)){
          mu[j] <- mu.old; phi[j] <- phi.old; xi[j] <- xi.old
          p.acc[j] <- 0
        } else{
          # evaluate log posterior of proposal vs current
          lpstar <- log.post(y.j, mustar, phistar, xistar)
          lp <- log.post(y.j, mu.old, phi.old, xi.old)
          a <- min(1, exp(lpstar - lp))
          p.acc[j] <- a
          
          # accept or reject
          if(a > runif(1)) { 
            mu[j] <- mustar; phi[j] <- phistar; xi[j] <- xistar
          } else{ 
            mu[j] <- mu.old; phi[j] <- phi.old; xi[j] <- xi.old
          }
        }
      }
      
      # save GEV parameters in the output
      out$mu[i,j] <- mu[j]
      out$phi[i,j] <- phi[j]
      out$xi[i,j] <- xi[j]
      out$pi[i,j] <- pi[j]
      
      # 4) update the parameters responsible for the adaptive update of the proposal covariance matrix
      if(i<=100){
        Sigma[,,j] <- diag(3) + (tau[j]^2/i)*diag(3)
      } else{
        if (i==101) {
          theta.i <- cbind(out$mu[1:i,j], out$phi[1:i,j], out$xi[1:i,j])
          Sigma[,,j] <- cov(theta.i)
          thetaM[j,] <- colMeans(theta.i)
        } else{
          theta <- c(out$mu[i,j], out$phi[i,j], out$xi[i,j])
          thetaM2 <- ((thetaM[j,]*i)+theta)/(i+1)
          Sigma[,,j] <- (i-1)/i*Sigma[,,j] + thetaM[j,]%*%t(thetaM[j,])-(i+1)/i*thetaM2%*%t(thetaM2)+1/i*theta%*%t(theta) + (tau[j]^2/i)*diag(3)
          thetaM[j,] <- thetaM2  
        }
      }
      if(i>n0) tau[j] <- exp(log(tau[j]) + c*(p.acc[j] - p.star)/max(200, i/3))
      
    }
  }
  return(out)
}