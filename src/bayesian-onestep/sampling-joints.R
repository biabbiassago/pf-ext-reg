# Parameters Sampling Function

# Sampling vector of Mu
sample.mu <- function(
    current, # vector
    z, # vector
    sigma, # current sigma , double
    xi, # current xi, douuble
    alpha0, # current alpha0, double
    alpha1, # current alpha1, double
    beta0, # current beta0, double
    beta1, # current beta1, double
    fisher_info_inv,
    s,
    distance_mat
){
  # s is the number of locations
  
  
  #mle_mod <- extRemes::fevd(z)
  #MLE_MU <- mle_mod$results$par["location"]
  #fisher_info_inv <- solve(-mle_mod$results$hessian)
  
  SIGMA <- exp_cov(beta0, 1/beta1, distance_mat)
  
  M <- alpha0 + alpha1*unique(x) 
  # x is a covariate vector with an observation per each location.
  
  gev_inside <- function(z,mu_value,sigma,xi){
    return(
      (1+xi*(z-rep(mu_value,each=MONTHS)/sigma))^(-1/xi)
    )
  }
  
  
  ##### CHANGE THIS TO BE ON THE LOG SCALE
  ##### THIS CAN BE THE LOGLIK OF GEV PDF AND THEN THE SUM OF THE NORMAL PART (LOGGED it)
  # SEE https://cran.r-project.org/web/packages/metropolis/vignettes/metropolis-vignette.html
  target <- function(mu_value){
    return(prod(
      prod(
            (1/sigma)*
            abs(gev_inside(z,mu_value,sigma,xi))^(xi+1)*
            exp(-gev_inside(z,mu_value,sigma,xi))
        )
      )*exp(-1/2*t(mu_value-M)%*%solve(SIGMA)%*%(mu_value-M))
    )
  }
  
  # error 
  
  proposed <- MASS::mvrnorm(
    1,
    mui,
    SIGMA
  )
  
  mh_ratio <- target(proposed)/target(current)
  if(runif(1)<mh_ratio){
    mu = proposed # accept new candidate
  } else {
    mu = current # reject
  }
  return(mu)
}
sample.xi <- function(
){
  return(0)
  #gibbs sampler
}
sample.sigma <- function(
     
){
  # gibbs sampler
  return(1.5) 
}


# other parameters. 
sample.alpha0 <- function(
){
  return(0)
}
sample.alpha1 <- function(
  # parameters needed    
){
  return(1.2)
  
}
sample.beta0 <- function(
  # parameters needed    
){
  return(0.2)
}
sample.beta0 <- function(
    # parameters needed    
){
  return(0.6)
}



