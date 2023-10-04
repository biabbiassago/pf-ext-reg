# Parameters Sampling Function

# Sampling vector of Mu
sample.mu <- function(
    current, # vector
    z, # vector
    sigmai, # current sigma , double
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
  
  mle_mod <- extRemes::fevd(z)
  MLE_MU <- mle_mod$results$par["location"]
  fisher_info_inv <- solve(-mle_mod$results$hessian)
  
  SIGMA <- exp_cov(beta0, 1/beta1, distance_mat)
  
  M <- alpha0 + alpha1*x 
  # x is a covariate vector with an observation per each location.
  
  
  target <- function(mu_value){
    (1/sqrt((2*pi)^s*det(SIGMA)))*exp(-1/2*(mui-m))
  }
  
  proposed <- mvnorm(
    s,
    M+solve(SIGMA%*%solve(fisher_info_inv + SIGMA))%*%(MLE_MU-M),
    SIGMA - SIGMA%*%solve(fisher_info_inv + SIGMA)%*%SIGMA
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
    current,  
){
  
  #gibbs sampler
}
sample.sigma <- function(
    current,  
){
  # gibbs sampler
}


# other parameters. 
sample.alpha0 <- function(
  current,
  alpha1
){
  target <- # posterior of alpha0
  
  
  
}
sample.alpha1 <- function(
  # parameters needed    
){
  current,
  
}
sample.beta0 <- function(
  # parameters needed    
){
  
}



