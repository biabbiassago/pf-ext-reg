library(MASS)

z <- # matrix with value per station 
x <- # vector of covariates
loc <- # locations
s <- # number of stations
i <- rep(1,50) # 

  
## This specifies the number of MCMC iteration
N.iter <- 1000


# define exponential covariance function
exp_cov <- function(sigma2,phi, distp){
  return(
    (sigma2*exp(-(distp/phi)))
  )
}


max_model <- function(z,x,loc,s,i){
  distance_mat <- fields::rdist(loc)
  
  for(i in 1:N.iter){
    
  
    pars.samples <- matrix(NA,N.iter,6)
    
    init.mui <- rep(0.5,s)
    
    # alpha are the mean parameters
    init.alpha0 <- 0
    init.alpha1 <- 0
    
    # beta are the covariance parameters
    init.beta0 <- 0.1
    init.beta1 <- 0.1
    init.sigma <- 0.2
    init.xi <- 0.5
    
    # initial sigma and initial M
    
    if(i==1){
      mui <- init.mui
      alpha0 <- init.alpha0
      alpha1 <- init.alpha1
      beta0 <- init.beta0
      beta1 <- init.beta1
      xi <- init.xi
      sigma <- init.sigma
    }
    if(i > 1){
      
      mui <- sample.mui(mui,sigma,xi, z, M, SIGMA, s)
      sigma <- sample.sigma(..)
      xi <- sample.xi(..)
      
      
      alpha0 <- sample.alpha0(alpha0, alpha1)
      alpha1 <- sample.alpha1(alpha1, alpha0)
      
      beta0 <- sample.beta0(..)
      beta1 <- sample.beta1(..)
      
      
    }
    pars.samples[i,] <- c(mui,alpha0,alpha1,beta0,beta1,xi,sigma)
    print(i)
  }
  return(list(
    "phi"=phi,
    "xi"=xi,
    "alpha0"=alpha0,
    "alpha1"=alpha1,
    "beta0"=beta0,
    "beta1"=beta1,
    "xi"=xi
  ))
  
}





