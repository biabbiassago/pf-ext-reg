library(MASS)
# get data running sim-timeconst.R
z <- df$value
x <- df$x
loc <- true_coords
s <- s
#i <- rep(1,MONTHS)  

  
## This specifies the number of MCMC iteratiion
N.iter <- 10


# define exponential covariance function
exp_cov <- function(sigma2,phi, distp){
  return(
    (sigma2*exp(-(distp/phi)))
  )
}


max_model <- function(z,x,loc,s,N.iter){
  
  # distance matrix from the coordinates
  distance_mat <- fields::rdist(loc)
  
  init.mui <- rep(0.5,s)
  
  # alpha are the mean parameters
  init.alpha0 <- 0
  init.alpha1 <- 0
  
  # beta are the covariance parameters
  init.beta0 <- 0.1
  init.beta1 <- 0.1
  init.sigma <- 0.2
  init.xi <- 0.5
  for(iter in 1:N.iter){
    
  
    pars.samples <- list()
    
    
    # initial sigma and initial M
    
    if(iter==1){
      mui <- init.mui
      alpha0 <- init.alpha0
      alpha1 <- init.alpha1
      beta0 <- init.beta0
      beta1 <- init.beta1
      xi <- init.xi
      sigma <- init.sigma
    }
    
    if(iter > 1){
      
      mui <- sample.mui(mui,z,sigma,xi,alpha0,alpha1,beta0,beta1,s,distance_mat)
      sigma <- sample.sigma()
      xi <- sample.xi()
      
      alpha0 <- sample.alpha0()
      alpha1 <- sample.alpha1()
      
      beta0 <- sample.beta0()
      beta1 <- sample.beta1()
      
      
    }
    #pars.samples[iter,] <- list("mui"=mui,"alpha0"=alpha0,"alpha1"=alpha1,"beta0"=beta0,"beta1"=beta1,"xi"=xi,"sigma"=sigma)
    print(iter)
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





