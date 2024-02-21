library(VGAM)
# Parameters Sampling Function

#### Sampling vector of Mu ####
sample.mui <- function(
    current, # vector
    z, # vector
    sigma, # current sigma , double
    xi, # current xi, double
    alpha0, # current alpha0, double
    alpha1, # current alpha1, double
    beta0, # current beta0, double
    beta1, # current beta1, double
    s,
    distance_mat
){
  # s is the number of locations
  
  
  #mle_mod <- extRemes::fevd(z)
  #MLE_MU <- mle_mod$results$par["location"]
  #fisher_info_inv <- solve(-mle_mod$results$hessian)
  
  SIGMA <- exp_cov(beta0, beta1, distance_mat)
  SIGMA_inv <- solve(SIGMA)
  M <- alpha0 + alpha1 * unique(x)
  
  acc_mui <- 0
  

  
  
  # x is a covariate vector with an observation per each location.
  zmat <- matrix(z, ncol = MONTHS)
  # gev_inside <- function(zmat,mu_value,sigma,xi){
  #   #helper function
  #   return(
  #     (1+xi*t(t(zmat)-mu_value)/sigma)^(-1/xi)
  #   )
  # }
  
  
  target <- function(mu_value) {
    # the part of the full conditional that contains mu : gev likelihood plus gp prior on mu
    tvec <- c()
    #  gev with different mu value for each location. 
    for (i in 1:dim(zmat)[1]) {
      tvec[i] <-
        sum(dgev(
          zmat[i, ],
          location = mu_value[i],
          scale = sigma,
          shape = xi,
          log = TRUE
        ))
    }
    
    return(
      # add prior to sum 
      sum(tvec) - 1 / 2 * t(mu_value - M) %*% SIGMA_inv %*% (mu_value - M)
    )
  }
  # propose from a mvn centered at current mu values. 
  proposed <- MASS::mvrnorm(1,
                            current,
                            exp_cov(0.01, 1, distance_mat))
  
  mh_ratio <- exp(target(proposed) - target(current))
  if (runif(1) < mh_ratio) {
    mu <- proposed # accept new candidate
    acc_mui <- acc_mui + 1
  } else {
    mu <- current # reject
  }
  return(list("mui" = mu, "acc_mui" = acc_mui))
}

#### Sample Beta0 ####
sample.beta0 <- function(
    current, # double, current value of beta0
    mui, # vector
    x, # vector of my x
    beta1, # double
    alpha_vec, # vector dim 2
    distance_mat
){
  acc_beta0 <- 0 
  # I take unique(x) bc my vector of x has the same repeated covariate for each measure of each location
  M <- alpha_vec[1] + alpha_vec[2]*unique(x) 
  R <- exp(-(distance_mat/beta1))
  R_inv <- solve(R)
  
  shape = s/2 + 3
  rate =  1/2*t(mui-M)%*%R_inv%*%(mui-M)-1/2
  
  # inverse gamma
  beta0 <- 1/rgamma(1,shape=shape,rate=rate)
  return(list("beta0"=beta0,"acc_beta0"=acc_beta0))
}


#### Sample Beta1 ####
sample.beta1 <- function(
  current,
  beta0,
  alpha0, 
  alpha1,
  mui, 
  distance_mat
){
  
  acc_beta1 <- 0

  M <- alpha0 + alpha1*unique(x) 
  n <- length(z)
  
  zmat <- matrix(z,ncol=MONTHS)
  target <- function(beta1_value){
    
    SIGMA <- exp_cov(beta0, beta1_value, distance_mat)
    log_moddet_sigma <- determinant(SIGMA, logarithm = TRUE)$modulus[1]
    # prior on the beta lognormal(-1/2,1/2)
    tgt <- - (1/2)*log_moddet_sigma - 1/2*t(mui-M)%*%solve(SIGMA)%*%(mui-M) + dlnorm(beta1_value,-1/2,1/2)
    
    return(tgt)
  }
  
  proposed <- rlnorm(
    1,
    log(current),
    0.5
  )
  
  # check that this is correct
    
  # dmvnorm(y,mean,sigma.star,log=T) + dlnorm(beta1.star.....log=T) + log(1/beta1.current)
  
  ## PR0POSAL IS LOG-NORMAL
  ## HENCE IT IS NOT SYMMETRIC. Adjust for ratio proposed/current.
  mh_ratio <- exp(target(proposed)-target(current))*proposed/current

  if(runif(1)<mh_ratio){
    beta1 <- proposed
    acc_beta1 <- acc_beta1 + 1
  } else {
    beta1 <- current# reject
  }

  return(list("beta1"=beta1,"acc_beta1"=acc_beta1))
}
#### Sample XI ####
sample.xi <- function(
    current,z,sigma,mu_value
){
  
  # currently only return 0 as we are keeping XI stable. 
  
  return(0)
}

#### Sample SIGMA ####
sample.sigma <- function(
  current,
  mui,
  z,
  xi
){
  zmat <- matrix(z,ncol=MONTHS)
  acc_sigma <- 0
  
  target <- function(sigma_value){
    tvec <- c()
    for(i in 1:dim(zmat)[1]){
      tvec[i] <- sum(dgev(zmat[i,],location=mui[i],scale=sigma_value,shape=xi,log = TRUE))
    }
  return(
      sum(tvec) + dlnorm(sigma_value,0,0.5,log=TRUE)
  )
  }
  
  
  proposed <- rlnorm(
    1,
    log(current),
    0.2
  )
  
  # mh_ratio <- exp(target(proposed)-target(current))*(proposed/current)
  # do MH RATIO ON LOG SCALE 
  log_mh_ratio <- target(proposed)-target(current) + log(proposed) - log(current)
  mh_ratio <- exp(log_mh_ratio)
  if(runif(1)<mh_ratio){
    sigma <- proposed # accept new candidate
    acc_sigma <- acc_sigma + 1
  } else {
    sigma <- current # reject
  }
  return(list("sigma"=sigma,"acc_sigma"=acc_sigma))
}

#### Sample ALPHA vector)####
sample.alpha <- function(
  current_vec,
  mui,
  beta0,
  beta1,
  x,
  distance_mat
){
  
  # current SIGMA
  acc_alpha <- 0
  SIGMA <- exp_cov(beta0, beta1, distance_mat)
  logsigmadet <- determinant(SIGMA,logarithm = TRUE)$modulus
  SIGMA_inv <- solve(SIGMA)
  
  #current V - this is the var-cov matrix of my prior on alpha. I set it to be the 2x2 id matrix for now. 
  V <- matrix(c(1,0,0,1),nrow=2,byrow=TRUE)
  logvdet <- determinant(V, logarithm = TRUE)$modulus
  V_inv <- solve(V)
  
  
  
  target <- function(alpha_vec){
    # current M
    M <- alpha_vec[1] + alpha_vec[2]*unique(x)
    tgt <- -(1/2)*logsigmadet - (1/2)*logvdet -1/2*t(mui-M)%*%SIGMA_inv%*%(mui-M)-
      (1/2)*t(alpha_vec)%*%V_inv%*%(alpha_vec)
   return(tgt) 
  }
  
  proposed <- MASS::mvrnorm(
    1,
    current_vec,
    matrix(c(0.2,0,0,0.2),nrow=2,byrow = T)
  )
  
  mh_ratio <- exp(target(proposed)-target(current_vec))
  if(runif(1)<mh_ratio){
    alpha_vec <- proposed # accept new candidate
    acc_alpha <- acc_alpha + 1
  } else {
    alpha_vec <- current_vec # reject
  }
  return(list("alpha_vec"=alpha_vec,"acc_alpha"=acc_alpha))
}



