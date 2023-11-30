library(MASS)
source(here::here("src/bayesian-onestep/sampling-joints.R"))
gev_inside <- function(zmat,mu_value,sigma,xi){
  #helper function
  return(
    (1+xi*t(t(zmat)-mu_value)/sigma)^(-1/xi)
  )
}
max_model <- function(z,x,loc,s,N.iter){
  
  # distance matrix from the coordinates
  distance_mat <- fields::rdist(loc)
  
  init.mui <- true_mus
  
  # alpha are the mean parameters
  init.alpha0 <- 0
  init.alpha1 <- 1.2
  
  # beta are the covariance parameters
  init.beta0 <- 0.8
  init.beta1 <- 0.6
  init.sigma <- 1.5
  init.xi <- 0
  acc_rates <- matrix(
    NA,
    ncol=5,
    nrow=N.iter-1,
    dimnames = list(1:(N.iter-1),c("mui","sigma","beta0","beta1","alpha"))
  )
  
  mu.samples <- matrix(NA,ncol=s,nrow=N.iter-1)
  beta0.samples <- matrix(NA,ncol=1,nrow=N.iter-1)
  beta1.samples <- matrix(NA,ncol=1,nrow=N.iter-1)
  xi.samples <- matrix(NA,ncol=1,nrow=N.iter-1)
  sigma.samples <- matrix(NA,ncol=1,nrow=N.iter-1)
  alpha.samples <- matrix(NA,ncol=2,nrow=N.iter-1)
  
  
  for(iter in 1:N.iter){
    
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
      
      tmp <- sample.mui(mui,z,sigma,xi,alpha0,alpha1,beta0,beta1,s,distance_mat)
      mui <- tmp$mui
      #mui <- true_mus
      mu.samples[iter-1,] <- mui
      acc_rates[iter-1,"mui"] <- tmp$acc_mui
      
      sigma_info <- sample.sigma(sigma, mui,z,xi)
      #sigma <- 1.5
      sigma.samples[iter-1,] <- sigma_info$sigma
      acc_rates[iter-1,"sigma"] <- sigma_info$acc_sigma
      xi <- sample.xi()
      xi.samples[iter-1,] <- xi
      
      alpha_info <- sample.alpha(c(alpha0,alpha1),mui,beta0,beta1,x,distance_mat)
      alpha_vec <- alpha_info$alpha_vec
      alpha0 <- alpha_vec[1]
      alpha1 <- alpha_vec[2]
      #alpha0 <- 0
      #alpha1 <- 1.2
      alpha.samples[iter-1,] <- alpha_vec
      acc_rates[iter-1,"alpha"] <- alpha_info$acc_alpha
      
      beta0_info <- sample.beta0(beta0,mui,x,beta1,alpha_vec,distance_mat)
      #beta0 <- 0.8
      beta0.samples[iter-1,] <- beta0_info$beta0
      acc_rates[iter-1,"beta0"] <- beta0_info$acc_beta0
      
      beta1_info <- sample.beta1(beta1, beta0, alpha0 ,alpha1, mui, distance_mat)
      #beta1 <- 0.6
      beta1.samples[iter-1,] <- beta1_info$beta1
      acc_rates[iter-1,"beta1"] <- beta1_info$acc_beta1
    }

    if(iter %% 1000 ==0){print(iter)}
    
    # saving file every 1e4 iterations to check if model taking too long
    if(iter %% 1e4 == 0){
      saveRDS(
      list(
        "mui"=mu.samples,
        "xi"=xi,
        "sigma"=sigma.samples,
        "alpha_vec"=alpha.samples,
        "beta0"=beta0.samples,
        "beta1"=beta1.samples,
        "acc_rates"=acc_rates
        ), here::here(paste0("src/bayesian-onestep/tmp/",iter,"iter.rds"))
      )
    }
    
    
  }
  return(list(
    "mui"=mu.samples,
    "xi"=xi,
    "sigma"=sigma.samples,
    "alpha_vec"=alpha.samples,
    "beta0"=beta0.samples,
    "beta1"=beta1.samples,
    "acc_rates"=acc_rates
  ))
}


# Other functions that might be helpful
# auxiliary functions
mean_after_burnin <- function(param,start,end){
  mean(param[start:end])
}

chain_plotter <- function(model){
  par(mfrow=c(4,2))
  plot(model$mui[,1],type="l")
  abline(h=true_mus[1],col="red")
  title("mu1")
  plot(model$mui[,2],type="l")
  abline(h=true_mus[2],col="red")
  title("mu2")
  plot(model$sigma,type="l")
  abline(h=1.5,col="red")
  title("sigma")
  plot(model$xi,type="l")
  abline(h=0,col="red")
  title("xi")
  plot(model$beta0,type="l")
  abline(h=0.8,col="red")
  title("beta0")
  plot(model$beta1,type="l")
  abline(h=0.6,col="red")
  title("beta1")
  plot(model$alpha_vec[,1],type="l")
  abline(h=0,col="red")
  title("alpha0")
  plot(model$alpha_vec[,2],type="l")
  abline(h=1.2,col="red")
  title("alpha1")
}
