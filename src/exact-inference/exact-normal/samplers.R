### Function ####

sample_S_x <- function(littlen, k, sigma2_S, rho_S, tau2, beta,D, eta_coeff, cov_Sk_inv){
  
  # TODO : this could be made faster read the algorithm and ask question to Veronica 
  C <- cbind(diag(1,littlen),matrix(0,nrow=littlen,ncol=(k-littlen)))
  Sigma_y <- diag(tau2,nrow=littlen)
  
  # arguments of SN
  Sigma_star <- chol2inv(chol( t(C)%*% chol2inv(chol(Sigma_y)) %*%C + cov_Sk_inv))
  mu_star <- (Sigma_star %*% t(C) %*% chol2inv(chol(Sigma_y))) %*%(y - D%*%eta_coeff )
  
  block_diag_G <- rbind(
    cbind(diag(1,nrow=littlen),matrix(0,nrow=littlen,ncol=(k-littlen))),
    cbind(matrix(0,nrow=(k-littlen),ncol=littlen),-diag(1,(k-littlen)))
  )
  G <- (beta/sqrt(sigma2_S))*block_diag_G
  bigGamma <- diag(1,k) + G%*%Sigma_star%*%t(G)
  bigGamma_inv <- chol2inv(chol(bigGamma))
  # A is the lower diagonal matrix obtained from Cholesky decomposition of bigGamma
  A <- t(chol(bigGamma))
  littleGamma <- G%*%mu_star
  bigDeltaT <- G%*%Sigma_star
  
  #u_star <- rtmvnorm(1, mean =rep(0,k),lower=(-littleGamma),algorithm="gibbs")
  u_star <- rtmvnorm(1, mean =rep(0,k),lower=(-littleGamma),D=(A),algorithm="gibbs")
  u <- A%*%t(u_star)
  mean_z_star <- (t(bigDeltaT)%*%bigGamma_inv%*%u)
  var_z_star <- (Sigma_star-t(bigDeltaT)%*%bigGamma_inv%*%bigDeltaT)
  z_star <- MASS::mvrnorm(1,mu= mean_z_star, var_z_star)
  
  S_k_sample <- z_star + mu_star
  return(as.vector(S_k_sample))
}

sample_reg_params <- function(S_n,tau2,D, eta_prior_mean=0, eta_prior_var=1000){
  n <- length(S_n)
  p <- dim(D)[2]-1
  # Set eta prior N(0,100)
  eta_prior_mean <- rep(eta_prior_mean,p+1)  
  
  
  sigma_star_eta <- chol2inv(chol((1/tau2)*(t(D)%*%D) + (1/eta_prior_var)))
  mu_star_eta <- sigma_star_eta* ((1/tau2)* t(D)%*%(y-S_n) + (1/eta_prior_var)*(eta_prior_mean))
  
  return(rnorm(1,mu_star_eta,sqrt(sigma_star_eta)))
}

sample_tau2 <- function(y,S_n,D,eta,tau_prior_a=0.01, tau_prior_b=0.01){
  n <- length(y)

  
  tau_a <- n/2 + tau_prior_a
  tau_b <- (1/2)*sum(
    (y-S_n-as.vector(D%*%eta))^2
  ) + tau_prior_b
  
  return(1/rgamma(1,tau_a,tau_b))
}

sample_beta <- function(beta_cur, sigma2_S_cur, S_k, littlen, k, beta_prior_mean, beta_prior_var){
  beta_acc <- 0
  
  beta_prop <- rnorm(1,beta_cur,1.0)
  
  log_mh_ratio <-
    target_beta_log(
      beta_prop,
      sigma = sqrt(sigma2_S_cur),
      S_k = S_k,
      beta_prior_mean = beta_prior_mean,
      beta_prior_var = beta_prior_var,
      littlen = littlen,
      k = k
    ) -
    target_beta_log(
      beta_cur,
      sigma = sqrt(sigma2_S_cur),
      S_k = S_k,
      beta_prior_mean = beta_prior_mean,
      beta_prior_var = beta_prior_var,
      littlen = littlen,
      k = k
    )
  
  # accept or reject step
  if (runif(1) < exp(log_mh_ratio)) {
    beta_acc <- 1
    beta_cur <- beta_prop
  }
  return(list(beta=beta_cur,beta_acc=beta_acc))
}

target_beta_log <-
  function(beta,
           sigma,
           S_k,
           beta_prior_mean,
           beta_prior_var,
           littlen,
           k) {
    # vector of 1 and -1 that differentiates 1:littlen, (littlen+1):k
    I_k <- c(rep(1,littlen),rep(-1,k-littlen))
    return(
      sum(pnorm((beta/sigma)*I_k*S_k,log=T)) +
        dnorm(beta, beta_prior_mean, sqrt(beta_prior_var),log = T)
    )
  }

sample_sigma2 <- function(sigma2_cur,S_k, beta_cur, R_inv,sigma2_prior_a=0.05, sigma2_prior_b=2.5){
  ## Let the prior for Sigma2 as IG (a,b)
  sigma2_acc <- 0
  sigma2_prop <- stats::rlnorm(1,log(sigma2_cur),0.2)
  log_mh_ratio <- sigma2_mh_logratio(sigma2_cur,sigma2_prop, S_k, sigma2_prior_a, sigma2_prior_b,beta_cur,R_inv)
  
  comp <- runif(1) < exp(log_mh_ratio)
  if (comp) {
    sigma2_acc <- 1
    sigma2_cur <- sigma2_prop
  }
  return(list(sigma2=sigma2_cur,sigma2_acc=sigma2_acc))
}

sigma2_mh_logratio <- function(sigma2_cur, sigma2_prop, S_k,sigma2_prior_a, sigma2_prior_b, beta_cur,R_inv){
  k <- length(S_k)
  
  # vector of 1 and -1 that differentiates 1:littlen, (littlen+1):k
  I_k <- c(rep(1,littlen),rep(-1,k-littlen))
  
  term1 <- sum(
    pnorm(
      (beta_cur/sqrt(sigma2_prop))*I_k*S_k,
      log=T)
  ) - sum(
    pnorm(
      (beta_cur/sqrt(sigma2_cur))*I_k*S_k,
      log=T)
  )
  
  term2 <- (-length(S_k)/2 - sigma2_prior_a )*(log(sigma2_prop)-log(sigma2_cur))
  
  term3 <- -((1/2)*(t(S_k)%*%R_inv%*%S_k) + sigma2_prior_b) * (1/sigma2_prop - 1/sigma2_cur)
  
  return(term1+term2+term3)
}

sample_rho <- function(rho_cur, sigma2_cur, S_k, all_coords,rho_prior_a, rho_prior_b){
  rho_prior_a <- 20
  rho_prior_b <- 40
  rho_acc <- 0
  #print(rho_cur)
  rho_prop <- stats::rlnorm(1,log(rho_cur),0.4)
  #print(rho_prop)
  log_mh_ratio_rho <- rho_mh_logratio(rho_cur,rho_prop, sigma2_cur, S_k, rho_prior_a, rho_prior_b, all_coords)
  rinf_num <- runif(1)
  
  if (runif(1) < exp(log_mh_ratio_rho)) {
    rho_acc <- 1
    rho_cur <- rho_prop
  }
  #print(rho_acc)
  return(list(rho=rho_cur,rho_acc=rho_acc,rinf_num = rinf_num, ratio=exp(log_mh_ratio_rho)))
  
}

rho_mh_logratio <- function(rho_cur,rho_prop,sigma2_cur,S_k,rho_prior_a, rho_prior_b,all_coords){
  
  #R_prop <- matern(fields::rdist(all_coords),1/rho_prop,kappa=1)
  #R_cur <- matern(fields::rdist(all_coords),1/rho_cur,kappa=1)
  
  R_prop <- fields::Matern(fields::rdist(all_coords),range =rho_prop,smoothness = 1)
  R_cur <- fields::Matern(fields::rdist(all_coords),range=rho_cur,smoothness = 1)
  
  #R_prop <- exp_cov(1,1/rho_prop,fields::rdist(all_coords))
  #R_cur <- exp_cov(1,1/rho_prop,fields::rdist(all_coords))
  
  term1 <- (-1/2)*(det(R_prop,log=T) - det(R_cur,log=T))
  term2 <- (rho_prior_a) * (log(rho_prop)-log(rho_cur))
  term3 <- - (1/(2*sigma2_cur))*(t(S_k)%*%(chol2inv(chol(R_prop))-chol2inv(chol(R_cur)))%*%S_k) - rho_prior_b*(rho_prop-rho_cur)
  
  return(term1+term2+term3)
}



## semivariogram of S at different iterations. 
## average square difference for two closeby points -- they should not be very correlated when the rho is calculate
## monitora elementi di S per vedere come si comportano. 
# esperimento: rimuovi lo step della multivariate normale -- guarda solo a sampling sulla truncated normal e vedi cosa succede.


