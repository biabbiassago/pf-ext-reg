#install.packages("tmvtnorm")
library(tmvtnorm)
library(geoR)
library(MASS)
source(here::here("src/utils.R"))
source(here::here("src/exact-inference/samplers.R"))

nsims <-1000
area_B <- 1
# define priors
lambda_prior_a <- 0.001
lambda_prior_b <- 0.001

eta_prior_mean <- 0
eta_prior_var <- 1000
tau_prior_a <- 0.001
tau_prior_b <- 0.001

sigma2_prior_a <- 0.05
sigma2_prior_b <- 2.5
rho_prior_a <- 20
rho_prior_b <- 40

beta_prior_mean <- 0
beta_prior_var <- 1



start <- Sys.time()
# simulate data 
make_sim_data <- function(beta=2){
  # simulate _all_ the points of the homogeneous poisson process
  k_true <- c(rpois(1,150))
  true_pp <- spatstat.random::runifpoint(k_true)
  
  # simulate S_k
  dist.g <- fields::rdist(coords(true_pp))
  #S_mat <- exp_cov(3,0.15,dist.g) ## change this to matern
  #S_mat <- cov.spatial(dist.g,cov.model = "exponential",cov.pars=c(3,0.15),kappa=1)
  #R <- matern(dist.g,0.15,kappa=1)
  R <- fields::Matern(dist.g,range=0.15,smoothness = 1)
  #R <- exp_cov(1,1/0.15,dist.g)
  S_mat <- 3*R
  S_k_true <- mvrnorm(1,mu=rep(0,k_true),Sigma=S_mat)
  
  
 #S_k_true <- pracma::sqrtm(S_mat)$B%*%rnorm(k_true,0,1)
  
  probs <- pnorm(beta*S_k_true/sqrt(3))
  obs_coords <- coords(true_pp)[(sapply(probs,function(x){rbinom(1,1,x)})==1),]
  idx_keep <- as.integer(rownames(obs_coords))
  S_x <- S_k_true[idx_keep]
  y_all <- 4 + S_k_true+ rnorm(k_true,0,sqrt(0.1))
  y <- y_all[idx_keep]
  S_k <- c(S_x,S_k_true[!(1:k_true %in% idx_keep)])
  all_coords <- rbind(obs_coords,coords(true_pp)[!(1:k_true %in% idx_keep),])
  
  return(list(y=y,y_all=y_all,coords_keep=obs_coords,true_s_x=S_x,true_S_k=S_k,true_pp=true_pp,all_coords=all_coords))
}


sim_data <- make_sim_data()

## initialize data
obs_coords <- sim_data$coords_keep
names(obs_coords) <- c("x","y")
littlen <- dim(obs_coords)[1]
print(littlen)
y <- sim_data$y
# design matrix, intercept only
D <- matrix(1,nrow=littlen,1)
true_S_k <- as.vector(sim_data$true_S_k)
true_S_x <- sim_data$true_s_x  
true_all_coords <- sim_data$all_coords

## initialize values
eta <- c(mean(y))
lambda_star <- c(100)
tau2 <- c(0.1)
beta <- c(0)  
beta_acc_rate <- c(1)
sigma2_S <- c(2)
sigma2_acc_rate <- c(1)
sigma2_mh_ratio <- c(1)
sigma2_mh_rinf <- c(1)
rho_S <- c(0.3)
rho_acc_rate <- c(1)

k_star <- c(rpois(1,lambda_star[1]))
pp <- spatstat.random::runifpoint(k_star[1])

all_coords <- rbind(obs_coords,coords(pp))
k <- c(dim(all_coords)[1])

#cov_S_k <-cov.spatial(fields::rdist(all_coords),cov.model = "matern",cov.pars=c(3,0.15),kappa=1)

#R <- matern(fields::rdist(all_coords),1/0.15,kappa=1)
R <- fields::Matern(fields::rdist(all_coords),range=0.15,smoothness = 1)
#R <- exp_cov(1,1/0.15,fields::rdist(all_coords))
cov_S_k <- 3*R

S_k <- pracma::sqrtm(cov_S_k)$B%*%rnorm(dim(all_coords)[1],0,1)
S_n <- matrix(NA,nrow=littlen,ncol=nsims)
S_n[,1] <- S_k[1:littlen]



for(i in 2:nsims){
  if((i %% 1000) == 0 || i == 500){
    print(i)
    
    print(paste0("sigma2_S: ", sigma2_S[i-1]))
    print(paste0("eta: ", eta[i-1]))
    print(paste0("lambda_star ",lambda_star[i-1]))
    print(paste0("rho_S ",rho_S[i-1]))
    print(paste0("beta ",beta[i-1]))
    print(paste0("tau2",tau2[i-1]))
    
    tmp <-   list(
      lambda_star=lambda_star,
      k=k,
      eta=eta,
      tau2=tau2,
      beta=beta,
      S_n = S_n,
      sigma2_S = sigma2_S,
      rho_S = rho_S,
      sim_data = sim_data,
      other_summaries = list(beta_acc_rate=beta_acc_rate,sigma2_acc_rate=sigma2_acc_rate,rho_acc_rate=rho_acc_rate,sigma2_mh_ratio=sigma2_mh_ratio,sigma2_mh_rinf=sigma2_mh_rinf)
    )
    saveRDS(tmp, here::here("outputs/mcmc-exact/tmp.rds"))
  }
  # Step 1
  lambda_star[i] <- rgamma(1, lambda_prior_a + k[i-1], lambda_prior_b + area_B)
  
  # STEP 2: Simulate Discarded Locations
  k_star[i] <- rpois(1,lambda_star[i]*area_B)
  #k_star[i] <-dim(true_all_coords)[1]
  #if(k_star[i]==0){k_star[i]=length(obs_locations)+1} 
  pp <- spatstat.random::runifpoint(k_star[i])
  
  # cov_full <- cov.spatial(
  #   fields::rdist(rbind(coords(pp),all_coords)),
  #   cov.model = "matern",
  #   cov.pars=c(sigma2_S[i-1],rho_S[i-1]),
  #   kappa=1
  # )
  
  #R_full <- matern(fields::rdist(rbind(coords(pp),all_coords)),1/rho_S[i-1],kappa=1)
  R_full <- fields::Matern(fields::rdist(rbind(coords(pp),all_coords)),range=rho_S[i-1],smoothness = 1)
  #R_full <- exp_cov(1,1/rho_S[i-1],fields::rdist(rbind(coords(pp),all_coords)))
  cov_full <- sigma2_S[i-1]*R_full
  
  cov_12 <-  cov_full[(1:k_star[i]), (k_star[i]+1):(k[i-1]+k_star[i])]
  cov_21 <- t(cov_12)
  cov_11 <- cov_full[1:k_star[i],1:k_star[i]]
  cov_22 <- cov_full[(k_star[i]+1):(k_star[i]+k[i-1]),(k_star[i]+1):(k_star[i]+k[i-1])]
  inv_cov_22 <- chol2inv(chol(cov_22))
  S_x_tilde_mean <- cov_12%*%inv_cov_22%*%(S_k)
  S_x_tilde_cov <- cov_11 - cov_12%*%inv_cov_22%*%cov_21
  S_x_tilde <- mvrnorm(1, mu = S_x_tilde_mean, Sigma=S_x_tilde_cov)
  
  # calculate discard probabilities 
  probs <- pnorm(-beta[i-1]*S_x_tilde/sqrt(sigma2_S[i-1]))
  # select "discarded" locations according to these probabilities
  disc_coords <- coords(pp)[(sapply(probs,function(x){rbinom(1,1,x)})==T),]
  all_coords <- rbind(obs_coords, disc_coords)
  
  #all_coords <- true_all_coords
  k[i] <- dim(all_coords)[1]
  
  # cov_Sk <- cov.spatial(
  #   fields::rdist(all_coords),
  #   cov.model = "matern",
  #   cov.pars=c(sigma2_S[i-1],rho_S[i-1]),kappa=1)
  
  #R <- matern(fields::rdist(all_coords),1/rho_S[i-1],kappa=1)
  R <-fields::Matern(fields::rdist(all_coords),range=rho_S[i-1],smoothness = 1)
  R_inv <- chol2inv(chol(R))
  cov_Sk_inv <-(1/sigma2_S[i-1])*R_inv 
  
  # Step 3: sample the random effects S_x
  S_k  <- sample_S_x(littlen,k[i],sigma2_S[i-1],rho_S[i-1],tau2[i-1],beta[i-1],D, eta[i-1],cov_Sk_inv)

  #S_k <- true_S_k
  S_n[,i] <- S_k[1:littlen]

  # Step 4: sample eta an tau
  eta[i] <- sample_reg_params(S_n[,i], tau2[i-1],D,eta_prior_mean, eta_prior_var)
  tau2[i] <- sample_tau2(y,S_n[,i],D,eta[i],tau_prior_a, tau_prior_b)
  #tau2[i] <- 0.1
  #eta[i] <- 4
  
  # Step 5: sample the variance parameter
  sigma2_tmp <- sample_sigma2(sigma2_S[i-1],S_k,beta_cur =beta[i-1],R_inv,sigma2_prior_a, sigma2_prior_b)
  sigma2_S[i] <- sigma2_tmp$sigma2
  sigma2_acc_rate[i] <- sigma2_tmp$sigma2_acc
  sigma2_mh_ratio <- sigma2_tmp$ratio
  sigma2_mh_rinf <- sigma2_tmp$rinf_num
  #sigma2_S[i] <- 3
  
  # Step 6 : sample the range parameter
  rho_tmp <- sample_rho(rho_S[i-1],sigma2_S[i],S_k =S_k,all_coords = all_coords,rho_prior_a, rho_prior_b)
  rho_S[i] <- rho_tmp$rho
  rho_acc_rate[i] <- rho_tmp$rho_acc
  #rho_S[i] <- 0.15
  
  # Step 7
  beta_tmp <- sample_beta(beta[i-1],sigma2_S[i], S_k, littlen, k[i], beta_prior_mean, beta_prior_var)
  beta[i] <- beta_tmp$beta
  beta_acc_rate[i] <- beta_tmp$beta_acc
  #beta[i] <-2
}


## save to a list.
tmp_results <-   list(
  lambda_star=lambda_star,
  k=k,
  eta=eta,
  tau2=tau2,
  beta=beta,
  S_n = S_n,
  sigma2_S = sigma2_S,
  rho_S = rho_S,
  sim_data = sim_data,
  other_summaries = list(beta_acc_rate=beta_acc_rate,sigma2_acc_rate=sigma2_acc_rate,rho_acc_rate=rho_acc_rate),
  priors = list(
    lambda=paste0("Gamma(",lambda_prior_a,",",lambda_prior_b,")"),
    tau2=paste0("IG(",tau_prior_a,",",tau_prior_b,")"),
    eta=paste0("n(",eta_prior_mean,",",eta_prior_var,")"),
    sigma2=paste0("IG(",sigma2_prior_a,",",sigma2_prior_b,")"),
    rho=paste0("Gamma(",rho_prior_a,",",rho_prior_b,")"),
    beta="n(0,1)")
)
saveRDS(tmp_results, here::here("outputs/mcmc-exact/nov4-all-large-priorrho20-40.rds"))
end <- Sys.time()
print(end-start)
