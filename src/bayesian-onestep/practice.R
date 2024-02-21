############# TESTING THE MCMC:
library(MASS)
source("src/bayesian-onestep/sampling-joints.R")
source("src/bayesian_onestep/bayesian_max_model.R")
# get data running sim-timeconst.R
z <- df$value
x <- df$x
loc <- true_coords
true_mus <- to_save$true_mus
s <- 100
#i <- rep(1,MONTHS)  


## This specifies the number of MCMC iteratiion
N.iter <- 10000


# define exponential covariance function
exp_cov <- function(sigma2,phi, distp){
  return(
    (sigma2*exp(-(distp/phi)))
  )
}


c1 <- max_model(z,x,loc,s,N.iter=1000)
chain_plotter(c1)
apply(c1$acc_rates,2,mean)


# Run again with 10000 iterations
c2 <- max_model(z,x,loc,s,N.iter=50000)
chain_plotter(c2)
saveRDS(c2,"outputs/mcmc_nov22/c2.rds")
apply(c2$acc_rates,2,mean)
saveRDS(apply(c2$acc_rates,2,mean),"outputs/mcmc_nov22/acc_rates.rds")



# Run again with 10000 iterations
c3 <- max_model(z,x,loc,s,N.iter=50000)
chain_plotter(c3)
saveRDS(c2,"outputs/mcmc_nov22/c3.rds")
apply(c3$acc_rates,2,mean)



posterior_mean_after_burnin <- list(
  "mui1" = c("true"=true_mus[1],"est"=mean_after_burnin(c2$mui[,1])),
  "mui2" = c("true"=true_mus[2],"est"=mean_after_burnin(c2$mui[,2])),
  "sigma" = c("true"=1.5,"est"=mean_after_burnin(c2$sigma)),
  "beta0" = c("true"=0.8,"est"=mean_after_burnin(c2$beta0)),
  "beta1" = c("true"=0.6,"est"=mean_after_burnin(c2$beta1)),
  "alpha0"= c("true"=0,"est"=mean_after_burnin(c2$alpha_vec[,1])),
  "alpha1"= c("true"=1.2,"est"=mean_after_burnin(c2$alpha_vec[,2]))
)

saveRDS(posterior_mean_after_burnin,"outputs/mcmc_nov22/post_mean.rds")

mean_after_burnin(c2$mui[,1])
# -0.0504331
mean_after_burnin(c2$mui[,2]) # -0.0504331



c_fixed_sigma <- max_model(z,x,loc,s,N.iter=2000)
chain_plotter(c_fixed_sigma)


c_fixed_mu <- max_model(z,x,loc,s,N.iter=2000)
chain_plotter(c_fixed_mu)


# Error sometimes?, beta0 get too large?


par(mfrow=c(1,1))
plot(c1$sigma,type="l")
abline(h=0.2,col="red")
title("sigma alone")



mean_after_burnin(c1$mui[,1])
mean_after_burnin(c1$mui[,2])
mean_after_burnin(c1$beta0)
mean_after_burnin(c1$beta1)

mean_after_burnin(c1$alpha_vec[,1])
mean_after_burnin(c1$alpha_vec[,2])





### Using Coda to Evaluate how long to run the model for
library(coda)
c1 <- max_model(z,x,loc,s,N.iter=5000)
test <-
  cbind(
    c1$mui,
    c1$sigma,
    c1$beta0,
    c1$beta1,
    c1$alpha_vec
  )
colnames(test) <- c(paste0("mui",seq(1,100)),"sigma","beta0","beta1","alpha0","alpha1")

c1_mcmc <- coda::mcmc(test)
samplesize_stats <- raftery.diag(c1_mcmc, q=0.025, r=0.005, s=0.95, converge.eps=0.001)
max(samplesize_stats$resmatrix[,"N"]) #460436











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


### First, look at models where the proposal distribution is 
### MVN(current, SIGMA) ,where 
# SIGMA = exp_cov(beta0, 1/beta1, distance_mat)
# 
# m0 <- max_model(z,x,loc,s,N.iter=1000)
# m0_mui <- apply(m0$mui,2,mean)
# print(m0$acc_mui)
# # Acceptance Rate: 0.01001001
# testing_plots(m0)
# 
# m1 <- max_model(z,x,loc,s,N.iter=10000)
# m1_mui <- apply(m1$mui,2,mean)
# print(m1$acc_mui)
# testing_plots(m1)
# #[1] 0.00390039.  ->Acceptance rate.
# 
# 
# m2 <- max_model(z,x,loc,s,N.iter=10000)
# m2_mui <- apply(m2$mui,2,mean)
# print(m2$acc_mui)
# #[1] 0.00310031.  ->Acceptance rate.
# testing_plots(m2)
# 
# ### Then look at models where the proposal distribution is 
# ### MVN(current, SIGMA) ,where 
# # SIGMA = exp_cov(0.1,2.0,distance_mat) 
# 
# m3 <- max_model(z,x,loc,s,N.iter=1000)
# m3_mui <- apply(m3$mui,2,mean)
# print(m3$acc_mui)
# testing_plots(m3)
# 
# 
# 
# m5 <- max_model(z,x,loc,s,N.iter=10000)
# m5_mui <- apply(m5$mui,2,mean)
# print(m5$acc_mui)
# testing_plots(m5)
# 
# 
# m4 <- max_model(z,x,loc,s,N.iter=100000)
# m4_mui <- apply(m4$mui,2,mean)
# print(m4$acc_mui)
# testing_plots(m4)
# 
# 
# 
# 
# 
# # Chnanged the beta sampler
# 
# 
# m6 <- max_model(z,x,loc,s,N.iter=5000)
# est_beta0 <-mean(m6$beta0)
# plot(m6$beta0,type="l")
# abline(h=0.2)
# testing_plots(m6)
# 
# 
# 
# m7 <- max_model(z,x,loc,s,N.iter=1000)
# est_beta0 <-mean(m6$beta0)
# plot(m6$beta0,type="l")
# abline(h=0.2,col="red")
# title("Beta0 samples 1000 iters")
# testing_plots(m6)
# 
# 
# ##Helpfer function
# testing_plots <- function(model){
#   par(mfrow=c(2,2))
#   plot(model$mui[,1],type="l")
#   abline(h=true_mus[1],col="red")
#   plot(model$mui[,2],type="l")
#   abline(h=true_mus[2],col="red")
#   plot(model$mui[,3],type="l")
#   abline(h=true_mus[3],col="red")
#   plot(model$mui[,4],type="l")
#   abline(h=true_mus[4],col="red")
#   ttl <- paste0("N.iter ",nrow(model$mui) + 1 ," first 4 mui_s (red line is truth)")
#   title(ttl, line = -2, outer = TRUE)
# }
# 
