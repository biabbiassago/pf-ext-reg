
temp_res <- readRDS("outputs/mcmc-exact/oct29-all-large.rds")

make_traceplots(temp_res,50000)
print_summary(temp_res,burnin=10000)



##### SUMMARIZING FNCTS #####



make_traceplots <- function(mod,nsims){
  par(mfrow=c(4,2))
  plot(1:nsims,mod$eta,"l"); title("eta; truth=4")
  abline(h=4,col="darkred",lw=2)
  plot(1:nsims,mod$lambda_star,"l"); title("Lambda* ; truth=100.")
  abline(h=150,col="darkred",lw=2)
  plot(1:nsims,mod$beta,"l"); title("beta; truth=2")
  abline(h=2,col="darkred",lw=2)
  plot(1:nsims,mod$tau2,"l"); title("tau2; truth=0.1")
  abline(h=0.1,col="darkred",lw=2)
  plot(1:nsims,mod$sigma2_S,"l"); title("sigma2; truth=3");
  abline(h=3,col="darkred",lw=2)
  plot(1:nsims,mod$rho_S,"l"); title("rho_S; truth=0.15");
  abline(h=0.15,col="darkred",lw=2)
  plot(1:nsims,mod$S_n[1,],"l"); title("S_n loc 1")
  abline(h=mod$sim_data$true_s_x[1],col="darkred",lw=2)
  plot(1:nsims,mod$S_n[40,],"l"); title("S_n loc 40")
  abline(h=mod$sim_data$true_s_x[40],col="darkred",lw=2)
}   

print_summary <- function(mod,burnin=1000){
  
  summ_mat <- matrix(NA, 6,4)
  pars <- c("eta","lambda_star","tau2","beta","sigma2_S","rho_S")
  rownames(summ_mat) <- pars
  colnames(summ_mat) <- c("mean","2.5%","97.5%","truth")
  truth <- c("eta"=4,"lambda_star"=150,"tau2"=0.1,"beta"=2,"sigma2_S"=3,"rho_S"=0.15)
  
  for(p in pars){
    summ_mat[p,] <- c(mean_and_95_ci(mod[[p]],burnin),"truth"=truth[p])
  }
  return(summ_mat)
}
mean_and_95_ci <- function(x,burnin){
  n <- length(x)
  est <- mean(x[burnin:n])
  low_ci <- quantile(x[burnin:n],0.025)
  up_ci <- quantile(x[burnin:n],0.975)
  
  return(round(c("mean"=est,low_ci,up_ci),3))
  
}



#traceplots 
par(mfrow=c(3,2))
plot(1:nsims,eta,"l"); title("eta; truth=4")
abline(h=4,col="darkred",lw=2)
plot(1:nsims,lambda_star,"l"); title("Lambda^*; truth=150")

abline(h=150,col="darkred",lw=2)
plot(1:nsims,beta,"l"); title("beta; truth=2")
abline(h=2,col="darkred",lw=2)
plot(1:nsims,tau2,"l"); title("tau2; truth=0.1")
abline(h=0.1,col="darkred",lw=2)
plot(1:nsims,S_n[1,],"l"); title("S_n at loc 1");
abline(h=true_S_x[1],col="darkred",lw=2)
plot(1:nsims,sigma2_S,"l"); title("sigma2; truth=3");
abline(h=3,col="darkred",lw=2)
par(mfrow=c(1,1))

plot(1:nsims,sigma2_S,"l"); title("sigma2; truth=3");
abline(h=3,col="darkred",lw=2)

plot(1:nsims,rho_S,"l"); title("rho_S; truth=0.15");
abline(h=3,col="darkred",lw=2)


par(mfrow=c(1,2))
hist(eta,main="Posterior of Eta");
hist(lambda_star,main="Posterior of Lambda Star")





make_priorvpost <- function(mod,burnin,nsims){
  
  pars <- c("eta","lambda_star","tau2","beta")
  truth <- c("eta"=4,"lambda_star"=100,"tau2"=0.1,"beta"=2)
  priors <- list(
    eta = rnorm(1000,0,10),
    lambda_star = 1/rgamma(1000,0.001,0.001),
    tau2 = 1/rgamma(1000,0.001,0.001),
    beta = rnorm(1000,0,1)
  )
  
  par(mfrow=c(2,2))
  for(p in pars){
    plot(density(mod[[p]][burnin:nsims]),main="");  
    lines(density(priors[[p]]),col="blue")
    abline(v=truth[p],col="darkred")
    title(paste0(p))
    legend("topright",legend=c("post","prior","truth"),col=c("black","blue","darkred"),pch=c(16,16,16))
  }
  
}       



traceplot <- function(x){plot(1:length(x),x,"l")}


resx <- print_summary(tmp_results)
saveRDS(resx,here::here("outputs/mcmc-exact/tableresuptobeta.rds"))


order_Sk <- sort(S_k)

probs <- (2/sqrt(3)*order_Sk)


#   to save results
# tmp_results <-   list(
#   lambda_star=lambda_star,
#   k=k,
#   eta=eta,
#   tau2=tau2,
#   S_n = S_n
# )
# saveRDS(tmp_results, here::here("outputs/mcmc-exact/oct18-lambdaetatau2Sn.rds"))




# plot(coords(pp)[,1],coords(pp)[,2],col="red")
# points(all_coords[,1],all_coords[,2],col="blue")
# points(obs_coords[,1],obs_coords[,2],col="green")