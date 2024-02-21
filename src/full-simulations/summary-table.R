library(tidyverse)


get_in_ci <- function(eg){
  rs <- matrix(NA,nrow=7)
  rs[1] <- sum((eg[["mu","truth"]] > eg[["mu","95ci_low"]]) & (eg[["mu","truth"]] < eg[["mu","95ci_hi"]]))/length(eg[["mu","est"]])
  
  rs[2] <-  (eg[["alpha0","truth"]] > eg[["alpha0","95ci_low"]]) & (eg[["alpha0","truth"]] < eg[["alpha1","95ci_hi"]])
  
  rs[3] <-  (eg[["alpha1","truth"]] > eg[["alpha1","95ci_low"]]) & (eg[["alpha1","truth"]] < eg[["alpha1","95ci_hi"]])
  
  
  rs[4] <-  (eg[["beta0","truth"]] > eg[["beta0","95ci_low"]]) & (eg[["beta0","truth"]] < eg[["beta1","95ci_hi"]])
  rs[5] <-  (eg[["beta1","truth"]] > eg[["beta1","95ci_low"]]) & (eg[["beta1","truth"]] < eg[["beta1","95ci_hi"]])
  
  rs[6] <- (eg[["gamma0","truth"]] > eg[["gamma0","95ci_low"]]) & (eg[["gamma0","truth"]] < eg[["gamma0","95ci_hi"]])
  rs[7] <- (eg[["gamma1","truth"]] > eg[["gamma1","95ci_low"]]) & (eg[["gamma1","truth"]] < eg[["gamma1","95ci_hi"]])
  
  rs[8] <-   sum((eg[["sigma2","truth"]] > eg[["sigma2","95ci_low"]]) & (eg[["sigma2","truth"]] < eg[["sigma2","95ci_hi"]]))/length(eg[["sigma2","est"]])
  rs[9] <-  (eg[["xi","truth"]] > eg[["xi","95ci_low"]]) & (eg[["xi","truth"]] < eg[["xi","95ci_hi"]])
  return(rs)
}

get_se <- function(eg){
  rs <- matrix(NA,nrow=7)
  rs[1] <- mean((eg[["mu","est"]] - eg[["mu","truth"]])^2)
  
  rs[2] <-  (eg[["alpha0","est"]] - eg[["alpha0","truth"]])^2
  
  rs[3] <- (eg[["alpha1","est"]] - eg[["alpha1","truth"]])^2
  
  
  rs[4] <- (eg[["beta0","est"]] - eg[["beta0","truth"]])^2
  rs[5] <- (eg[["beta1","est"]] - eg[["beta1","truth"]])^2
  
  rs[6] <- (eg[["gamma0","est"]] - eg[["gamma0","truth"]])^2
  rs[7] <- (eg[["gamma1","est"]] - eg[["gamma1","truth"]])^2
  
  
  rs[8] <-  mean((eg[["sigma2","est"]] - eg[["sigma2","truth"]])^2)
  rs[9] <-  (eg[["xi","est"]] - eg[["xi","truth"]])^2
  return(rs)
}

## USE THIS WHEN YOU HAVE ALL FILES WITH SIMS
make_summary_table <- function(S,MONTHS,inst){
path_ <- paste0("outputs/full-simulations/",inst,"/")
sim_setup <- paste0("npts",S,"stations",MONTHS,"months")
df_list <- list.files(path = here::here(path_),pattern = sim_setup,full.names=TRUE) %>%
  map(readRDS)

nms <- c("location","alpha0","alpha1","beta0","beta1","gamma0","gamma1","sigma2","xi")

cov_list <- lapply(df_list, get_in_ci)
y <- do.call(cbind,cov_list)
cvg<- apply(y,1,mean)

se_list <- lapply(df_list, get_se)
z <- do.call(cbind,se_list)
rmse <- apply(z,1,function(x) sqrt(mean(x)))

summ <- cbind(cvg,rmse)
rownames(summ) <- nms
return(summ)
}

make_summary_table(50,90,"non-pref-two-step")


