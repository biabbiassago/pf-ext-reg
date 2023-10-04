options(scipen=999)
library(dplyr)


N<-400
make_sim_results <- function(rslt){
  #avg_over_sims <- lapply(rslt, function(x){apply(x,2,mean,na.rm=TRUE)})
  
  mean_over_sims <- lapply(rslt, function(x){apply(x,2,mean,na.rm=TRUE)})
  return(mean_over_sims)
}


get_no_sampling_results <- function(){
  # true_data <- read.csv(paste0("data/sim-gev-2-PHIX",PHIX,".csv"))
  # 
  # mean_by_station <- true_data %>%
  #   group_by(station) %>%
  #   summarize(
  #     mean_value = mean(value),
  #     x_coords = first(x_coords),
  #     y_coords = first(y_coords),
  #     x = first(x)
  #   )
  # 
  #   SAMPLE_SIZE <- N; b=0
    # no_sampling_results <- get_regression_estimates(b,mean_by_station,true_data)
    no_sampling_results <- readRDS(here::here(paste0("outputs/25no_sampling_results-PHIX",PHIX,".rds")))
    # no_sampling_results <- list(
    #   "beta0"=no_sampling_results$beta0_est,
    #   "beta1"=no_sampling_results$beta1_est,
    #   "sigma"=no_sampling_results$sigma_est,
    #   "phi"=no_sampling_results$phi_est,
    #   "variogsigma"=no_sampling_results$variog_cov_pars[1],
    #   "variogphi"=no_sampling_results$variog_cov_pars[2]
    # )
return(no_sampling_results)
}
# 
# PHIX<-0.5
# x <- get_no_sampling_results()
# saveRDS(x,"outputs/no_sampling_results-PHIX0.5.rds")

compute_bias <- function(rslt_summ, no_sampling_results){
  bias <- rslt_summ
  for(i in names(rslt_summ)[1:4]){
    bias[[i]]<- rslt_summ[[i]]- no_sampling_results[i,"no_samp"]
  }
  return(bias[1:4])
}




print_results <- function(){
  row_names_set <- c("beta0","beta1","sigma2","phi","krig-likfit-sigma2","krig-likfit-phi")
  print("Results for Simulation: \n")
  print(paste("PHIX=",PHIX," and SAMPLE SIZE=",SAMPLE_SIZE))
  
  
  rslt <- readRDS(
      here::here(paste0(
        "outputs/reg-est-loc-2023-07-25-sim-gev-4sampsize",
        SAMPLE_SIZE,
        "phi",
        PHIX,
        ".rds"
      )
    )
  )
  
  rslt_summ <- make_sim_results(rslt)
  rslt_summ_tbl <- data.frame(t(sapply(rslt_summ,c)))
  rownames(rslt_summ_tbl) <- row_names_set
  print("Summary of results (mean):")
  
  print(rslt_summ_tbl)

  no_samp <- get_no_sampling_results()
  
  print("Estimates without Sampling:")
  no_samp_tbl <- data.frame(sapply(no_samp[1:4],c))
  rownames(no_samp_tbl) <- row_names_set[1:4]
  colnames(no_samp_tbl) <- "no_samp"
  print(no_samp_tbl)
  
  
  bias <- compute_bias(rslt_summ, no_samp_tbl)
  print("Bias on the estimates")
  bias_tbl <- data.frame(t(sapply(bias,c)))
  rownames(bias_tbl) <- row_names_set[1:4]
  print(bias_tbl)
}



