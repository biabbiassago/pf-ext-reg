today<-"0512"


source("src/utils.R")
source("src/2.sample.R")
source("src/3.evr_bystation.R")

nsims <- 100
proposed_b <- c(5,3,0)

loc_rmse_fullsim <- function(b, mean_by_station, true_data){
  # given a data-set and a b "strength of preferential sampling - run"
  xb <- prefsamp_station(b,mean_by_station, true_data)
  evr_pred_df <- main_est_loc(xb$sample_df,true_data)
  return(evr_pred_df$rmse)
}

### ------------------- ####
#### RUN SIM #####
true_data <- read.csv("data/sim1.csv")
mean_by_station <- true_data %>%
  group_by(station) %>%
  summarize(
    mean_value = mean(value),
    x_coords = first(x_coords),
    y_coords = first(y_coords),
    x = first(x)
)


rmse_results <- matrix(NA,nrow=nsims,ncol=3)
for(b in proposed_b){
  print(b)
  rmse_results[,which(proposed_b == b)] <- replicate(
    nsims, loc_rmse_fullsim(b,mean_by_station,true_data)
  )
}
colnames(rmse_results) <- proposed_b
### ------------------- ####

### SUMMARIZE RESULTS

par(mfrow=c(1,3))
for(b in proposed_b){
  hist(rmse_results[,as.character(b)],main=paste("B=",b),xlab="RMSE")
  abline(v=mean(rmse_results[,as.character(b)]),col="red")
}

apply(rmse_results,2,mean,na.rm=TRUE)
apply(rmse_results,2,max,na.rm=TRUE)
apply(rmse_results,2,min,na.rm=TRUE)
apply(rmse_results,2,var,na.rm=TRUE)

#saveRDS(rmse_results, paste0("outputs/rmse_loc",today,".rds"))

## Q???


# EXAMPLE ERROR
# [1] 0
# model   psill      range
# 1   Exp 35.4051 -0.8453095
# Error in load.variogram.model(object$model[[name]], c(i - 1, i - 1), max_dist = max_dist) : 
#   variogram range can never be negative
