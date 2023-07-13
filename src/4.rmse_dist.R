library(extRemes)

SIM_TYPE <- "sim3"

today <- Sys.Date()

source(here::here("src/utils.R"))
source(here::here("src/2.sample.R"))
source(here::here("src/3.evr_bystation.R"))

nsims <- 100
proposed_b <- c(5,3,1,0)

loc_rmse_fullsim <- function(b, mean_by_station, true_data){
  # given a data-set and a b "strength of preferential sampling - run"
  full_parms_est <- estimate_parms(true_data)
  xb <- prefsamp_station(b,mean_by_station, true_data)
  evr_pred_df <- main_est_loc(xb$sample_df,true_data,full_parms_est)
  return(evr_pred_df$rmse)
}

### ------------------- ####
#### RUN SIM #####
true_data <- read.csv(paste0("data/",SIM_TYPE,".csv"))

mean_by_station <- true_data %>%
  group_by(station) %>%
  summarize(
    mean_value = mean(value),
    x_coords = first(x_coords),
    y_coords = first(y_coords),
    x = mean(x)
)


rmse_results <- matrix(NA,nrow=nsims,ncol=length(proposed_b))
for(b in proposed_b){
  print(b)
  rmse_results[,which(proposed_b == b)] <- replicate(
    nsims, loc_rmse_fullsim(b,mean_by_station,true_data)
  )
}
colnames(rmse_results) <- proposed_b
### ------------------- ####

### SUMMARIZE RESULTS

par(mfrow=c(2,2))
for(b in proposed_b){
  hist(
    rmse_results[,as.character(b)],
    main=paste("B=",b),
    xlab="RMSE"
  )
}

meanres <- apply(rmse_results,2,mean,na.rm=TRUE)
maxres <- apply(rmse_results,2,max,na.rm=TRUE)
minres <- apply(rmse_results,2,min,na.rm=TRUE)
varres<- apply(rmse_results,2,var,na.rm=TRUE)

resrmse <- data.frame(meanres,maxres,minres,varres)

#kableExtra::kbl(resrmse,col.names=c("Mean","Max","Min","Var"),format="latex",booktabs=TRUE,digits = 3)
#kableExtra::kbl(resrmse,col.names=c("Mean","Max","Min","Var"),format="html",booktabs=TRUE,digits = 3)

resrmse %>%
  kableExtra::kbl(
    col.names=c("Mean","Max","Min","Var"),
    booktabs=TRUE, 
    digits = 3,
    caption= "Mean, Max, Min and Var of the RMSE by value of b, the parameter that controls preferential sampling",
    label="rmse",
    format="latex") %>% 
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")


saveRDS(rmse_results, paste0("outputs/rmse_loc-",today,"-",SIM_TYPE,".rds"))

#rmse_results <- readRDS("outputs/rmse_loc2023-06-14.rds")
