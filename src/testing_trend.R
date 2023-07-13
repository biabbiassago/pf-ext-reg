

source(here::here("src/utils.R"))
source(here::here("src/sim-gev/sim-timeconst.R"))
source(here::here("src/2.sample.R"))
source(here::here("src/3.evr_bystation.R"))

nsims <- 100
today <- Sys.Date()
proposed_b <- c(5,3,1,0)
SIM_TYPE <- "sim2-gev"

true_data <- read.csv(here::here("data/sim-gev-2.csv"))

mean_by_station <- true_data %>%
  group_by(station) %>%
  summarize(
    mean_value = mean(value),
    x_coords = first(x_coords),
    y_coords = first(y_coords),
    x = first(x)
)


fit_evr_by_station_trend <- function(i,dat){
  
  stat_dat <- dat %>% filter(station==i)
  fit<- fevd(value, location.fun = ~ 1, data=stat_dat)
  return(fit$results$par)
}


TRUE_MUS <- mu


# sim_beta_trend <- function(b, mean_by_station, true_data){
#   # given a data-set and a b "strength of preferential sampling - run"
#   xb <- prefsamp_station(b,mean_by_station, true_data)
#   pred_df <- make_pred_trend_df(xb$sample_df)
#   y <- predict_location(pred_df)
#   lin_fit <- lm(y~X)
#   beta_est <- lin_fit$coefficients[2]
#   return(beta_est)
# }
# 



# 
make_pred_trend_df <- function(cur_sample){

  parms_est <- sapply(
    unique(cur_sample$station),
    function(x) fit_evr_by_station_trend(x,cur_sample)
  )

  station_mus <- parms_est["location",]
  pred_df <- cur_sample %>%
    group_by(station) %>%
    summarise(
      mean_value = mean(value),
      station=first(station),
      x_coords=first(x_coords),
      y_coords=first(y_coords),
    ) %>% mutate(
      loc_par = station_mus
    )
  return(pred_df)
}
predict_location <- function(pred_df){
  # use kriging to predict
  
  coordinates(pred_df) <- ~ x_coords + y_coords
  
  # TODO : better? e.g. predict fit.variogram parameters using the likfit function to adapt
  
  # parametric variogram. Intial values are computed as to mimic the process of
  # "looking" at a semivariogram
  suppressMessages({
    loc_vg <- geoR::variog(data=pred_df$loc_par,coords=pred_df@coords)
    loc_fit <- geoR::variofit(loc_vg, cov.model = "exp", fix.nugget = T)
    
    loc_predicted <- krige.conv(
      data = pred_df$loc_par,
      coords=pred_df@coords,
      locations = true_coords,
      krige = krige.control(type.krige = "OK",obj.model = loc_fit)
    )
  })
  # true value for the actually sampled stations
  return(loc_predicted$predict)
}



SAMPLE_SIZE <- 100
get_regression_estimates <- function(b, mean_by_station, true_data){
  
  # take a sample and get sampled stations
  xb <- prefsamp_station(b=b, mean_by_station, true_data)
  sampled_stations <- sort(xb$stations)
  pred_df <- make_pred_trend_df(cur_sample = xb$sample_df) # create df 
  
  # do kriging, using predict_location function
  # if full sample size no kriging.
  if(SAMPLE_SIZE != N){
    y <- predict_location(pred_df)
    y[sampled_stations] <- pred_df$loc_par
  }
  else{
    y <- pred_df$loc_par
  }
  
  #fit ML regression
  ## start creating dataset as needed
  geo_df <- as.geodata(
    data.frame(
      y=y,
      x=X,
      x_coords=mean_by_station$x_coords,
      y_coords=mean_by_station$y_coords),
    data.col=c(1,2),
    coords.col = c(3,4)
  )
  #use likfit to get ML estimates
  lin_fit <- likfit(
    geo_df,
    coords = geo_df$coords,
    data=geo_df$data[,1],
    cov.model="exponential",
    trend = ~ geo_df$data[,2],
    ini=c(SIGMA2,PHI),
    fix.nugget = TRUE,
    message=FALSE
  )
  # get parameters of interest
  beta0_est <- lin_fit$beta[1]
  beta1_est <- lin_fit$beta[2]
  sigma_est <- lin_fit$sigmasq
  phi_est <- lin_fit$phi
  
  return(
    list(
      "beta0_est" = beta0_est,
      "beta1_est" = beta1_est,
      "sigma_est" = sigma_est,
      "phi_est" = phi_est
    )
  )
}


SAMPLE_SIZE <- 100
reg_results_beta0 <- matrix(NA,nrow=nsims,ncol=length(proposed_b))
reg_results_beta1 <- matrix(NA,nrow=nsims,ncol=length(proposed_b))
reg_results_sigma <- matrix(NA,nrow=nsims,ncol=length(proposed_b))
reg_results_phi <- matrix(NA,nrow=nsims,ncol=length(proposed_b))
for(b in proposed_b){
  for(i in 1:nsims){
    print(paste0("For b=",b," iteration : ",i,"/",nsims))
    tmp_res <- get_regression_estimates(b,mean_by_station,true_data)
    reg_results_beta0[i,which(proposed_b == b)] <- tmp_res$beta0_est
    reg_results_beta1[i,which(proposed_b == b)] <- tmp_res$beta1_est
    reg_results_sigma[i,which(proposed_b == b)] <- tmp_res$sigma_est
    reg_results_phi[i,which(proposed_b == b)] <- tmp_res$phi_est
  }
}

reg_results_byB <- lapply(list(
  "beta0"=reg_results_beta0,
  "beta1"=reg_results_beta1,
  "sigma"=reg_results_sigma,
  "phi"=reg_results_phi
  ),
  FUN= function(x){colnames(x)<-proposed_b;x}
)

saveRDS(
  reg_results_byB,
  paste0("outputs/reg-est-loc-",today,"-",SIM_TYPE,".rds")
)

avg_over_sims <- lapply(reg_results_byB, function(x){apply(x,2,mean,na.rm=TRUE)})
var_over_sims <- lapply(reg_results_byB, function(x){apply(x,2,var,na.rm=TRUE)})


#FULL SET
SAMPLE_SIZE <- 225
no_sampling_results <- get_regression_estimates(b,mean_by_station,true_data)
names(no_sampling_results) <- names(reg_results_byB)
SAMPLE_SIZE <- 100


PLOTTING_PAR <- "beta1"
par(mfrow=c(2,2))
for(b in proposed_b){
  hist(reg_results_byB[[PLOTTING_PAR]][,1],
       main=paste0(PLOTTING_PAR," for b=",b),
       xlab="beta1 est")
}

