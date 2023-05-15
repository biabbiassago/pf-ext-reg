
library(patchwork)
library(sp)
library(gstat)
source("src/utils.R")

# true_data <- read.csv("data/sim1.csv")
# cur_sample <- xb3$sample_df

fit_evr_by_station <- function(i,dat){
  
  # fit a EVR to each station with no predictors yet.
  stat_dat <- dat %>% filter(station==i)
  fit<- fevd(value~1,data=stat_dat)
  # location: 8.166, scale = 1.05, shape = -0.117
  return(fit$results$par)
}

make_pred_df <- function(cur_sample){
  parms_est <- sapply(
    unique(cur_sample$station),
    function(x) fit_evr_by_station(x,cur_sample)
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
      loc_par = parms_est["location",],
      scale_par = parms_est["scale",],
      shape = parms_est['shape',]
    )
  return(pred_df)
}
predict_location <- function(pred_df){
  coordinates(pred_df) <- ~ x_coords + y_coords
  # use kriging to predict
  loc_vg <- variogram(loc_par~1, pred_df) # calculates sample variogram values
  
  # TODO : predict fit.variogram parameters using the likfit function to adapt
  
  
  loc_fit <- fit.variogram(loc_vg, model=vgm(99, "Exp",0.1)) # fit model
  # initial value for s can be cov(loc_par,loc_par)
  
  coordinates(coords) <- ~Var1 + Var2
  loc_predicted <- krige(loc_par ~ 1, pred_df,coords, model=loc_fit)
  loc_predicted[pred_df$station,"var1.pred"] <- pred_df$loc_par
  return(loc_predicted)
}
estimate_parms <- function(true_data){
  sapply(
  unique(true_data$station),
  function(x) fit_evr_by_station(x,true_data)
  )
}
make_results_df <- function(true_data, cur_sample, loc_predicted, full_parms_est){
  df<- true_data %>%
    group_by(station) %>%
    summarise(
      mean_value = mean(value),
      station=first(station),
      x_coords=first(x_coords),
      y_coords=first(y_coords),
    ) %>% mutate(
      loc_par_est = full_parms_est["location",],
      loc_par_pred = loc_predicted$var1.pred,
      scale_par = full_parms_est["scale",],
      shape = full_parms_est['shape',],
      sampled = ifelse(station %in% cur_sample$stations,1,0)
    )
  return(df)
}
plot_actvpred_loc <- function(actvpred_df){
  # takes a df in format actvpred and makes plot to show actual vs predicted and sampled 
  # units
  p1 <- ggplot(actvpred_df, aes(x = x_coords, y = y_coords, fill = loc_par_est)) + 
    geom_raster() +
    xlab("x") + 
    ylab("y") + 
    scale_fill_viridis_c(option="viridis",direction=-1) +
    theme_bw()
  
  p2 <- ggplot(actvpred_df, aes(x=x_coords, y=y_coords, fill = loc_par_pred)) +   geom_raster() +
    xlab("x") + 
    ylab("y") + 
    scale_fill_viridis_c(option="viridis",direction=-1) +
    theme_bw()
  
  p3 <- ggplot(actvpred_df, aes(x=x_coords, y=y_coords, fill = sampled)) +   geom_raster() +
    xlab("x") + 
    ylab("y") + 
    theme_bw()
  
  p1+p2+p3
}

rmse_actvpred_loc <- function(x){
  #"takes a df in format actvpred and returns for location RMSE
  sum(
  (x$loc_par_est-x$loc_par_pred)^2)/N
}
main_est_loc <- function(cur_sample, true_data){
  # main function to apply to each cur sample
  pred_df <- make_pred_df(cur_sample)
  loc_predicted <- predict_location(pred_df)
  full_parms_est <- estimate_parms(true_data)
  actvpred_df <- make_results_df(true_data,cur_sample,loc_predicted,full_parms_est)
  return(
    list(
      "df" = actvpred_df,
      "rmse"= rmse_actvpred_loc(actvpred_df)
    )
  )
}

