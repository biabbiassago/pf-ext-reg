
library(geoR)
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

# param_variog <- function(emp_variog){
#   # from the empirical semi-variogram, get starting initial values
#   # for fitting the parametric semi-variogram
# 
#   init_nugget <- emp_variog[1,"gamma"]
#   init_range <- emp_variog[which.max(emp_variog$gamma),"dist"]
#   init_psill <- max(emp_variog$gamma) - init_nugget
# 
#   param_variog <- gstat::fit.variogram(
#     emp_variog , model=gstat::vgm(psill=init_psill,"Exp",range=init_range,nugget=init_nugget),
#     fit.method=2
#   )
# 
#   return(param_variog)
# }

predict_location <- function(pred_df){
  # use kriging to predict
  
  coordinates(pred_df) <- ~ x_coords + y_coords
  #loc_vg <- variogram(loc_par~1, pred_df) # calculates sample variogram values
  
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
  loc_predicted$predict[pred_df$station] <- pred_df$loc_par
  return(loc_predicted$predict)
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
      loc_par_pred = loc_predicted,
      scale_par = full_parms_est["scale",],
      shape = full_parms_est['shape',],
      sampled = ifelse(station %in% cur_sample$stations,1,0)
    )
  return(df)
}
plot_actvpred_loc <- function(actvpred_df,b){
  # takes a df in format actvpred and makes plot to show actual vs predicted and sampled 
  # units
  p1 <- ggplot(actvpred_df, aes(x = x_coords, y = y_coords, fill = loc_par_pred)) + 
    geom_raster() +
    xlab("x") + 
    ylab("y") + 
    scale_fill_viridis_c(option="viridis",direction=-1,name="mu") +
    theme_bw() + 
    ggtitle(paste("Interpolated mu for B=",b))
  
  # p2 <- ggplot(actvpred_df, aes(x=x_coords, y=y_coords, fill = loc_par_est)) +   geom_raster() +
  #   xlab("x") + 
  #   ylab("y") + 
  #   scale_fill_viridis_c(option="viridis",direction=-1) +
  #   theme_bw()
  # 
  # p3 <- ggplot(actvpred_df, aes(x=x_coords, y=y_coords, fill = sampled)) +   geom_raster() +
  #   xlab("x") + 
  #   ylab("y") + 
  #   theme_bw()
  
  return(p1)
}

rmse_actvpred_loc <- function(x){
  #"takes a df in format actvpred and returns for location RMSE
  sqrt(
    sum(
  (x$loc_par_est-x$loc_par_pred)^2)/N
  )
}
main_est_loc <- function(cur_sample, true_data, full_params_est){
  # main function to apply to each cur sample
  pred_df <- make_pred_df(cur_sample)
  loc_predicted <- predict_location(pred_df)
  actvpred_df <- make_results_df(true_data,cur_sample,loc_predicted,full_parms_est)
  return(
    list(
      "df" = actvpred_df,
      "rmse"= rmse_actvpred_loc(actvpred_df)
    )
  )
}


# plots

b5 <-  main_est_loc(xb5$sample_df, true_data, full_params_est)$df
pi5 <- plot_actvpred_loc(main_est_loc(xb5$sample_df, true_data, full_params_est)$df,b=5)
pi3 <- plot_actvpred_loc(main_est_loc(xb3$sample_df, true_data, full_params_est)$df,b=3)
pi1 <- plot_actvpred_loc(main_est_loc(xb1$sample_df, true_data, full_params_est)$df,b=1)
pi0 <- plot_actvpred_loc(main_est_loc(xb0$sample_df, true_data, full_params_est)$df,b=0)
pi_est <- ggplot(b5, aes(x=x_coords, y=y_coords, fill = loc_par_est)) +   geom_raster() +
  xlab("x") +
  ylab("y") +
  scale_fill_viridis_c(option="viridis",direction=-1,name="Est. mu") +
  ggtitle("Full Dataset") +
  theme_bw()

  
  (pi5 + pi3)/ (pi1 + pi0)/(plot_spacer() + pi_est)
  
  