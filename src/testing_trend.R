
source(here::here("src/utils.R"))
source(here::here("src/2.sample.R"))
source(here::here("src/3.evr_bystation.R"))

print(paste("Phi is",PHIX,"sample size is",SAMPLE_SIZE))

nsims <- 100
today <- Sys.Date()
proposed_b <- c(5,3,1,0)
SIM_TYPE <- "sim3-gev"
true_data <- read.csv(paste0("data/sim-gev-3-PHIX",PHIX,".csv"))
  
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
  fit<- fevd(stat_dat$value)
  return(fit$results$par)
}


#TRUE_MUS <-mu


make_pred_trend_df <- function(cur_sample){

  parms_est <- sapply(
    unique(cur_sample$station),
    function(x) fit_evr_by_station_trend(x,cur_sample$sample_df)
  )

  station_mus <- parms_est["location",]
  pred_df <- cur_sample$sample_df %>%
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
predict_location <- function(pred_df,krig=TRUE){
  # use kriging to predict
  
  coordinates(pred_df) <- ~ x_coords + y_coords
  
  
  # parametric variogram. Intial values are computed as to mimic the process of
  # "looking" at a semivariogram
  suppressWarnings({
    loc_vg <- geoR::variog(data=pred_df$loc_par,coords=pred_df@coords, messages=FALSE)
    loc_fit <- geoR::variofit(loc_vg, cov.model = "exp", fix.nugget = T,messages=FALSE)
    if(krig==TRUE){
      loc_predicted <- krige.conv(
        data = pred_df$loc_par,
        coords=pred_df@coords,
        locations = true_coords,
        krige = krige.control(type.krige = "OK",obj.model = loc_fit),
      )
    }
    else{
      loc_predicted <- NULL
    }

  })
  cov_pars <- loc_fit$cov.pars
  names(cov_pars) <- c("sigma2mu","phimu")
  # true value for the actually sampled stations
  return(list("predicted" = loc_predicted$predict, "cov_pars" = cov_pars))
}


get_regression_estimates <- function(b, mean_by_station, true_data){
  
  # take a sample and get sampled stations
  xb <- prefsamp_station(b=b, mean_by_station, true_data)
  sampled_stations <- sort(xb$stations)
  pred_df <- make_pred_trend_df(cur_sample = xb) # create df 
  
  # do kriging, using predict_location function
  # if full sample size no kriging.
  if(SAMPLE_SIZE != N){
    variog_res <- predict_location(pred_df)
    y <- variog_res$predicted
    y[sampled_stations] <- pred_df$loc_par
  }
  else{
    variog_res <- predict_location(pred_df,krig=FALSE)
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
    ini=c(SIGMA2X,PHIX),
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
      "phi_est" = phi_est,
      "variog_cov_pars" = variog_res$cov_pars
    )
  )
}


reg_results_beta0 <- matrix(NA,nrow=nsims,ncol=length(proposed_b))
reg_results_beta1 <- matrix(NA,nrow=nsims,ncol=length(proposed_b))
reg_results_sigma <- matrix(NA,nrow=nsims,ncol=length(proposed_b))
reg_results_phi <- matrix(NA,nrow=nsims,ncol=length(proposed_b))
reg_results_variogsigma2 <- matrix(NA,nrow=nsims,ncol=length(proposed_b))
reg_results_variogphi <- matrix(NA,nrow=nsims,ncol=length(proposed_b))

for(b in proposed_b){
  for(i in 1:nsims){
    print(paste0("For b=",b," iteration : ",i,"/",nsims))
    tmp_res <- get_regression_estimates(b,mean_by_station,true_data)
    reg_results_beta0[i,which(proposed_b == b)] <- tmp_res$beta0_est
    reg_results_beta1[i,which(proposed_b == b)] <- tmp_res$beta1_est
    reg_results_sigma[i,which(proposed_b == b)] <- tmp_res$sigma_est
    reg_results_phi[i,which(proposed_b == b)] <- tmp_res$phi_est
    reg_results_variogsigma2[i,which(proposed_b == b)] <- tmp_res$variog_cov_pars["sigma2mu"]
    reg_results_variogphi[i,which(proposed_b == b)] <- tmp_res$variog_cov_pars["phimu"]
  }
}

reg_results_byB <- lapply(list(
  "beta0"=reg_results_beta0,
  "beta1"=reg_results_beta1,
  "sigma"=reg_results_sigma,
  "phi"=reg_results_phi,
  "variogsigma"=reg_results_variogsigma2,
  "variogphi"=reg_results_variogphi
  ),
  FUN= function(x){colnames(x)<-proposed_b;x}
)


# #FULL SET
# SAMPLE_SIZE <- N
# no_sampling_results <- get_regression_estimates(b,mean_by_station,true_data)


saveRDS(
  reg_results_byB,
  paste0(
    "outputs/reg-est-loc-",
    today,"-",
    SIM_TYPE,
    "sampsize",
    SAMPLE_SIZE,
    "phi",
    PHIX,
    ".rds")
)

# avg_over_sims <- lapply(reg_results_byB, function(x){apply(x,2,mean,na.rm=TRUE)})
# var_over_sims <- lapply(reg_results_byB, function(x){apply(x,2,var,na.rm=TRUE)})
# 
# 
# 
# PLOTTING_PAR <- "beta1"
# par(mfrow=c(2,2))
# for(b in proposed_b){
#   hist(reg_results_byB[[PLOTTING_PAR]][,1],
#        main=paste0(PLOTTING_PAR," for b=",b),
#        xlab="beta1 est")
# }
# 
# 
# 




