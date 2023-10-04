# PLOTTING
library(patchwork)
library(tidyverse)


SAMPLE_SIZE <- 280
PHIX <- 0.2

source(here::here("src/2.sample.R"))
source(here::here("src/utils.R"))
source(here::here("src/3.evr_bystation.R"))


sims <- read_rds(paste0("data/sim-gev-4-PHIX",PHIX,".rds"))
true_data <- sims$true_data
true_mus <- sims$true_mus
mean_by_station <- true_data %>%
  group_by(station) %>%
  summarize(
    mean_value = mean(value),
    x_coords = first(x_coords),
    y_coords = first(y_coords),
    x = first(x),
    eta = first(eta)
  )


fit_evr_by_station_trend <- function(i,dat){
  
  stat_dat <- dat %>% filter(station==i)
  fit <- fevd(stat_dat$value, method="GMLE")
  return(fit$results$par)
}
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
    
    
    value_mean_model <- likfit(
      coords = coordinates(pred_df),
      data=pred_df$mean_value,
      cov.model="exponential",
      ini=c(SIGMA2,PHI),
      fix.nugget = TRUE,
      message=FALSE
    )
    
    
    # 
    # loc_vg <- geoR::variog(data=pred_df$loc_par,coords=pred_df@coords, messages=FALSE)
    # loc_fit <- geoR::variofit(loc_vg, cov.model = "exp", fix.nugget = T,messages=FALSE)
    
    
    if(krig==TRUE){
      loc_predicted <- krige.conv(
        data = pred_df$loc_par,
        coords=pred_df@coords,
        locations = true_coords,
        krige = krige.control(type.krige = "OK",obj.model = value_mean_model),
      )
    }
    else if (krig==FALSE){
      loc_predicted <- NULL
    }
    
  })
  cov_pars <- value_mean_model$cov.pars
  names(cov_pars) <- c("sigma2mu","phimu")
  # true value for the actually sampled stations
  return(list("predicted" = loc_predicted$predict, "cov_pars" = cov_pars))
}

est_mus <- sapply(1:N, function(x) fit_evr_by_station_trend(x,dat=true_data)["location"])


kriged_mus <- matrix(NA, nrow=400, ncol=3)
colnames(kriged_mus) <- c("3","1","0")



in_sample_ind <- data.frame("3"=rep(NA,400),"1"=rep(NA,400),"0"=rep(NA,400))


for(B in c(3,1,0)){
  print(B)
  xb <- prefsamp_station(b=B,mean_by_station,true_data)
  pred_df <- make_pred_trend_df(cur_sample = xb) # create df 
  variog_res <- predict_location(pred_df)
  kriged_mus[,as.character(B)]<- variog_res$predicted
  kriged_mus[xb$stations,as.character(B)] <- pred_df$loc_par
  in_sample_ind[,as.character(B)] <- factor(sapply(1:N,function(x){ifelse(x %in% xb$stations,1,0)}),levels=c(0,1))
}



df_plot <- data.frame(
  true_mus,
  est_mus,
  "krig_mus3" = kriged_mus[,"3"],
  "in_sample_ind3"= in_sample_ind[,"3"],
  "krig_mus1" = kriged_mus[,"1"],
  "in_sample_ind1"= in_sample_ind[,"1"],
  "krig_mus0" = kriged_mus[,"0"],
  "in_sample_ind0"= in_sample_ind[,"0"],
  "x" = mean_by_station$x,
  "eta" = mean_by_station$eta,
  "x_coords"=mean_by_station$x_coords,"y_coords"=mean_by_station$y_coords)

rownames(df_plot) <- 1:N


p1 <- ggplot(df_plot, aes(x=x_coords, y=y_coords, fill = true_mus)) + 
  geom_raster() +
  xlab("x") +
  ylab("y") +
  scale_fill_viridis_c(option="viridis",direction=-1,name="True Mu Value",limits=c(-2,3)) +
  ggtitle(paste0("PHIX:",PHIX," true values of mu")) + 
  theme_bw()

p2 <- ggplot(df_plot, aes(x=x_coords, y=y_coords, fill = est_mus)) + 
  geom_raster() +
  xlab("x") +
  ylab("y") +
  scale_fill_viridis_c(option="viridis",direction=-1,name="Est Mu Value (GEV reg)",limits=c(-2,3)) +
  ggtitle(paste0("PHIX:",PHIX," estimated mu (full sample)" )) + 
  theme_bw()



kriged_plots <- function(B){
  p3 <- ggplot(df_plot, aes(x=x_coords, y=y_coords, fill = kriged_mus[,B])) + 
    geom_raster() +
    xlab("x") +
    ylab("y") +
    scale_fill_viridis_c(option="viridis",direction=-1,name="kriged mu Values",limits=c(-2,3)) +
    ggtitle(paste0("PHIX:",PHIX," kriged mu, b=",B)) + 
    theme_bw()
  
  p4 <- ggplot(df_plot, aes(x=x_coords, y=y_coords, fill = in_sample_ind[,B])) + 
    geom_raster() +
    xlab("x") +
    ylab("y") +
    scale_fill_grey(start=1,end=0.2,name="In Sampled Ind") + 
    ggtitle(paste0("Indicator if in sample, b=",B, "SAMPLE SIZE = ", SAMPLE_SIZE)) + 
    theme_bw()
  
  return(list(p3,p4))
  
}

plots3 <- kriged_plots("3")
plots1 <- kriged_plots("1")
plots0 <- kriged_plots("0")

p9 <- ggplot(df_plot, aes(x=x_coords, y=y_coords, fill = x)) + 
  geom_raster() +
  xlab("x") +
  ylab("y") +
  scale_fill_viridis_c(option="inferno",direction=-1,name="True X Values") +
  ggtitle(paste0("PHIX:",PHIX," Values of X")) + 
  theme_bw()

p10 <-ggplot(df_plot, aes(x=x_coords, y=y_coords, fill = mean_by_station$eta)) + 
  geom_raster() +
  xlab("x") +
  ylab("y") +
  scale_fill_viridis_c(option="inferno",direction=-1,name="True Eta Values") +
  ggtitle(paste0("PHIX:",PHIX," Values of Eta")) + 
  theme_bw()


final_plot <- (p1 + p2)  / (p9 + p10) / (plots3[[1]] + plots3[[2]]) / (plots1[[1]] + plots1[[2]]) / (plots0[[1]] + plots0[[2]])


ggsave(filename = paste0("plot0801-PHIX",PHIX,"SAMPLESIZE",SAMPLE_SIZE,".png"), width = 9, height=12, units = "in")
