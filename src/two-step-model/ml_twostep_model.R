library(sp)
library(geoR)
source(here::here("src/3.evr_bystation.R"))

two_step_model <- function(full_dat){
    tmp <- full_dat
    df <- full_dat$true_data
    true_coords <- cbind(unique(df$x_coords),unique(df$y_coords))
    pred_df <- make_pred_df(df)
    variog_res <- predict_location(pred_df,true_coords)
    mean_by_station <- make_mean_by_station_df(df)
    geo_df <- as.geodata(
      data.frame(
        y=variog_res,
        x=mean_by_station$x,
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
      print.pars=FALSE,
      message=FALSE
    )
    # get parameters of interest
    mu_est <- pred_df$loc_par
    mu_025 <- pred_df$loc_025ci
    mu_975 <- pred_df$loc_975ci
    
    alpha0_est <- lin_fit$beta[1]
    alpha0_low <- alpha0_est - qnorm(0.025,lower.tail=F)*sqrt(lin_fit$beta.var[1,1])
    alpha0_hi <-  alpha0_est + qnorm(0.025,lower.tail=F)*sqrt(lin_fit$beta.var[1,1])
    
    
    alpha1_est <- lin_fit$beta[2]
    alpha1_low <- alpha1_est - qnorm(0.025,lower.tail=F)*sqrt(lin_fit$beta.var[2,2])
    alpha1_hi <-  alpha1_est + qnorm(0.025,lower.tail=F)*sqrt(lin_fit$beta.var[2,2])
    
    
    beta0_est <- lin_fit$sigmasq
    beta0_low <- NA
    beta0_hi <- NA
    
    
    beta1_est <- lin_fit$phi
    beta1_low <- NA
    beta1_hi <- NA
    
    
    rslts <- as.list(numeric(7*4))
    dim(rslts) <- c(7,4)
    
    # mus
    rslts[[1,1]] <- mu_est
    rslts[[1,2]] <- mu_025
    rslts[[1,3]] <- mu_975
    rslts[[1,4]] <- tmp$true_mus
    
    # alpha0
    rslts[[2,1]] <- alpha0_est
    rslts[[2,2]] <- alpha0_low
    rslts[[2,3]] <- alpha0_hi
    rslts[[2,4]] <- 0
    
    #alpha1
    rslts[[3,1]] <- alpha1_est
    rslts[[3,2]] <- alpha1_low
    rslts[[3,3]] <- alpha1_hi
    rslts[[3,4]] <- tmp$data_parms$ALPHA1
    
    #beta0
    rslts[[4,1]] <- beta0_est
    rslts[[4,2]] <- beta0_low
    rslts[[4,3]] <- beta0_hi
    rslts[[4,4]] <- tmp$data_parms$SIGMA2ETA
      
    #beta1
    rslts[[5,1]] <- beta1_est
    rslts[[5,2]] <- beta1_low
    rslts[[5,3]] <- beta1_hi
    rslts[[5,4]] <- tmp$data_parms$PHIETA
      
    #sigma2 of the gev (meaningless in this case)
    rslts[[6,1]] <- NA
    rslts[[6,2]] <- NA
    rslts[[6,3]] <- NA
    rslts[[6,4]] <-  tmp$data_parms$SCALE
      
    #xi
    rslts[[7,1]] <- NA
    rslts[[7,2]] <- NA
    rslts[[7,3]] <- NA
    rslts[[7,4]] <-  tmp$data_parms$SHAPE
    
    colnames(rslts) <- c("est","95ci_low","95ci_hi","truth")
    rownames(rslts) <- c("mu","alpha0","alpha1","beta0","beta1","sigma2","xi")
    return(rslts)
}

#results <- two_step_model(full_dat)
