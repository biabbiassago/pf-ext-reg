" Simulate a Spatial Process for the max using a GEV "
library(tidyverse)
library(MASS)
library(fields)
library(extRemes)
library(VGAM)
source(here::here("src/utils.R"))

make_true_data <- function(n_stations,months_rep){
  sim_locs <- sim_random_locs(n_stations)
  x_coords <- sim_locs$x_coords
  y_coords <- sim_locs$y_coords
  distance_mat <- sim_locs$distance_mat
  X <- MASS::mvrnorm(
    1,
    mu = rep(1,n_stations),
    Sigma=make_Sigma_x(distance_mat,n_stations)
  )
  
  eta <- MASS::mvrnorm(
    1,
    mu = rep(0,n_stations),
    Sigma=make_Sigma_eta(distance_mat,n_stations)
  )
  
  mu <- ALPHA1 * X + eta
  
  nu <- MASS::mvrnorm(
    1,
    mu = rep(0,n_stations),
    Sigma=make_Sigma_nu(distance_mat,n_stations)
  )
  sigma2 <- exp(nu)
  
  df <- data.frame(
    apply(
      cbind(mu,sigma2),1,
      function(x) rgev(months_rep,location=x["mu"],scale=x["sigma2"],shape=SHAPE)
      )
   ) %>%
    pivot_longer(
      everything(),
      names_to = "station"
    ) %>%
    mutate(
      station=parse_number(station),
    ) %>%
    arrange(station) %>%
    mutate(
      measurements= rep(1:months_rep,n_stations),
      x = rep(X,each=months_rep),
      eta = rep(eta,each=months_rep),
      nu = rep(nu, each=months_rep),
      x_coords = rep(x_coords,each=months_rep),
      y_coords= rep(y_coords,each=months_rep),
    )
  
  data_params = list(
    "SIGMA2X"=SIGMA2X,
    "SIGMA2ETA"=SIGMA2ETA,
    "SIGMA2NU" = SIGMA2NU,
    "PHIX"=PHIX,
    "PHIETA"=PHIETA,
    "PHINU"=PHINU,
    "ALPHA0"=0,
    "ALPHA1"=ALPHA1,
    "SHAPE"=SHAPE
  )
  
  true_data = list(
    "true_data" = df, 
    "true_mus" = mu,
    "true_sigmas" = sigma2,
    "data_parms" = data_params
  )
  
  return(
    true_data  
  )
  if(save==TRUE){
    write_rds(
      true_data,
      file=here::here(paste0("data/sim-gev",now(),".rds"))
    )
  }
}
