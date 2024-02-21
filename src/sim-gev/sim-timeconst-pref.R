" Simulate a Spatial Process for the max using a GEV with Preferential Sampling"
library(tidyverse)
library(MASS)
library(fields)
library(extRemes)
library(VGAM)
source(here::here("src/utils.R"))
source(here::here("src/sampling_probs.R"))

make_true_data_pref <- function(n_stations,months_rep,b=1,save=FALSE){
  # sim_locs <- sim_grid(0.05)
  # grid_dim <- length(sim_locs$x_coords)
  # x_coords <- sim_locs$x_coords
  # y_coords <- sim_locs$y_coords
  # 
  # distance_mat <- sim_locs$distance_mat
  
  INIT_LOCS <- 2500
  sim_locs <- sim_random_locs(INIT_LOCS)
  x_coords <- sim_locs$x_coords
  y_coords <- sim_locs$y_coords
  distance_mat <- sim_locs$distance_mat
  
  
  X <- MASS::mvrnorm(
    1,
    mu = rep(1,INIT_LOCS),
    Sigma=make_Sigma_x(distance_mat,INIT_LOCS)
  )
  
  eta <- MASS::mvrnorm(
    1,
    mu = rep(0,INIT_LOCS),
    Sigma=make_Sigma_eta(distance_mat,INIT_LOCS)
  )
  
  mu <- ALPHA1 * X + eta

  df <- data.frame(
    sapply(
      mu,
      function(x) rgev(months_rep,location=x,scale=SCALE,shape=SHAPE))
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
      measurements= rep(1:months_rep,INIT_LOCS),
      x = rep(X,each=months_rep),
      eta = rep(eta,each=months_rep),
      x_coords = rep(x_coords,each=months_rep),
      y_coords= rep(y_coords,each=months_rep),
  )
  
  df_mean <- make_mean_by_station_df(df)
  
  probs <- generate_sampling_probs(b=b,df_mean,type = "quadratic")
  s_i <- sample(
    1:INIT_LOCS,
    size=n_stations,
    replace=FALSE,
    prob=probs$pi
  )
  sample_df <- df %>% filter(station %in% s_i)

  data_params = list(
    "PHIX" = PHIX,
    "PHIETA" = PHI,
    "SIGMA2X" = SIGMA2X,
    "SIGMA2ETA" = SIGMA2,
    "ALPHA0" = 0,
    "ALPHA1" = ALPHA1,
    "SCALE" = SCALE,
    "SHAPE" = SHAPE
  )
  
  true_data = list(
    "true_data" = sample_df, 
    "true_mus" = mu[s_i],
    "sampled_stations"=s_i,
    "data_parms" = data_params
  )
  
  
  if(save==TRUE){
    write_rds(
      true_data,
      file=here::here(paste0("data/sim-gev-pref",now(),".rds"))
    )
  }
  return(
    true_data  
  )
}
