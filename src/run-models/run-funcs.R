make_stan_data <- function(dat) {
  return(
    list(
      s = dat$data_parms$s,
      dgrid = dat$data_parms$dgrid,
      N = (dat$data_parms$s) * (dat$data_parms$months_rep),
      cell_size = 1 / ((dat$data_parms$DIMXY) ^ 2),
      grid_idx = dat$data_parms$station_grid_id,
      months = dat$data_parms$months_rep,
      z = dat$true_data$value,
      DMat = dat$distance_mats$DMat,
      DFullMat = dat$distance_mats$DFullMat,
      locs = dat$data_parms$locs_counts
    )
  )
}

crude_init_values <- function(fd, s, chains) {
  mle <- SpatialExtremes::gevmle(fd$true_data$value)
  mle_loc <- mle["loc"]
  mle_scale <- mle["scale"]
  mle_shape <- mle["shape"]
  vals <-
    list(
      "beta0" = 0.5,
      "beta1" = 0.5,
      "gamma0" = 0.1,
      "gamma1" = 0.1,
      "eta" = rep(0, s),
      "nu" = rep(0, s),
      "alpha0" = mle_loc,
      "rho0" = mle_scale,
      "xi" = mle_shape
    )
  init_vals <- replicate(chains, vals, FALSE)
  return(init_vals)
}

prepare_stan_data_no_pref <- function(sim_dat) {
  #x_mat <-
  #cbind(rep(1, length(sim_dat$true_mus)), unique(sim_dat$true_data$x))
  
  unique_coords <- sim_dat$true_data %>%
    group_by(station) %>%
    summarize(x_coords = first(x_coords),
              y_coords = first(y_coords))
  
  loc <- cbind(unique_coords$x_coords, unique_coords$y_coords)
  
  distance_mat <- as.matrix(dist(loc))
  
  s <- length(sim_dat$true_mus)
  months <- dim(sim_dat$true_data)[1] / length(sim_dat$true_mus)
  
  stan_data <- list(
    s = s,
    months = months,
    N = s * months,
    z = sim_dat$true_data$value,
    #x = x_mat,
    DMat = distance_mat
  )
  return(stan_data)
}