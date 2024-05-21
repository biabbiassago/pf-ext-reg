" Simulate a Spatial Process for the max using a GEV "
library(sf)
library(tidyverse)
library(MASS)
library(fields)
source(here::here("src/utils.R"))
source(here::here("src/rep-gev/reparametrized_gev.R"))


make_true_data_pp <-
  function(n_stations,
           months_rep,
           shape,
           a,
           pref_on,
           dimxy = DIMXY,
           save = FALSE) {

  
    pref_field <- geoR::grf(
      10000,
      grid = "reg",
      nsim = 1,
      cov.model = "exponential",
      cov.pars = c(PREF_VAR,PREF_SCALE),
      message=FALSE
    )
    
    phi <- geoR::sample.geodata(
      pref_field,
      size=n_stations,
      prob=exp(pref_field$data)
    )

    
    DMat <-
      fields::rdist(phi$coords) # distance matrix between the points
    
    x_coords <- phi$coords[,"x"]
    y_coords <- phi$coords[,"y"]
    dgrid <- dimxy ^ 2

    
    locs_counts_mat <- quadratcount(as.ppp(phi$coords,W=owin()), nx = dimxy, ny = dimxy)
    locs_counts <- as.vector(t(locs_counts_mat))
    
    # this is the distance matrix that has both the centroids and 
    # the observed points, so that we can draw samples from the MVN normal
    # that includes all of these things. 
    # using SF package right now to do it but i think there is a better way that does
    # not rely on it. 
    q <- st_as_sf(as.ppp(phi$coords,W=owin())) %>% filter(label != "window")
    ysf <- st_make_grid(q, n = c(dimxy, dimxy)) %>% st_sf() %>%
      mutate(id = assign_regular_ids(dimxy)) %>% arrange(id)
    DFullMat <- st_distance(c(st_centroid(ysf)$geometry,q$geom))
    
    eta <- MASS::mvrnorm(1,
                         mu = rep(0, n_stations),
                         Sigma = make_Sigma_eta(DMat, n_stations))
    
    nu <- MASS::mvrnorm(1,
                        mu = rep(0, n_stations),
                        Sigma = make_Sigma_nu(DMat, n_stations))
    if (pref_on == "median") {
      qmedian <- ALPHA0 + eta + a * phi$data
      iqr <- exp(RHO0 + nu)
      rep_pars <- as.matrix(data.frame(qmedian, iqr, shape))
      gev_pars <- as.matrix(data.frame(rep2gev(qmedian, iqr, shape)))
    }
    else if (pref_on == "loc") {
      mu <- ALPHA0 + eta + a * phi$data
      sigma <- exp(RHO0 + nu)
      rep_pars <- as.matrix(data.frame(gev2rep(mu, sigma, shape)))
      gev_pars <- as.matrix(data.frame(mu, sigma, shape))
    }
    else{
      stop("pref_on (Where preferential sampling is on) must be on `median` or `loc` only.")
    }
    
    
    df <- data.frame(apply(gev_pars, 1,
                           function(x)
                             rgev(
                               months_rep,
                               loc = x["mu"],
                               scale = x["sigma"],
                               shape = x["xi"]
                             ))) %>%
      pivot_longer(everything(),
                   names_to = "station") %>%
      mutate(station = parse_number(station),) %>%
      arrange(station) %>%
      mutate(
        measurements = rep(1:months_rep, n_stations),
        eta = rep(eta, each = months_rep),
        nu = rep(nu, each = months_rep),
        x_coords = rep(x_coords, each = months_rep),
        y_coords = rep(y_coords, each = months_rep),
      )
    
    data_params = list(
      "SIGMA2ETA" = SIGMA2ETA,
      "SIGMA2NU" = SIGMA2NU,
      "PHIETA" = PHIETA,
      "PHINU" = PHINU,
      "ALPHA0" = ALPHA0,
      "RHO0" = RHO0,
      "PREF_VAR"=PREF_VAR,
      "PREF_SCALE"=PREF_SCALE,
      "s" = n_stations,
      "dgrid" = dgrid,
      "months_rep" = months_rep,
      "locs_counts" = locs_counts,
      "DIMXY" = dimxy
    )
    
    distance_mats = list("DMat" = DMat,
                         "DFullMat " = DFullMat)
    
    true_data = list(
      "true_data" = df,
      "true_medians" = rep_pars[, "qmedian"],
      "true_iqrs" = rep_pars[, "iqr"],
      "true_mus" = gev_pars[, "mu"],
      "true_sigmas" = gev_pars[, "sigma"],
      "true_xi" = shape,
      "true_phi" = phi,
      "data_parms" = data_params,
      "distance_mats" = distance_mats
    )
    
    return(true_data)
    if (save == TRUE) {
      write_rds(true_data,
                file = here::here(paste0(
                  "data/sim-gev/gen_pp_", now(), ".rds"
                )))
    }
  }

## auxiliary function
assign_regular_ids <- function(dimxy) {
  start <- dimxy * dimxy
  idx <- c()
  for (iter in seq(0, start - dimxy, dimxy)) {
    idx = append(idx, ((start - iter) - (dimxy - 1)):(start - iter))
  }
  return(idx)
}

