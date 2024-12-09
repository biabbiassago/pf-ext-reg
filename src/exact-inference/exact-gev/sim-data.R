library(SpatialExtremes)
library(fields)
library(geoR)
library(MASS)
library(spatstat)
source(here::here("src/rep-gev/reparametrized_gev.R"))



make_sim_data_max <-
  function(beta = 2,
           true_lambda_star = 300,
           true_eta = 2.5,
           true_nu = 1,
           true_sigma2_S = 3,
           true_rho_S = 0.2, 
           true_sigma2_T = 0.8, # will be ignored if `varying_range=F`
           true_rho_T = 0.3, # will be ignored if `varying_range=F`
           true_xi=0.1,
           varying_range = F # decides whether to add random effects to the range parameter
    ) {
    k_true <- c(rpois(1, true_lambda_star))
    true_pp <- spatstat.random::runifpoint(k_true)
    
    # simulate S_k
    dist.g <- fields::rdist(coords(true_pp))
    R_S <-
      fields::Matern(dist.g, range = true_rho_S, smoothness = 1)
    S_mat <- true_sigma2_S * R_S
    S_k_unordered <- mvrnorm(1, mu = rep(0, k_true), Sigma = S_mat)
    
    # Select which coordinates will be kept
    probs <- pnorm(beta * S_k_unordered / sqrt(true_sigma2_S))
    obs_coords <-
      coords(true_pp)[(sapply(probs, function(x) {
        rbinom(1, 1, x)
      }) == 1),]
    idx_keep <- as.integer(rownames(obs_coords))
    S_x <- S_k_unordered[idx_keep]
    n <- dim(obs_coords)[1]
    
    # some clean up (reorder coordinates and random effects as c(observed, discarded))
    # just to make it easier to keep track of things after.
    S_k <- c(S_x, S_k_unordered[!(1:k_true %in% idx_keep)])
    all_coords <-
      rbind(obs_coords, coords(true_pp)[!(1:k_true %in% idx_keep),])
    
    # Set the median parameter
    q_all <- true_eta + S_k
    q <- q_all[1:n]
    
    # Simulating the range parameter and the spatial random effects associated with it
    if(varying_range ==T){
      R_T <- fields::Matern(fields::rdist(all_coords),
                            range = true_rho_T,
                            smoothness = 1)
      T_mat <- true_sigma2_T * R_T
      T_k <- mvrnorm(1, mu = rep(0, k_true), Sigma = T_mat)
      log_r_all <- true_nu + T_k
    }
    else if(varying_range == F){
      log_r_all <- true_nu
      T_k <- NULL
    }

    # Shape parameter set in the arguments (true_xi)
    
    # Simulate data. Keep only selected coordinates
    y_all <- rgevrep(k_true, q_all, exp(log_r_all), true_xi)
    y <- y_all[1:n]
    
    return(
      list(
        y = y,
        obs_coords = obs_coords,
        littlen = littlen,
        y_all = y_all,
        true_q = q_all,
        true_log_r = log_r_all,
        true_S_k = S_k,
        true_T_k = T_k,
        all_coords = all_coords,
        params = c(
          "true_lambda_star"=true_lambda_star,
          "true_eta"=true_eta,
          "true_nu"=true_nu,
          "true_sigma2_S"=true_sigma2_S,
          "true_rho_S"=true_rho_S,
          "true_sigma2_T"=true_sigma2_T,
          "true_rho_T"=true_rho_T,
          "true_xi"=true_xi
        )
      )
    )
  }
