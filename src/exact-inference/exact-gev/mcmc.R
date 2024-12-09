source(here::here("src/exact-inference/exact-gev/samplers.R"))
source(here::here("src/exact-inference/exact-gev/constants.R"))

gev_exact_mcmc <- function(y, obs_coords, nsims, out_file_loc) {
  # Initialize Chains
  littlen <- dim(obs_coords)[1]
  lambda_star <- c(100)
  k_star <- c(rpois(1, lambda_star[1]))
  pp <- spatstat.random::runifpoint(k_star[1])
  all_coords <- rbind(obs_coords, coords(pp))
  k <- c(dim(all_coords)[1])
  beta <- c(0)
  beta_acc_rate <- c(1)
  ## GEV Level 2 pars
  eta <- c(median(y))
  eta_acc_rate <- c(1)
  nu <- c(log(IQR(y)))
  nu_acc_rate <- c(1)
  xi <- c(0.5)
  xi_acc_rate <- c(1)
  ## GEV Random Effects on Median
  R <-
    fields::Matern(fields::rdist(all_coords),
                   range = 0.3,
                   smoothness = 1)
  cov_S_k <- 2 * R
  R_inv <- solve(R)
  S_k <-
    pracma::sqrtm(cov_S_k)$B %*% rnorm(dim(all_coords)[1], 0, 1)
  S_n <- matrix(NA, nrow = littlen, ncol = nsims)
  S_n[, 1] <- S_k[1:littlen]
  ## GEV Level 3 pars
  sigma2_S <- c(2)
  sigma2_S_acc_rate <- c(1)
  rho_S <- c(0.3)
  rho_S_acc_rate <- c(1)
  
  
  
  #### START RUNNING THE MCMC ITERATIONS #####
  for (i in 2:nsims) {
    #print(i)
    if ((i %% 1000) == 0) {
      cur_info_print(i)
      tmp <-   list(
        lambda_star = lambda_star,
        k = k,
        eta = eta,
        nu = nu,
        xi = xi,
        beta = beta,
        S_n = S_n,
        sigma2_S = sigma2_S,
        rho_S = rho_S,
        sim_data = sim_data
      )
      saveRDS(tmp, here::here("outputs/mcmc-exact/tmp.rds"))
    }
    
    # Step 1 : Sample Lambda Star
    lambda_star[i] <-
      sample_lambda_star(
        k[i - 1],
        lambda_prior_a = PRIOR_LAMBDA_A,
        lambda_prior_b = PRIOR_LAMBDA_B,
        area_B = 1
      )
    
    # Step 2: Simulate Discarded Locations
    all_coords_prev <- all_coords
    all_coords <-
      sample_all_coords(lambda_star[i],
                        beta[i - 1],
                        sigma2_S[i - 1],
                        rho_S[i - 1],
                        S_k,
                        obs_coords,
                        all_coords,
                        area_B)
    k[i] <- dim(all_coords)[1]
    
    #  # Step 3: sample the random effects S_x
    S_k <-
      sample_Sk_ess(
        S_k,
        beta[i - 1],
        sigma2_S[i - 1],
        rho_S[i - 1],
        eta[i - 1],
        nu[i - 1],
        xi[i - 1],
        y,
        littlen,
        lambda_star[i],
        obs_coords,
        all_coords_prev,
        all_coords,
        area_B = 1
      )
    S_n[, i] <- S_k[1:littlen]
    R_S <-
      fields::Matern(fields::rdist(all_coords),
                     range = rho_S[i - 1],
                     smoothness = 1)
    R_S_inv <- chol2inv(chol(R_S))
    
    ## Step 4: GEV level 2 parameters
    
    # Median Parameter eta
    eta_tmp <-
      sample_eta(eta[i - 1],
                 nu[i - 1],
                 xi[i - 1],
                 S_n[, i],
                 y,
                 prior_eta_mean = PRIOR_ETA_MEAN,
                 prior_eta_var = PRIOR_ETA_VAR)
    eta[i] <- eta_tmp$eta
    eta_acc_rate[i] <- eta_tmp$eta_acc
    
    # Scale parameter Xi
    nu_tmp <-
      sample_nu(nu[i - 1],
                eta[i - 1],
                xi[i - 1],
                S_n[, i],
                y,
                prior_nu_mean = PRIOR_NU_MEAN,
                prior_nu_var = PRIOR_NU_VAR)
    nu[i] <- nu_tmp$nu
    nu_acc_rate[i] <- nu_tmp$nu_acc
    
    # Shape parameter Xi
    xi_tmp <-
      sample_xi(xi[i - 1],
                eta[i],
                nu[i],
                S_n[, i],
                y,
                prior_xi_mean = PRIOR_XI_MEAN,
                prior_xi_var = PRIOR_XI_VAR)
    xi[i] <- xi_tmp$xi
    xi_acc_rate[i] <- xi_tmp$xi_acc
    
    # Step 5 : Level 3 Spatial pars
    sigma2_S_tmp <-
      sample_sigma2_S(
        sigma2_S[i - 1],
        S_k,
        beta[i - 1],
        R_S_inv,
        littlen,
        k,
        sigma2_S_prior_a = PRIOR_SIGMA2_S_A,
        sigma2_S_prior_b = PRIOR_SIGMA2_S_B
      )
    sigma2_S[i] <- sigma2_S_tmp$sigma2_S
    sigma2_S_acc_rate[i] <- sigma2_S_tmp$sigma2_S_acc
    
    rho_S_tmp <-
      sample_rho_S(
        rho_S[i - 1],
        sigma2_S[i],
        S_k,
        all_coords,
        rho_prior_a = PRIOR_RHO_S_A,
        rho_prior_b = PRIOR_RHO_S_B
      )
    rho_S[i] <- rho_S_tmp$rho_S
    rho_S_acc_rate[i] <- rho_S_tmp$rho_S_acc
    
    # Step 6:  Preferential par beta
    beta_tmp <-
      sample_beta(
        beta[i - 1],
        sigma2_S[i],
        S_k,
        littlen,
        k[i],
        beta_prior_mean = PRIOR_BETA_MEAN,
        beta_prior_var = PRIOR_BETA_VAR
      )
    beta[i] <- beta_tmp$beta
    beta_acc_rate[i] <- beta_tmp$beta_acc
  }
  
  set_priors = list(
    lambda = paste0("Gamma(", PRIOR_LAMBDA_A, ",", PRIOR_LAMBDA_B, ")"),
    eta = paste0("n(", PRIOR_ETA_MEAN, ",", PRIOR_ETA_VAR, ")"),
    nu = paste0("n(", PRIOR_NU_MEAN, ",", PRIOR_NU_VAR, ")"),
    xi = paste0("logn(", PRIOR_XI_MEAN, ",", PRIOR_XI_VAR, ")"),
    sigma2 = paste0("IG(", PRIOR_SIGMA2_S_A, ",", PRIOR_SIGMA2_S_B, ")"),
    rho = paste0("Gamma(", PRIOR_RHO_S_A, ",", PRIOR_RHO_S_B, ")"),
    beta = paste0("n(", PRIOR_BETA_MEAN, ",", PRIOR_BETA_VAR, ")")
  )
  
  ## save to a list.
  tmp_results <-   list(
    lambda_star = lambda_star,
    k = k,
    eta = eta,
    nu = nu,
    xi = xi,
    beta = beta,
    S_n = S_n,
    sigma2_S = sigma2_S,
    rho_S = rho_S,
    sim_data = sim_data,
    other_summaries = list(
      eta_acc_rate = eta_acc_rate,
      nu_acc_rate = nu_acc_rate,
      xi_acc_rate = xi_acc_rate,
      beta_acc_rate = beta_acc_rate,
      sigma2_S_acc_rate = sigma2_S_acc_rate,
      rho_S_acc_rate = rho_S_acc_rate,
      priors = set_priors
    )
  )
  # Save Results to disk.
  saveRDS(tmp_results, out_file_loc)
  return(tmp_results)
}

cur_info_print <- function(i) {
  print(paste0("sigma2_S: ", sigma2_S[i - 1]))
  print(paste0("eta: ", eta[i - 1]))
  print(paste0("nu: ", nu[i - 1]))
  print(paste0("xi: ", xi[i - 1]))
  print(paste0("lambda_star ", lambda_star[i - 1]))
  print(paste0("rho_S ", rho_S[i - 1]))
  print(paste0("beta ", beta[i - 1]))
  print(paste0("Sn_1: ", S_n[1, i - 1]))
  print(paste0("Sn_2: ", S_n[2, i - 1]))
}
