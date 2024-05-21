source(here::here("src/utils.R"))
return_levels <- function(mu, sigma, xi, year) {
  return(mu - (sigma / xi) * (1 - (-log(1 - 1 / year)) ^ (-xi)))
}

true_return_levels <- function(fd, nyear) {
  mus <- fd$true_mus
  sigmas <- fd$true_sigmas
  xi <- fd$true_xi
  return(return_levels(mus, sigmas, xi, nyear))
}

make_traceplots_main_pars <- function(modfit, base = FALSE) {
  if (base == FALSE) {
    pt <- traceplot(modfit, c("alpha0", "rho0", "xi", "lambda","a")) 
  }
  if (base == TRUE) {
    pt <- traceplot(modfit, c("alpha0", "rho0", "xi"))
  }
  return(pt)
}


make_traceplots_spatial_pars <- function(modfit, base = FALSE) {
  if (base == FALSE) {
    pt <- traceplot(modfit, c("beta0","beta1","gamma1","spref0","spref1","eta[1]","nu[1]","phi_proc[1]","phi_proc[102]"))
  }
  if (base == TRUE) {
    pt <- traceplot(modfit,  c("beta0","beta1","gamma1","eta[1]","nu[1]"))
  }
  return(pt)
}


make_traceplots_trans_pars <- function(modfit, base = FALSE) {
  pt <- traceplot(modfit, c("mu[1]","mu[10]","q[1]","q[10]","sigma[1]","sigma[10]","iqr[1]","iqr[10]","ret_10[1]","ret_10[10]"))
  return(pt)
}



return_levels_rmse <- function(modfit, true_data, nyear = 10) {
  truth <- true_return_levels(true_data, nyear)
  est <- summary(modfit, "ret_10")$summary[, "mean"]
  return(sqrt(mean((truth - est) ^ 2)))
}

medians_rmse <- function(modfit, true_data) {
  est <- summary(modfit, "q")$summary[, "mean"]
  return(sqrt(mean((
    true_data$true_medians - est
  ) ^ 2)))
}
summary_table <-
  function(modfit,
           true_data,
           base = FALSE,
           print = TRUE,
           return = FALSE) {
    pt1 <- summary(modfit, c("alpha0", "rho0", "xi"))$summary[, c("mean", "2.5%", "50%", "97.5%", "n_eff")]
    
    pt1 <- cbind(
      pt1,
      "truth" = c(
        true_data$data_parms$ALPHA0,
        true_data$data_parms$RHO0,
        true_data$true_xi
      )
    )
    
    pt2 <-
      summary(modfit, c("beta0", "beta1", "gamma1"))$summary[, c("mean", "2.5%", "50%", "97.5%", "n_eff")]
    pt2 <-
      cbind(
        pt2,
        "truth" = c(
          true_data$data_parms$SIGMA2ETA,
          true_data$data_parms$PHIETA,
          true_data$data_parms$PHINU
        )
      )
    
    if (base == FALSE) {
      pt3 <-
        summary(modfit, c("a", "spref0", "spref1", "lambda"))$summary[, c("mean", "2.5%", "50%", "97.5%", "n_eff")]
      pt3 <- cbind(pt3, "truth" = c(1, PREF_VAR, PREF_SCALE, true_data$data_parms$s))
    }
    
    if (print == TRUE) {
      print(pt1)
      cat("-------------------------------------------------------------------")
      cat("\n")
      print(pt2)
      if (base == FALSE) {
        cat("-------------------------------------------------------------------")
        cat("\n")
        print(pt3)
      }
    }
    if (return == TRUE) {
      return(ifelse(base == TRUE, rbind(pt1, pt2), rbind(pt1, pt2, pt3)))
    }
  }


## plotting the true data
