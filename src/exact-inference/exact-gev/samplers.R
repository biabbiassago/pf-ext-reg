source(here::here("src/bGEV/bGEVcode.R"))

#### Poisson Process Part ####
sample_lambda_star <-
  function(k_cur,
           lambda_prior_a,
           lambda_prior_b,
           area_B = 1) {
    return(rgamma(1, lambda_prior_a + k_cur, lambda_prior_b + area_B))
  }

sample_all_coords <-
  function(lambda_star,
           beta,
           sigma2_S,
           rho_S,
           S_k,
           obs_coords,
           all_coords,
           area_B = 1) {
    # Simulate a homgeneous pp with intensity lambda star
    k <- dim(all_coords)[1]
    k_star <- rpois(1, lambda_star * area_B)
    pp <- spatstat.random::runifpoint(k_star)
    # Generate a set of S_x_tilde of this PP conditional on the previous S_k
    R_full <-
      fields::Matern(fields::rdist(rbind(coords(pp), all_coords)),
                     range = rho_S,
                     smoothness = 1)
    cov_full <- sigma2_S * R_full
    cov_12 <-  cov_full[(1:k_star), (k_star + 1):(k + k_star)]
    cov_21 <- t(cov_12)
    cov_11 <- cov_full[1:k_star, 1:k_star]
    cov_22 <-
      cov_full[(k_star + 1):(k_star + k), (k_star + 1):(k_star + k)]
    inv_cov_22 <- chol2inv(chol(cov_22))
    S_x_tilde_mean <- cov_12 %*% inv_cov_22 %*% (S_k)
    S_x_tilde_cov <- cov_11 - cov_12 %*% inv_cov_22 %*% cov_21
    S_x_tilde <-
      mvrnorm(1, mu = S_x_tilde_mean, Sigma = S_x_tilde_cov)
    
    # calculate discard probabilities
    probs <- pnorm(-beta * S_x_tilde / sqrt(sigma2_S))
    # select "discarded" locations according to these probabilities
    disc_coords <-
      coords(pp)[(sapply(probs, function(x) {
        rbinom(1, 1, x)
      }) == T), ]
    all_coords <- rbind(obs_coords, disc_coords)
    return(all_coords)
  }

#### Spatial Random Effects ####

sample_Sk_ess <- function(S_k_cur,
                          beta,
                          sigma2_S,
                          rho_S,
                          eta,
                          nu,
                          xi,
                          y,
                          littlen,
                          lambda_star,
                          obs_coords,
                          all_coords_prev,
                          all_coords_new,
                          area_B = 1) {
  k_new <- dim(all_coords_new)[1]
  R <-
    fields::Matern(fields::rdist(all_coords_new),
                   range = rho_S,
                   smoothness = 1)
  SIGMA <- sigma2_S * R
  R_inv <- chol2inv(chol(R))
  ellipse_ess <- as.vector(rmvnorm(1, Sigma = SIGMA))
  
  u <- runif(1)
  
  # Conditional on the S_k vector at the previous, we get the values at the new locations.
  S_k_cur_atnew <-
    get_Sk_atnew(S_k_cur,
                 rho_S,
                 sigma2_S,
                 all_coords_prev,
                 all_coords_new,
                 littlen)
  logy <-
    lik_Sk(S_k_cur_atnew,
           R_inv,
           beta,
           sigma2_S,
           eta,
           nu,
           xi,
           y,
           littlen,
           lambda_star,
           k_new) + log(u)
  
  # define the bracket
  theta <- runif(1, 0, 2 * pi)
  theta_min <- theta - 2 * pi
  theta_max <- theta
  
  ess_sk_targ <- function(theta) {
     # this auxiliary function generates a "proposal" S_k and the associated likelihood value to compare to the previous one.
    S_k_new <- S_k_cur_atnew * cos(theta) + ellipse_ess * sin(theta)
    # p(data|S_k_new)) (logscale)
    tmp_targ <-
      lik_Sk(S_k_new,
             R_inv,
             beta,
             sigma2_S,
             eta,
             nu,
             xi,
             y ,
             littlen,
             lambda_star,
             k_new)
    return(list("s" = S_k_new, "t" = tmp_targ))
  }
  ess_step <- ess_sk_targ(theta)
  tmp_targ <- ess_step$t
  
  while (tmp_targ < logy) {
    if (theta < 0) {
      theta_min <- theta
    }
    else{
      theta_max <- theta
    }
    theta <- runif(1, theta_min, theta_max)
    ess_step <- ess_sk_targ(theta)
    tmp_targ <- ess_step$t
  }
   # only once tmp_targ > log(y) , exit the while loop and return the new value of S_k
  return(ess_step$s)
}

get_Sk_atnew <-
  function(S_k_cur,
           rho_S,
           sigma2_S,
           all_coords_prev,
           all_coords_new,
           littlen) {
    k_prev <- dim(all_coords_prev)[1] # k[i-1]
    k_new <- dim(all_coords_new)[1] # k[i]
    
    # discarded coordinates at iteration i
    disc_coords <- all_coords_new[(littlen + 1):k_new,]
    
    # Here i generate S_k_cur AT THE NEW LOCATIONS conditional on the previous ones
    big_coords <- rbind(all_coords_prev, disc_coords)
    Sigma_full <-
      sigma2_S * fields::Matern(fields::rdist(big_coords),
                                range = rho_S,
                                smoothness = 1)
    k_big <- dim(big_coords)[1]
    Sigma_11 <- Sigma_full[1:k_prev, 1:k_prev]
    Sigma_11_inv <- chol2inv(chol(Sigma_11))
    
    Sigma_12 <- Sigma_full[(1:k_prev), (k_prev + 1):k_big]
    Sigma_21 <- Sigma_full[(k_prev + 1):k_big, 1:k_prev]
    Sigma_22 <- Sigma_full[(k_prev + 1):k_big, (k_prev + 1):k_big]
    
    cond_mean <- Sigma_21 %*% Sigma_11_inv %*% S_k_cur
    cond_Sigma <- Sigma_22 - Sigma_21 %*% Sigma_11_inv %*% Sigma_12
    #
    disc_Sk_for_current <- MASS::mvrnorm(1, cond_mean, cond_Sigma)
    
    # recall: the first littlen elements remain the same (same, observed locations).
    S_k_cur_atnew <- c(S_k_cur[1:littlen], disc_Sk_for_current)
    return(S_k_cur_atnew)
  }

sample_Sk_mh <- function(S_k_cur,
                         R_inv_cur,
                         beta,
                         sigma2_S,
                         rho_S,
                         eta,
                         nu,
                         xi,
                         y,
                         littlen,
                         obs_coords,
                         all_coords_prev,
                         all_coords_cur) {
  # TODO
  
  
}

lik_Sk <-
  function(S_k_tmp,
           R_inv_tmp,
           beta,
           sigma2_S,
           eta,
           nu ,
           xi ,
           y,
           littlen,
           lambda_star,
           k) {
    I_k <- c(rep(1, littlen), rep(-1, k - littlen))
    t1 <- sum(pnorm((beta / sqrt(sigma2_S)) * I_k * S_k_tmp,
                    log = T))
    tmp <-
      cbind(
        "y" = y,
        "q" = eta + S_k_tmp[1:littlen],
        "sb" = rep(exp(nu), length(y)),
        "xi" = rep(xi, length(y))
      )
    
    t2 <- sum(apply(tmp, 1, function(x) {
      dbgev2(x["y"], x["q"], x["sb"], x["xi"], log = T)
    }))
    
    #t3 <- (-lambda_star) + k*(lambda_star) - sum(log(seq(1:k)))
    return(t1 + t2)
  }

target_Sk <-
  function(S_k_tmp,
           R_inv_tmp,
           beta,
           sigma2_S,
           eta,
           nu,
           xi,
           y,
           littlen,
           k) {
    I_k <- c(rep(1, littlen), rep(-1, k - littlen))
    t1 <- sum(pnorm((beta / sqrt(sigma2_S)) * I_k * S_k_tmp,
                    log = T))
    tmp <-
      cbind(
        "y" = y,
        "q" = eta + S_k_tmp[1:littlen],
        "sb" = rep(exp(nu), length(y)),
        "xi" = rep(xi, length(y))
      )
    
    t2 <- sum(apply(tmp, 1, function(x) {
      dbgev2(x["y"], x["q"], x["sb"], x["xi"], log = T)
    }))
    
    
    t3 <-
      (-1 / (2 * sigma2_S)) * t(S_k_tmp) %*% R_inv_tmp %*% S_k_tmp
    return(t1 + t2 + t3)
  }


#### GEV Part ####

sample_eta <-
  function(eta_cur,
           nu,
           xi,
           S_n,
           y,
           prior_eta_mean = 0,
           prior_eta_var = 100) {
    eta_acc <- 0
    eta_prop <- rnorm(1, eta_cur, 0.7)
    log_mh_ratio <-
      target_eta(eta_prop, nu, xi, S_n, y, prior_eta_mean, prior_eta_var) -  target_eta(eta_cur, nu, xi, S_n, y, prior_eta_mean, prior_eta_var)
    
    if (runif(1) < exp(log_mh_ratio)) {
      eta_acc <- 1
      eta_cur <- eta_prop
    }
    return(list(eta = eta_cur, eta_acc = eta_acc))
  }

target_eta <-
  function(eta,
           nu,
           xi,
           S_n,
           y,
           prior_eta_mean,
           prior_eta_var) {
    tmp <-
      cbind(
        "y" = y,
        "q" = eta + S_n,
        "sb" = rep(exp(nu), length(y)),
        "xi" = rep(xi, length(y))
      )
    
    t1 <- sum(apply(tmp, 1, function(x) {
      dbgev2(x["y"], x["q"], x["sb"], x["xi"], log = T)
    }))
    
    t2 <- (-1 / 2) * (eta - prior_eta_mean) ^ 2 / prior_eta_var
    return(t1 + t2)
  }
sample_nu <-
  function(nu_cur,
           eta,
           xi,
           S_n,
           y,
           prior_nu_mean = 0,
           prior_nu_var = 100) {
    nu_acc <- 0
    nu_prop <- rnorm(1, nu_cur, 0.4)
    log_mh_ratio <-
      target_nu(nu_prop, eta, xi, S_n, y, prior_nu_mean, prior_nu_var) - target_nu(nu_cur, eta, xi, S_n, y, prior_nu_mean, prior_nu_var)
    
    if (runif(1) < exp(log_mh_ratio)) {
      nu_acc <- 1
      nu_cur <-  nu_prop
    }
    return(list(nu = nu_cur, nu_acc = nu_acc))
  }
target_nu <-
  function(nu,
           eta,
           xi,
           S_n,
           y,
           prior_nu_mean,
           prior_nu_var) {
    tmp <-
      cbind(
        "y" = y,
        "q" = eta + S_n,
        "sb" = rep(exp(nu), length(y)),
        "xi" = rep(xi, length(y))
      )
    
    t1 <- sum(apply(tmp, 1, function(x) {
      dbgev2(x["y"], x["q"], x["sb"], x["xi"], log = T)
    }))
    
    t2 <- (-1 / 2) * (nu - prior_nu_mean) ^ 2 / prior_nu_var
    return(t1 + t2)
  }


sample_xi <-
  function(xi_cur,
           eta,
           nu,
           S_n,
           y,
           prior_xi_mean = -1,
           prior_xi_var = 1) {
    xi_acc <- 0
    xi_prop <- rlnorm(1, log(xi_cur), 1)
    # target + adjust for asymmetrical proposal
    num <-
      target_xi(xi_prop, eta, nu, S_n, y, prior_xi_mean, prior_xi_var)
    den <-
      target_xi(xi_cur, eta, nu, S_n, y, prior_xi_mean, prior_xi_var)
    log_mh_ratio <-
      num - den - log(xi_cur) + log(xi_prop)
    if (runif(1) < exp(log_mh_ratio)) {
      xi_acc <- 1
      xi_cur <- xi_prop
    }
    return(list(xi = xi_cur, xi_acc = xi_acc))
    
  }
target_xi <-
  function(xi_tmp,
           eta,
           nu,
           S_n,
           y,
           prior_xi_mean,
           prior_xi_var) {
    # Remember xi must always be positive.
    tmp <-
      cbind(
        "y" = y,
        "q" = eta + S_n,
        "sb" = rep(exp(nu), length(y)),
        "xi" = rep(xi_tmp, length(y))
      )
    
    t1 <- sum(apply(tmp, 1, function(x) {
      dbgev2(x["y"], x["q"], x["sb"], x["xi"], log = T)
    }))
    t2 <- dlnorm(xi_tmp, prior_xi_mean, prior_xi_var)
    return(t1 + t2)
  }



#### Spatial Hyper-parameters Part ####

sample_sigma2_S <-
  function(sigma2_S_cur,
           S_k,
           beta,
           R_S_inv,
           littlen,
           k,
           sigma2_S_prior_a = 1,
           sigma2_S_prior_b = 5) {
    ## Let the prior for Sigma2 as IG (a,b)
    sigma2_S_acc <- 0
    sigma2_S_prop <- stats::rlnorm(1, log(sigma2_S_cur), 0.2)
    log_mh_ratio <-
      sigma2_mh_logratio(
        sigma2_S_cur,
        sigma2_S_prop,
        S_k,
        beta,
        R_S_inv,
        littlen,
        k,
        sigma2_S_prior_a,
        sigma2_S_prior_b
      )
    
    comp <- runif(1) < exp(log_mh_ratio)
    if (comp) {
      sigma2_S_acc <- 1
      sigma2_S_cur <- sigma2_S_prop
    }
    return(list(sigma2_S = sigma2_S_cur, sigma2_S_acc = sigma2_S_acc))
  }

sigma2_mh_logratio <-
  function(sigma2_S_cur,
           sigma2_S_prop,
           S_k,
           beta,
           R_S_inv,
           littlen,
           k,
           sigma2_S_prior_a,
           sigma2_S_prior_b) {
    k <- length(S_k)
    
    # vector of 1 and -1 that differentiates 1:littlen, (littlen+1):k
    I_k <- c(rep(1, littlen), rep(-1, k - littlen))
    
    term1 <- sum(pnorm((beta / sqrt(sigma2_S_prop)) * I_k * S_k,
                       log = T)) - sum(pnorm((beta / sqrt(sigma2_S_cur)) * I_k * S_k,
                                             log = T))
    
    term2 <-
      (-length(S_k) / 2 - sigma2_S_prior_a) * (log(sigma2_S_prop) - log(sigma2_S_cur))
    
    term3 <-
      -((1 / 2) * (t(S_k) %*% R_S_inv %*% S_k) + sigma2_S_prior_b) * (1 / sigma2_S_prop - 1 /
                                                                        sigma2_S_cur)
    
    return(term1 + term2 + term3)
  }


sample_rho_S <-
  function(rho_cur,
           sigma2_S,
           S_k,
           all_coords,
           rho_prior_a = 2000,
           rho_prior_b = 4000) {
    rho_acc <- 0
    #print(rho_cur)
    rho_prop <- stats::rlnorm(1, log(rho_cur), 0.3)
    #print(rho_prop)
    log_mh_ratio_rho <-
      rho_S_mh_logratio(rho_cur,
                        rho_prop,
                        sigma2_S,
                        S_k,
                        rho_prior_a,
                        rho_prior_b,
                        all_coords)
    rinf_num <- runif(1)
    
    if (runif(1) < exp(log_mh_ratio_rho)) {
      rho_acc <- 1
      rho_cur <- rho_prop
    }
    #print(rho_acc)
    return(list(rho_S = rho_cur, rho_S_acc = rho_acc))
    
  }

rho_S_mh_logratio <-
  function(rho_cur,
           rho_prop,
           sigma2_S,
           S_k,
           rho_prior_a,
           rho_prior_b,
           all_coords) {
    #R_prop <- matern(fields::rdist(all_coords),1/rho_prop,kappa=1)
    #R_cur <- matern(fields::rdist(all_coords),1/rho_cur,kappa=1)
    
    R_S_prop <-
      fields::Matern(fields::rdist(all_coords),
                     range = rho_prop,
                     smoothness = 1)
    R_S_cur <-
      fields::Matern(fields::rdist(all_coords),
                     range = rho_cur,
                     smoothness = 1)
    
    #R_prop <- exp_cov(1,1/rho_prop,fields::rdist(all_coords))
    #R_cur <- exp_cov(1,1/rho_prop,fields::rdist(all_coords))
    
    term1 <-
      (-1 / 2) * (det(R_S_prop, log = T) - det(R_S_cur, log = T))
    term2 <- (rho_prior_a) * (log(rho_prop) - log(rho_cur))
    term3 <-
      -(1 / (2 * sigma2_S)) * (t(S_k) %*% (chol2inv(chol(R_S_prop)) - chol2inv(chol(R_S_cur))) %*%
                                 S_k) - rho_prior_b * (rho_prop - rho_cur)
    
    return(term1 + term2 + term3)
  }




#### Preferential Parameter Part ####
sample_beta <-
  function(beta_cur,
           sigma2_S,
           S_k,
           littlen,
           k,
           beta_prior_mean = 0,
           beta_prior_var = 1) {
    beta_acc <- 0
    beta_prop <- rnorm(1, beta_cur, 1.0)
    
    log_mh_ratio <-
      target_beta_log(
        beta_prop,
        sigma = sqrt(sigma2_S),
        S_k = S_k,
        littlen = littlen,
        k = k,
        beta_prior_mean = beta_prior_mean,
        beta_prior_var = beta_prior_var
      ) -
      target_beta_log(
        beta_cur,
        sigma = sqrt(sigma2_S),
        S_k = S_k,
        littlen = littlen,
        k = k,
        beta_prior_mean = beta_prior_mean,
        beta_prior_var = beta_prior_var
      )
    
    # accept or reject step
    if (runif(1) < exp(log_mh_ratio)) {
      beta_acc <- 1
      beta_cur <- beta_prop
    }
    return(list(beta = beta_cur, beta_acc = beta_acc))
  }


target_beta_log <-
  function(beta,
           sigma,
           S_k,
           littlen,
           k,
           beta_prior_mean,
           beta_prior_var) {
    # vector of 1 and -1 that differentiates 1:littlen, (littlen+1):k
    I_k <- c(rep(1, littlen), rep(-1, k - littlen))
    return(sum(pnorm((beta / sigma) * I_k * S_k, log = T)) +
             dnorm(beta, beta_prior_mean, sqrt(beta_prior_var), log = T))
  }