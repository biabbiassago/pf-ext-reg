library(rstan)
# make sure you are using the right stan code

CODE_PATH <- here::here("src/bayesian-onestep/stan-code/fit1-full-gev.stan")

stan_model <- "
// This Stan program defines the model
functions{
  // functions for GEV loglikelihood
  real vector_comp_less(vector g,real thresh){
    // returns TRUE if any elements is ABOVE the threshold
    // sample CANNOT be from xi less0
    int N = rows(g);
    real tmp;
    real comp;
    comp=0;
    for(k in 1:N){
      tmp = g[k] < thresh; 
      comp = comp + tmp;
    }
    return N != comp;
  }
  real vector_comp_more(vector g,real thresh){
    // returns TRUE if any elements is BELOW the threshold
    // sample CANNOT be from xi more0
    int N = rows(g);
    real tmp;
    real comp; // number of elements that are above the threshold
    comp=0;
    for(k in 1:N){
      tmp = g[k] > thresh; 
      comp = comp + tmp;
    }
    return N != comp;
  }
  
  real gev_lpdf(vector x, real mu, real sigma, real xi) {
    vector[rows(x)] t;
    vector[rows(x)] lp;
    int N;
    real thresh;
    thresh = (mu-sigma)/xi;
    N = rows(x);
    real comp_less;
    comp_less = vector_comp_less(x,thresh);
    
    real comp_more;
    comp_more = vector_comp_more(x,thresh);
    
    if(comp_less && xi<0){
      //reject(\"sample cannot be generated from such xi\",xi);
      for(n in 1:N){
        lp[n] = 0;
      }
      
    }
    
    else if(comp_more && xi>0){
      //reject(\"sample cannot be generated from such xi\",xi);
      for(n in 1:N){
        lp[n] = 0;
      }
    }
    
    else{
      for(n in 1:N){
        t[n] = abs(xi) < 1e-10 ? exp((mu - x[n]) / sigma) : pow(1 + xi * ((x[n] - mu ) / sigma), -1/xi);
        lp[n] = -log(sigma) + (xi + 1) * log(t[n]) - t[n];
        
      }
      
    }
    
    return sum(lp);
  }
  
  
  // exponential covariance function def
  matrix exp_cov(matrix DMat, real sigma_sq, real scale) {
    int S = dims(DMat)[1];
    matrix[S, S] K;
    for (i in 1:(S-1)) {
      K[i, i] = sigma_sq;
      for (j in (i + 1):S) {
        K[i, j] = sigma_sq * exp(- DMat[i,j] / scale);
        K[j, i] = K[i, j];
      }
    }
    K[S, S] =sigma_sq;
    return K;
  }
}
data {
  int s;
  int months;
  int N;
  vector[N] z;
  matrix[s,2] x; 
  matrix[s, s] DMat; // Distance matrix
  
  matrix[2,2] V;
  
}
parameters {
  vector[s] mu;
  vector<lower=0>[s] sigma2;
  
  real<lower=0.01> beta0;
  real<lower=0.01> beta1;
  real<lower=0.01> gamma0;
  real<lower=0.01> gamma1;
  
  row_vector[2] alpha;
  real xi;
}
model {
  
  int pos; 
  
  matrix[s,s] SIGMAETA;
  matrix[s,s] SIGMANU;
  
  row_vector[s] m;
  vector[N] locs;
  vector[N] scales;
  
  SIGMAETA = exp_cov(DMat, beta0, beta1);
  SIGMANU = exp_cov(DMat,gamma0,gamma1);
  
  matrix[s,s] L_SE = cholesky_decompose(SIGMAETA);
  matrix[s,s] L_SN = cholesky_decompose(SIGMANU);
  
  // priors
  alpha ~ multi_normal(rep_vector(0,2), V);
  // priors on variance and range of the location par
  target += lognormal_lpdf(beta0 | -2, 0.5);
  target += lognormal_lpdf(beta1 | -2,0.5);
  // priors on variance and range of the scale par
  target += lognormal_lpdf(gamma0 | -5, 0.2);
  target += lognormal_lpdf(gamma1 | -5, 0.2);
  
  xi ~ normal(0,0.25);
  
  // data model
  m = alpha*x';
  
  target += multi_normal_cholesky_lpdf(mu |m,L_SE);
  target += multi_normal_cholesky_lpdf(log(sigma2) |rep_vector(0,s),L_SN);
  
  // likelihood
  pos = 1;
  for (i in 1:s) {
    segment(z, pos, months) ~ gev(mu[i], sigma2[i],xi);
    pos = pos + months;
  }
}
"


options(mc.cores = parallel::detectCores())
# helper function
prepare_stan_data <- function(sim_dat) {
  x_mat <-
    cbind(rep(1, length(sim_dat$true_mus)), unique(sim_dat$true_data$x))
  
  unique_coords <- sim_dat$true_data %>%
    group_by(station) %>%
    summarize(x_coords = first(x_coords),
              y_coords = first(y_coords))
  
  loc <- cbind(unique_coords$x_coords, unique_coords$y_coords)
  
  distance_mat <- as.matrix(dist(loc))
  V <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  
  s <- length(sim_dat$true_mus)
  months <- dim(sim_dat$true_data)[1] / length(sim_dat$true_mus)
  
  stan_data <- list(
    s = s,
    months = months,
    N = s * months,
    z = sim_dat$true_data$value,
    x = x_mat,
    DMat = distance_mat,
    #noise = rnorm(s, 0, 0.5),
    V = V
  )
  return(stan_data)
}
make_params_table_summary <- function(stan_fit1) {
  # get parameters of interest
  mu_est <- summary(stan_fit1, pars = "mu")$summary[, "mean"]
  mu_025 <- summary(stan_fit1, pars = "mu")$summary[, "2.5%"]
  mu_975 <- summary(stan_fit1, pars = "mu")$summary[, "97.5%"]
  
  alpha0_est <- summary(stan_fit1, pars = "alpha")$summary[1, "mean"]
  alpha0_low <- summary(stan_fit1, pars = "alpha")$summary[1, "2.5%"]
  alpha0_hi <-  summary(stan_fit1, pars = "alpha")$summary[1, "97.5%"]
  
  
  alpha1_est <- summary(stan_fit1, pars = "alpha")$summary[2, "mean"]
  alpha1_low <- summary(stan_fit1, pars = "alpha")$summary[2, "2.5%"]
  alpha1_hi <-  summary(stan_fit1, pars = "alpha")$summary[2, "97.5%"]
  
  
  beta0_est <- summary(stan_fit1, pars = "beta0")$summary[, "mean"]
  beta0_low <- summary(stan_fit1, pars = "beta0")$summary[, "2.5%"]
  beta0_hi <- summary(stan_fit1, pars = "beta0")$summary[, "97.5%"]
  
  beta1_est <- summary(stan_fit1, pars = "beta1")$summary[, "mean"]
  beta1_low <- summary(stan_fit1, pars = "beta1")$summary[, "2.5%"]
  beta1_hi <- summary(stan_fit1, pars = "beta1")$summary[, "97.5%"]
  
  gamma0_est <- summary(stan_fit1, pars = "gamma0")$summary[, "mean"]
  gamma0_low <- summary(stan_fit1, pars = "gamma0")$summary[, "2.5%"]
  gamma0_hi <- summary(stan_fit1, pars = "gamma0")$summary[, "97.5%"]
  
  gamm1_est <- summary(stan_fit1, pars = "gamma1")$summary[, "mean"]
  gamma1_low <- summary(stan_fit1, pars = "gamma1")$summary[, "2.5%"]
  gamma1_hi <- summary(stan_fit1, pars = "gamma1")$summary[, "97.5%"]
  
  
  sigma2_est <- summary(stan_fit1, pars = "sigma2")$summary[, "mean"]
  sigma2_low <- summary(stan_fit1, pars = "sigma2")$summary[, "2.5%"]
  sigma2_hi <- summary(stan_fit1, pars = "sigma2")$summary[, "97.5%"]
  
  xi_est <- summary(stan_fit1, pars = "xi")$summary[, "mean"]
  xi_low <- summary(stan_fit1, pars = "xi")$summary[, "2.5%"]
  xi_hi <- summary(stan_fit1, pars = "xi")$summary[, "97.5%"]
  
  
  rslts <- as.list(numeric(9 * 4))
  dim(rslts) <- c(9, 4)
  
  # mus
  rslts[[1, 1]] <- mu_est
  rslts[[1, 2]] <- mu_025
  rslts[[1, 3]] <- mu_975
  rslts[[1, 4]] <- fd$true_mus
  
  # alpha0
  rslts[[2, 1]] <- alpha0_est
  rslts[[2, 2]] <- alpha0_low
  rslts[[2, 3]] <- alpha0_hi
  rslts[[2, 4]] <- fd$data_parms$ALPHA0
  
  #alpha1
  rslts[[3, 1]] <- alpha1_est
  rslts[[3, 2]] <- alpha1_low
  rslts[[3, 3]] <- alpha1_hi
  rslts[[3, 4]] <- fd$data_parms$ALPHA1
  
  #beta0
  rslts[[4, 1]] <- beta0_est
  rslts[[4, 2]] <- beta0_low
  rslts[[4, 3]] <- beta0_hi
  rslts[[4, 4]] <- fd$data_parms$SIGMA2ETA
  
  #beta1
  rslts[[5, 1]] <- beta1_est
  rslts[[5, 2]] <- beta1_low
  rslts[[5, 3]] <- beta1_hi
  rslts[[5, 4]] <- fd$data_parms$PHIETA
  
  
  #gamma0
  rslts[[6,1]] <- gamma0_est
  rslts[[6,2]] <- gamma0_low
  rslts[[6,3]] <- gamma0_hi
  rslts[[6,4]] <- fd$data_parms$SIGMA2NU
  
  #gamma1
  rslts[[7,1]] <- gamma1_est
  rslts[[7,2]] <- gamma1_low
  rslts[[7,3]] <- gamma1_hi
  rslts[[7,4]] <- fd$data_parms$PHINU
  
  #sigma2 of the gev
  rslts[[8, 1]] <- sigma2_est
  rslts[[8, 2]] <- sigma2_low
  rslts[[8, 3]] <- sigma2_hi
  rslts[[8, 4]] <- fd$true_sigmas
  
  #xi
  rslts[[9, 1]] <- xi_est
  rslts[[9, 2]] <- xi_low
  rslts[[9, 3]] <- xi_hi
  rslts[[9, 4]] <- fd$data_parms$SHAPE
  
  colnames(rslts) <- c("est", "95ci_low", "95ci_hi", "truth")
  rownames(rslts) <-
    c("mu", "alpha0", "alpha1", "beta0", "beta1", "gamma0","gamma1","sigma2", "xi")
  return(rslts)
}


# MAIN FUNCTION, FULL STAN MODEL.
one_step_model <- function(fd,
                           iter = 10000,
                           warmup = 1000,
                           chains = 3, test=FALSE) {
  stan_data <- prepare_stan_data(fd)
  
  stan_fit1 <- stan(
    #file = CODE_PATH, #idk why this is not working rn
    model_code = stan_model,
    data = stan_data,
    chains = chains,
    warmup = warmup,
    iter = iter,
    control = list(adapt_delta = 0.999, max_treedepth = 13),
  )
  if(test==TRUE){
    stan_fit1 <- stan(
      #file = CODE_PATH, #idk why this is not working rn
      model_code = stan_model,
      data = stan_data,
      chains = 1,
      warmup = 10,
      iter = 90,
      control = list(adapt_delta = 0.999, max_treedepth = 13),
    )
  }

  
  
  
  rslts <- make_params_table_summary(stan_fit1)
  return(rslts)
}
