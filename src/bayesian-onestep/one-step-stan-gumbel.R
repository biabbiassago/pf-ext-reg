library(rstan)
options(mc.cores = parallel::detectCores())
# helper function
prepare_stan_data <- function(sim_dat){
  x_mat <- cbind(rep(1,length(sim_dat$true_mus)),unique(sim_dat$true_data$x))
  unique_coords <- sim_dat$true_data %>%
    group_by(station) %>%
    summarize(x_coords = first(x_coords), y_coords = first(y_coords))
  loc <- cbind(unique_coords$x_coords,unique_coords$y_coords)
  
  distance_mat <- as.matrix(dist(loc))
  V <- matrix(c(1,0,0,1),nrow=2,byrow=TRUE)
  
  s <- length(sim_dat$true_mus)
  months <- dim(sim_dat$true_data)[1]/length(sim_dat$true_mus)
  
  stan_data <- list(
    s = s,
    months = months,
    N = s*months,
    z = sim_dat$true_data$value,
    x = x_mat,
    DMat = distance_mat,
    noise = rnorm(s,0,0.5),
    V = V
  )
  return(stan_data)
}
one_step_model <- function(fd,iter=10000,warmup=1000,chains=1){
    fit1 <- "
    functions {
      // exponential covariance function def
        matrix exp_cov(matrix DMat, real sigma_sq, real scale, vector noise) {
          int S = dims(DMat)[1];
          matrix[S, S] K;
          for (i in 1:(S-1)) {
            K[i, i] = sigma_sq + noise[i];
            for (j in (i + 1):S) {
              K[i, j] = sigma_sq * exp(- DMat[i,j] / scale);
              K[j, i] = K[i, j];
            }
          }
          K[S, S] =sigma_sq + noise[S] ;
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
      vector[s] noise;
    
      matrix[2,2] V;
      
    }
    parameters {
      vector[s] mu;
      vector<lower=0>[s] sigma2;
      real<lower=0.1> beta0;
      real<lower=0.1> beta1;
      vector[2] alpha;
    }
    model {
      int pos = 1;
      int pos_scale = 1;
      matrix[s,s] SIGMA;
      row_vector[s] m;
      vector[N] locs;
      vector[N] scales;
      SIGMA = exp_cov(DMat, beta0, beta1,noise);
      
    matrix[s,s] L_S = cholesky_decompose(SIGMA);
      
    // priors
    alpha ~ multi_normal(rep_vector(0,2), V);
    m = to_row_vector(alpha)*x';
    
    target += lognormal_lpdf(beta0 | 1,0.5);
    target += lognormal_lpdf(beta1 | -0.5,0.5);
      
    for(i in 1:s){
        target += inv_gamma_lpdf(sigma2[i]| 3,1);
    }  
    
    // data model
    target += multi_normal_cholesky_lpdf(mu |m,L_S);
    
    //gumbel likelihood
    pos = 1;
    for (i in 1:s) {
        segment(z, pos, months) ~ gumbel(mu[i], sigma2[i]);
        pos = pos + s;
      }
  }
  "
  stan_data <- prepare_stan_data(fd)
  stan_fit1 <- stan(
    model_code = fit1,
    data = stan_data,
    chains = chains,
    warmup = warmup,
    iter = iter,
    control = list(adapt_delta = 0.999, max_treedepth=13)
  )
  
  
  
  # get parameters of interest
  mu_est <- summary(stan_fit1,pars="mu")$summary[,"mean"]
  mu_025 <- summary(stan_fit1,pars="mu")$summary[,"2.5%"]
  mu_975 <- summary(stan_fit1,pars="mu")$summary[,"97.5%"]
  
  alpha0_est <- summary(stan_fit1,pars="alpha")$summary[1,"mean"]
  alpha0_low <- summary(stan_fit1,pars="alpha")$summary[1,"2.5%"]
  alpha0_hi <-  summary(stan_fit1,pars="alpha")$summary[1,"97.5%"]
  
  
  alpha1_est <- summary(stan_fit1,pars="alpha")$summary[2,"mean"]
  alpha1_low <- summary(stan_fit1,pars="alpha")$summary[2,"2.5%"]
  alpha1_hi <-  summary(stan_fit1,pars="alpha")$summary[2,"97.5%"]
  
  
  beta0_est <- summary(stan_fit1,pars="beta0")$summary[,"mean"]
  beta0_low <- summary(stan_fit1,pars="beta0")$summary[,"2.5%"]
  beta0_hi <- summary(stan_fit1,pars="beta0")$summary[,"97.5%"]
  
  beta1_est <- summary(stan_fit1,pars="beta1")$summary[,"mean"]
  beta1_low <- summary(stan_fit1,pars="beta1")$summary[,"2.5%"]
  beta1_hi <- summary(stan_fit1,pars="beta1")$summary[,"97.5%"]
  
  sigma2_est <- summary(stan_fit1,pars="sigma2")$summary[,"mean"]
  sigma2_low <- summary(stan_fit1,pars="sigma2")$summary[,"2.5%"]
  sigma2_hi <- summary(stan_fit1,pars="sigma2")$summary[,"97.5%"]
  
  
  
  
  rslts <- as.list(numeric(7*4))
  dim(rslts) <- c(7,4)
  
  # mus
  rslts[[1,1]] <- mu_est
  rslts[[1,2]] <- mu_025
  rslts[[1,3]] <- mu_975
  rslts[[1,4]] <- fd$true_mus
  
  # alpha0
  rslts[[2,1]] <- alpha0_est
  rslts[[2,2]] <- alpha0_low
  rslts[[2,3]] <- alpha0_hi
  rslts[[2,4]] <- fd$data_parms$ALPHA0
  
  #alpha1
  rslts[[3,1]] <- alpha1_est
  rslts[[3,2]] <- alpha1_low
  rslts[[3,3]] <- alpha1_hi
  rslts[[3,4]] <- fd$data_parms$ALPHA1
  
  #beta0
  rslts[[4,1]] <- beta0_est
  rslts[[4,2]] <- beta0_low
  rslts[[4,3]] <- beta0_hi
  rslts[[4,4]] <- fd$data_parms$SIGMA2ETA
  
  #beta1
  rslts[[5,1]] <- beta1_est
  rslts[[5,2]] <- beta1_low
  rslts[[5,3]] <- beta1_hi
  rslts[[5,4]] <- fd$data_parms$PHIETA
  
  #sigma2 of the gev 
  rslts[[6,1]] <- sigma2_est
  rslts[[6,2]] <- sigma2_low
  rslts[[6,3]] <- sigma2_hi
  rslts[[6,4]] <- fd$data_parms$SCALE
  
  #xi
  rslts[[7,1]] <- NA
  rslts[[7,2]] <- NA
  rslts[[7,3]] <- NA
  rslts[[7,4]] <- fd$data_parms$SHAPE
  
  colnames(rslts) <- c("est","95ci_low","95ci_hi","truth")
  rownames(rslts) <- c("mu","alpha0","alpha1","beta0","beta1","sigma2","xi")
  return(rslts)
}


