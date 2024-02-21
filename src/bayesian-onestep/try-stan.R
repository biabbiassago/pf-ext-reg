library(rstan)


# Stan Code

fit1 <- "
functions {
  // exponential covariance function def
    matrix exp_cov(matrix DMat, real sigma_sq, real scale) {
      int N = dims(DMat)[1];
      matrix[N, N] K;
      for (i in 1:(N-1)) {
        K[i, i] = sigma_sq;
        for (j in (i + 1):N) {
          K[i, j] = sigma_sq * exp(- DMat[i,j] / scale );
          K[j, i] = K[i, j];
        }
      }
      K[N, N] = sigma_sq ;
      return K;
    }
}
data {
  int s;
  int months;
  int N;
  vector[s*months] z;
  matrix[s,2] x; 
  matrix[s, s] DMat; // Distance matrix

  matrix[2,2] V;
  
}
parameters {
  vector[s] mu;
  real<lower=0> sigma2;
  real<lower=0> beta0;
  real<lower=0> beta1;
  vector[2] alpha;
}
model {
  
  matrix[N,N] SIGMA;
  row_vector[s] m;
  vector[N] locs;
  
  SIGMA = exp_cov(DMat, beta0, beta1);
  
  // priors
  alpha  ~ multi_normal(rep_vector(0,2), V);
  m = to_row_vector(alpha)*x';

  target += lognormal_lpdf(beta0 | 1,0.5);
  target += lognormal_lpdf(beta1 | -0.5,0.5);
  target += inv_gamma_lpdf(sigma2| 3,1);

  // data model
  target += multi_normal_lpdf(mu| m,SIGMA);
  
  // get same location for each station
  int pos = 1;
  for(i in 1:s){
    for(t in 1:months){
      locs[pos] = mu[i];
      pos += 1;
    }
  }
  z ~ gumbel(locs, sigma2); //likelihood
}
"

prepare_stan_data <- function(sim_data){
  x_mat <- cbind(rep(1,length(sim_dat$true_mus)),unique(sim_dat$true_data$x))
  
  loc <- expand.grid(
    unique(sim_dat$true_data$x_coords),
    unique(sim_dat$true_data$y_coords)
  )
  distance_mat <- fields::rdist(loc)
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
    V = V
  )
  return(stan_data)
}

# Make stan data

### Run the following:

#No x corr
#- 4 stations
corr <- "noxcorr"
sim_dat <- readRDS("data/bayesian-onestep/simdatasetnoxcorr_4_250months.rds")
stan_data <- prepare_stan_data(sim_dat)
stan_fit1 <- stan(
  model_code = fit1,
  data = stan_data,
  chains = 3,
  warmup = 1000,
  iter = 10000,
  control = list(adapt_delta = 0.999, max_treedepth=13)
)
saveRDS(stan_fit1, paste0("outputs/bayesian-onestep/stan-models/stanmod-",corr,"-s",stan_data$s,"-months",stan_data$months,".rds"))


# - 16 stations
start <- Sys.time()
sim_dat <- readRDS("data/bayesian-onestep/simdatasetnoxcorr_16_250months.rds")
corr <- "noxcorr"
stan_data <- prepare_stan_data(sim_dat)
stan_fit2 <- stan(
  model_code = fit1,
  data = stan_data,
  chains = 1,
  warmup = 1000,
  iter = 10000,
  control = list(adapt_delta = 0.999, max_treedepth=13)
)
saveRDS(stan_fit2, paste0("outputs/bayesian-onestep/stan-models/stanmod-",corr,"-s",stan_data$s,"-months",stan_data$months,".rds"))
end <- Sys.time()
tot_time <- start-end
saveRDS(
  list("tot_time"=tot_time),
  here::here(paste0("src/bayesian-onestep/tmp/stan_timesperiter",today(),".rds"))
)

# x also spatially related:
# 4 stations
corr <- "xcorr"
sim_dat <- readRDS("data/bayesian-onestep/simdataset_4_250months.rds")
stan_data <- prepare_stan_data(sim_dat)
stan_fit1 <- stan(
  model_code = fit1,
  data = stan_data,
  chains = 3,
  warmup = 1000,
  iter = 10000,
  control = list(adapt_delta = 0.999, max_treedepth=13)
)
saveRDS(stan_fit1, paste0("outputs/bayesian-onestep/stan-models/stanmod-",corr,"-s",stan_data$s,"-months",stan_data$months,".rds"))

# - 16 stations
sim_dat <- readRDS("data/bayesian-onestep/simdataset_16_250months.rds")
corr <- "xcorr"
stan_data <- prepare_stan_data(sim_dat)
stan_fit2 <- stan(
  model_code = fit1,
  data = stan_data,
  chains = 1,
  warmup = 1000,
  iter = 10000,
  control = list(adapt_delta = 0.999, max_treedepth=13)
)
saveRDS(stan_fit2, paste0("outputs/bayesian-onestep/stan-models/stanmod-",corr,"-s",stan_data$s,"-months",stan_data$months,".rds"))

