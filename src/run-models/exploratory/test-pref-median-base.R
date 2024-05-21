## run same but with no pref adj
# fit uncorrected model.
library(rstan)
source(here::here("src/run-models/run-funcs.R"))

dat_sample <- readRDS(here::here("data/test-median-pref/prefonmedian-0.898.rds"))
stan_data <- prepare_stan_data_no_pref(dat_sample)
random_data_id <- 0.898

CODE_PATH <- here::here("src/stan-files/fit-repgev.stan")
init_vals <- crude_init_values(dat_sample,dat_sample$data_parms$s,chains = 1)
#rm(dat_sample)
stan_fit_median_no_pref <- stan(
  file = CODE_PATH,
  #model_code = stan_model,
  data = stan_data,
  chains = 1,
  iter = 5000,
  init = init_vals,
  warmup = 2000,
  pars = c("SIGMAETA","SIGMANU","SIGMAPHI","L_SE","L_SN","L_SPHI"),
  include=FALSE,
  control= list(max_treedepth = 12)
)

saveRDS(
  stan_fit_median_no_pref,here::here(paste0("outputs/test-pref-on-median-0524/modelbaseiters5k-",random_data_id,".rds")))
