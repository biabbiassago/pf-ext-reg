library(rstan)
source(here::here("src/simulate-data/make_data_pp.R"))
source(here::here("src/run-models/run-funcs.R"))

set.seed(4651)
DIMXY <- 10
print(DIMXY)
random_data_id <- round(runif(1,0,1),3)
dat_sample <- make_true_data_pp(100,20,0.1,a=1,pref_on="median",dimxy=DIMXY)
print("data created")
saveRDS(dat_sample,here::here(paste0("data/test-median-pref/prefonmedian-",random_data_id,".rds")))
stan_data <- make_stan_data(dat_sample)
init_vals <- crude_init_values(dat_sample,dat_sample$data_parms$s,chains = 1)
#rm(dat_sample)
CODE_PATH <- here::here("src/stan-files/fit-repgev-pref.stan")

begin<-Sys.time()

stan_fit_median <- stan(
  file = CODE_PATH,
  #model_code = stan_model,
  data = stan_data,
  chains = 1,
  iter = 8000,
  #thin=10,
  init = init_vals,
  warmup = 2000,
  pars = c("phi_grid","SIGMAETA","SIGMANU","SIGMAPHI","L_SE","L_SN","L_SPHI"),
  include=FALSE
)

saveRDS(stan_fit_median,here::here(paste0("outputs/test-pref-on-median-0524/model2grid",DIMXY,"iters8k-",random_data_id,".rds")))

end <- Sys.time()
timetake <- end-begin
print(timetake)

