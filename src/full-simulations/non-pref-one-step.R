source(here::here("src/sim-gev/sim-timeconst.R"))
source(here::here("src/bayesian-onestep/one-step-stan.R"))

args <- commandArgs(trailingOnly=TRUE)
sim_number <- args[1]
S <- as.numeric(args[2])
MONTHS <- as.numeric(args[3])

# for testing
# source(here::here("src/sim-gev/sim-timeconst.R"))
# source(here::here("src/bayesian-onestep/one-step-stan.R"))
# 
# S <- 50
# MONTHS <- 30

fd <- make_true_data(S,MONTHS)
results <- one_step_model(fd,iter = 100000)


path_name <- paste0(
  "outputs/full-simulations/non-pref-one-step/npts",
  S,
  "stations",
  MONTHS,
  "months",
  sim_number,
  today(),
  ".rds"
)
saveRDS(results,here::here(path_name))


