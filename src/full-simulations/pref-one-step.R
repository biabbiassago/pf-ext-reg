args <- commandArgs(trailingOnly=TRUE)
sim_number <- args[1]
S <- as.numeric(args[2])
MONTHS <- as.numeric(args[3])

source(here::here("src/sim-gev/sim-timeconst-pref.R"))
source(here::here("src/bayesian-onestep/one-step-stan-gumbel.R"))
# S <- 200
# MONTHS <- 30
fd <- make_true_data_pref(S,MONTHS)
results <- one_step_model(fd)

path_name <- paste0(
  "outputs/full-simulations/pref-one-step/npts",
  S,
  "stations",
  MONTHS,
  "months",
  sim_number,
  ".rds"
)
saveRDS(results,here::here(path_name))