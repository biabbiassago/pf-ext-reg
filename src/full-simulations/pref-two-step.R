source(here::here("src/sim-gev/sim-timeconst-pref.R"))
source(here::here("src/two-step-model/ml_twostep_model.R"))

args <- commandArgs(trailingOnly=TRUE)
sim_number <- args[1]
S <- as.numeric(args[2])
MONTHS <- as.numeric(args[3])

fd <- make_true_data_pref(S,MONTHS)
results <- two_step_model(fd)

path_name <- paste0(
  "outputs/full-simulations/pref-two-step/npts",
  S,
  "stations",
  MONTHS,
  "months",
  sim_number,
  ".rds"
)
saveRDS(results,here::here(path_name))