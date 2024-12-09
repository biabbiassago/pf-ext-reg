source(here::here("src/exact-inference/exact-gev/sim-data.R"))

start <- Sys.time()
OUT_FILE_LOC <- here::here("outputs/mcmc-exact/gev-rewritetest.rds")

# Simulate data
sim_data <- make_sim_data_max(true_lambda_star = 250)
obs_coords <- sim_data$obs_coords
y <- sim_data$y

# Run model
mod <- gev_exact_mcmc(y,obs_coords,NSIMS,OUT_FILE_LOC)
end <- Sys.time()
print(end-start)