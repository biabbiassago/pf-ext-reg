source(here::here("src/exact-inference/exact-gev/sim-data.R"))
source(here::here("src/exact-inference/exact-gev/mcmc.R"))

args <- commandArgs(trailingOnly=TRUE)
qq <- ifelse(is.na(args[1]),"any",args[1])
FILENAME <- paste0("largecheck-gev-exact-",qq)
start <- Sys.time()
OUT_FILE_LOC <- here::here(paste0("outputs/mcmc-exact/",FILENAME,".rds"))

# Simulate data
sim_data <- make_sim_data_max(true_lambda_star = 300)
obs_coords <- sim_data$obs_coords
y <- sim_data$y

# Run model
mod <- gev_exact_mcmc(y,obs_coords,NSIMS,OUT_FILE_LOC)
end <- Sys.time()
print(end-start)