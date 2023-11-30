set.seed(4649)
# test with the recommended number of iterations from Raftery Diagnostics.

## This specifies the number of MCMC iteration
ITER <- 470000

# 1. Simulate data. 
start <- Sys.time()
#source(here::here("src/sim-gev/sim-timeconst.R"))
#source(here::here("src/utils.R"))
source(here::here("src/bayesian-onestep/bayesian_max_model.R"))

# to_save <- make_true_data()
# df <- to_save$true_data
# 
# writeRDS(
#   to_save,
#   file=here::here("data/bayesian-onestep/simdataset.rds")
# )

to_save <- readRDS(here::here("data/bayesian-onestep/simdataset.rds"))
df <- to_save$true_data
z <- df$value
x <- df$x
loc <- true_coords
true_mus <- to_save$true_mus
s <- 100


# 2. Run MCMC

m1 <- max_model(z,x,loc,s,N.iter=ITER)
saveRDS(m1,here::here(paste0("outputs/bayesian-onestep/m1.rds")))


end <- Sys.time()
tot_time <- start-end
# save how long it took
saveRDS(
  list("iters"=ITER,"tot_time"=tot_time),
  here::here(paste0("src/bayesian-onestep/tmp/timesperiter",today(),".rds"))
)

