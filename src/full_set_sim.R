# simulate whether for full dataset method works.
# example on how to run with 400 obs (o.05 gaps in the unit square, PHIX=0.5, SIGMA2_ETA = 1.0)
#Rscript src/full_set_sim.R 0.05 0.5 1.0


# PASS: GAP_SIZE, PHIX, SIGMA2, PROPORTION OF FULL SAMPLE


args <- commandArgs(trailingOnly = TRUE)
GAP_SIZE <- as.numeric(args[[1]])
PHIX <- as.numeric(args[[2]])
SIGMA2 <- as.numeric(args[[3]])

source(here::here("src/sim-gev/sim-timeconst.R"))
source(here::here("src/testing_trend.R"))
source(here::here("src/utils.R"))

b<-0
nsims<-50



if(args[[4]]=="all"){
  SAMPLE_SIZE <- N
} else {
  SAMPLE_SIZE <- as.numeric(args[[4]])*N
}


full_rslts <- matrix(NA,nrow=nsims,ncol=4)
for(i in 1:nsims){
    print(paste("Simulation",i,"out of",nsims))
    simdata <- make_true_data()
    true_data <- simdata$true_data
    mean_by_station <- make_mean_by_station_df(true_data)
    no_sampling_results <- get_regression_estimates(b,mean_by_station,true_data)
    full_rslts[i,] <- unlist(no_sampling_results[1:4])
}

saveRDS(
  list(
    "full_rslts"=full_rslts,
    "params"=
      c("PHIX"=PHIX,"SIGMA2X"=SIGMA2X, "PHI"=PHI,"SIGMA2"=SIGMA2,"N"=N)
  ),
  file = paste0("outputs/full_set_sim/sim-no-sampling-PHIX",PHIX,"N",N,"SIGMA2",SIGMA2,"SAMPSIZE",SAMPLE_SIZE,".rds")
)

# rslt_summ <- apply(full_rslts,2,mean)
# write.csv(
#   rslt_summ,
#   file = paste0(
#     "summary-sim-no-sampling-PHIX",
#     PHIX,
#     "N",
#     N,
#     ".csv")
# )


