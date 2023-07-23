phis <- c(0.2,0.5,1)
sample_sizes <- c(120,200,280)

for(cur_phi in phis){
  PHIX <- cur_phi
  source(here::here("src/sim-gev/sim-timeconst.R"))
  
  for(j in sample_sizes){
    SAMPLE_SIZE <- j
    source(here::here("src/testing_trend.R"))
  }
}