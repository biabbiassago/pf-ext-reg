"fit EVR by station"


fit_evr_by_station <- function(i,dat){
  stat_dat <- dat %>% filter(station==i)
  fit<- fevd(value~1,data=stat_dat)
  # location: 8.166, scale = 1.05, shape = -0.117
  return(fit$results$par)
}


sapply(
  unique(sample_df$station),
  function(x) fit_evr_by_station(x,sample_df)
)
