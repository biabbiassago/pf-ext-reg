library(tidyverse)
library(patchwork)
source("src/utils.R")
set.seed(3523)

##### Auxiliary function ####
generate_sampling_probs <- function(b,dat){
  n <- length(dat$x)

  ### check
  #pi <- sapply(BETA1*dat$x + eta, function(x) max(0,x))
  eps <- rnorm(n,0,1)
  pi <- abs(b*dat$mean_value+eps)
  pi <- pi/sum(pi)
  return(list("pi"=pi,"b"=b))
}
sample_by_stations <- function(dat,prob){
  ### Taking the Sample ####
  s_i <- sample(
      1:N,
      size=SAMPLE_SIZE,
      replace=FALSE,
      prob=prob$pi
  )
  
  # vector of sampled stations
  #r_s <- ifelse(1:N %in% s_i, 1, 0)
  
  sample_df <- dat %>% filter(station %in% s_i)
  return(list("sample_df"=sample_df,"stations"=s_i,"b"=prob$b))
}
##### Main Function ####

prefsamp_station <- function(b,dat_mean, dat_full){
  
  pi <- generate_sampling_probs(b,dat_mean)
  x <- sample_by_stations(dat_full,pi)
  return(x)
}



#### Plot function for checking ####
make_samp_plot <- function(cur_samp){
  
  sampled_stations <- cur_samp$stations
  
  df_plot <- mean_by_station  %>% mutate(
    sampled_measurements = case_when(
      station %in% sampled_stations ~ mean_value,
      TRUE ~ NA,
    )
  )
  
  p1 <- ggplot(df_plot, aes(x = x_coords, y = y_coords, fill = sampled_measurements)) + 
    geom_raster() +
    xlab("x") + 
    ylab("y") + 
    scale_fill_viridis_c(option="viridis",direction=-1) +
    theme_bw()
  
  p2 <- ggplot(mean_by_station, aes(x=x_coords, y=y_coords, fill = mean_value)) +   geom_raster() +
    xlab("x") + 
    ylab("y") + 
    scale_fill_viridis_c(option="viridis",direction=-1) +
    theme_bw()
  
  p1 + p2 + plot_annotation(title=sprintf("Sample with b=%i and full data",cur_samp$b))
  
}



####### RUN #########
#if you want to run just for one example use this below ~ generates samples

##### get data (mean by station) ####
# true_data <- read.csv("data/sim1.csv")
# mean_by_station <- true_data %>%
#   group_by(station) %>%
#   summarize(
#     mean_value = mean(value),
#     x_coords = first(x_coords),
#     y_coords = first(y_coords),
#     x = first(x)
#   )

# xb3 <- prefsamp_station(3,mean_by_station, true_data)
# 
# sample_df <- xb3$sample_df
# sampled_stations <- xb3$stations
# make_samp_plot(xb3)
# 
# 
# xb1 <- prefsamp_station(1,mean_by_station, true_data)
# make_samp_plot(xb1)
# 
# 
# xb0 <- prefsamp_station(0,mean_by_station, true_data)
# make_samp_plot(xb0)
