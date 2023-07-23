library(tidyverse)
library(patchwork)
source(here::here("src/utils.R"))
source(here::here("src/sampling_probs.R"))
set.seed(3523)

sample_by_stations <- function(dat,prob){
  ### Taking the Sample on the Ys####
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
      station %in% sampled_stations ~ 1,
      TRUE ~ 0,
    )
  )
  
  p1 <- ggplot(df_plot, aes(x = x_coords, y = y_coords, fill = factor(sampled_measurements))) + 
    geom_raster() +
    xlab("x") + 
    ylab("y") + 
    theme_bw() +
    scale_fill_discrete(guide = "none", type=c("lightgrey","tomato4")) + 
    ggtitle(paste("B=",cur_samp$b)) 
    
  return(p1)
  
}



###### RUN #########
#if you want to run just for one example use this below ~ generates samples
# 
# true_data <- read.csv(here::here("data/sim-gev-2.csv"))
# mean_by_station <- true_data %>%
#    group_by(station) %>%
#    summarize(
#     mean_value = mean(value),
#     x_coords = first(x_coords),
#     y_coords = first(y_coords),
#     x = mean(x)
# )
# 
# mean_by_station %>%
#   ggplot(aes(x=mean_value)) +
#   geom_histogram(bins=15, fill="white",color="black") +
#   theme_classic() +
#   ggtitle("Mean (monthly) value by station") +
#   xlab("Average max")
# 
# 
# xb5 <- prefsamp_station(5,mean_by_station, true_data)
# p1 <- make_samp_plot(xb5)
# 
# 
# xb3 <- prefsamp_station(3,mean_by_station, true_data)
# p2 <- make_samp_plot(xb3)
# #
# #
# xb1 <- prefsamp_station(1,mean_by_station, true_data)
# p3 <- make_samp_plot(xb1)
# #
# #
# xb0 <- prefsamp_station(0,mean_by_station, true_data)
# p4 <- make_samp_plot(xb0)
# 
# 
# 
# p5 <- ggplot(mean_by_station, aes(x=x_coords, y=y_coords, fill = mean_value)) +   geom_raster() +
#   xlab("x") +
#   ylab("y") +
#   scale_fill_viridis_c(option="viridis",direction=-1,name="Mean Y value") +
#   theme_bw()
# 

# (p1 + p2)/ (p3 + p4)/(plot_spacer() + p5 + plot_spacer())
