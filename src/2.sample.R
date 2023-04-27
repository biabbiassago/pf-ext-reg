library(tidyverse)
library(patchwork)

source("src/1.sim.R")

set.seed(3523)
#NSIMS<-100
SAMPLE_SIZE = 30

true_data <- read.csv("data/sim1.csv")

mean_by_station <- true_data %>%
  group_by(station) %>%
  summarize(
    mean_value = mean(value),
    x_coords = first(x_coords),
    y_coords = first(y_coords),
    x = first(x)
)

generate_sampling_probs <- function(b,dat){
  S_eta <- SIGMA2 * exp(-(distance_mat/PHI))
  n <- length(dat$x)
  eta <- MASS::mvrnorm(
    n = 1,
    mu = rep(0,n),
    Sigma=S_eta
  )
  ### check
  #pi <- sapply(BETA1*dat$x + eta, function(x) max(0,x))
  pi <- abs(b*dat$mean_value + b*eta)
  pi <- pi/sum(pi)
  return(pi)
}
pi <- generate_sampling_probs(3,mean_by_station)

sample_by_stations <- function(dat,prob){
  ### Taking the Sample ####
  s_i <- sample(
      1:N,
      size=SAMPLE_SIZE,
      replace=FALSE,
      prob=prob
  )
  
  # vector of sampled stations
  #r_s <- ifelse(1:N %in% s_i, 1, 0)
  
  sample_df <- dat %>% filter(station %in% s_i)
  return(list("sample_df"=sample_df,"stations"=s_i))
}

x <- sample_by_stations(true_data,pi)
sample_df <- x$sample_df
sampled_stations <- x$stations


### EXAMPLE PLOT ####

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

p1 + p2 + plot_annotation(title="Sample with b=3 and full data")

