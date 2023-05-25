" Simulate a Spatial Process for the MAX- with Exponential covariance fnct "
set.seed(123)
library(tidyverse)
library(MASS)
library(fields)
source("src/utils.R")

X <- mvrnorm(
  1,
  mu = rep(1,N),
  Sigma=Sigma_x
)

eta <- mvrnorm(
  1,
  mu = rep(0,N),
  Sigma=Sigma_eta
)

rep_data <- matrix(NA, nrow=N, ncol=MONTHS)

for(j in 1:MONTHS){
  
  # generate 30 OBS "per day" per station (100 stations)
  gau_process <- mvrnorm(
    OBS,
    mu=BETA1*X+eta,
    Sigma=diag(N)
  )
  # get the max by each station (100 stations) for that "month"
  max_obs <- apply(gau_process, 2, max)
  rep_data[,j] <- max_obs
}

df <- data.frame(
    x_coords = rep(x_coords, each=n_coords),
    y_coords= rep(y_coords, n_coords),
    rep_data,
    x = X
  ) %>%
  pivot_longer(
    cols=-c("x_coords","y_coords","x"),
    names_to = "measurements"
  ) %>%
  mutate(measurements = parse_number(measurements)) %>%
  mutate(station = rep(1:N,each=MONTHS))


# plot examples, by staton for all measurements 
# just a few stations
df %>%
  filter(station %in% 1:12) %>%
  ggplot(aes(x=value)) + 
  geom_density() + facet_wrap(~station) + theme_bw() + xlab("Y") + ggtitle("Distribution of monthly max by station \nTruth Data") 

# quick check for nine month
par(mfrow=c(3,3))
for(i in c(1:9)){
  fields::image.plot(
    x_coords,
    y_coords, 
    matrix(rep_data[,i], n_coords,n_coords),
    col=terrain.colors(100)
  ) + title(paste("True Data, measurements:",i))
}

write.csv(
  df,
  file="data/sim1.csv"
)




## this is for one station???
# library(extRemes)
# station1 <- df %>% filter(station==1)
# fit1 <- fevd(value~1,data=station1)
# # location      scale      shape 
# # 16.5786462  0.4848307 -0.3521654 
# plot(fit1)
# 
# 
# station1 <- df %>% filter(station==1)
# fit1 <- fevd(value~x,data=station1)
# 
# fitoverall <- fevd(value~1,data=df)
# # Estimated parameters:
# #  location     scale     shape 
# # -8.337301  5.059122 -0.070127 
# 
