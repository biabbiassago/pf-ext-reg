" Simulate a Spatial Process for the MAX- with Exponential covariance fnct "
set.seed(123)
library(tidyverse)
library(MASS)
library(fields)
source("src/utils.R")

X <- mvrnorm(
  MONTHS,
  mu = rep(1,N),
  Sigma=Sigma_x
) %>% as.vector()

eta <- mvrnorm(
  MONTHS,
  mu = rep(0,N),
  Sigma=Sigma_eta
) %>% as.vector()



# generate 24*30 OBS "per hour per day" per station (100 stations)
gau_process <- mvrnorm(
  OBS,
  mu=BETA1*X+eta,
  Sigma=diag(N*MONTHS)
)
max_obs <- apply(gau_process, 2, max)


df <- data.frame(
  x_coords = rep(x_coords,each=length(x_coords)*MONTHS),
  y_coords= rep(rep(y_coords, each = MONTHS),length(x_coords)),
  value=max_obs,
  x = X
) %>%
  mutate(measurements = rep(1:MONTHS,N)) %>%
  mutate(station = rep(1:N,each=MONTHS))


# plot examples, by station for all measurements 
# just a few stations
df %>%
  filter(station %in% 1:12) %>%
  ggplot(aes(x=value)) +
  geom_density() +
  facet_wrap(~station) +
  theme_bw() +
  xlab("y") +
  ggtitle("distribution of monthly max by station \ntruth data")
# 
# # quick check for nine month
# par(mfrow=c(3,3))
# for(i in c(1:9)){
#   fields::image.plot(
#     x_coords,
#     y_coords,
#     matrix(rep_data[,i], n_coords,n_coords),
#     col=terrain.colors(100)
#   ) + title(paste("True Data, measurements:",i))
# }

write.csv(
  df,
  file="data/sim3.csv"
)
# df <- read.csv("data/sim3.csv")




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