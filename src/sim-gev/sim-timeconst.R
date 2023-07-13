" Simulate a Spatial Process for the max using a GEV "
set.seed(123)
library(tidyverse)
library(MASS)
library(fields)
library(extRemes)
library(VGAM)
source(here::here("src/utils.R"))

SCALE <- 1.5
SHAPE <- 0

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

mu <- BETA1 * X + eta

df <- data.frame(
  sapply(
    mu, function(x) rgev(MONTHS,location=x,scale=SCALE,shape=SHAPE))
  ) %>%
  pivot_longer(
    everything(),
    names_to = "station"
) %>%
  mutate(
    station=parse_number(station),
  ) %>%
  arrange(station) %>%
  mutate(
    measurements= rep(1:MONTHS,N),
    x = rep(X,each=MONTHS),
    x_coords = rep(x_coords,each=length(x_coords)*MONTHS),
    y_coords= rep(rep(y_coords, each = MONTHS),length(y_coords)),
  )

df %>%
  filter(station %in% 13:20) %>%
  ggplot(aes(x=value)) +
  geom_density() +
  facet_wrap(~station) +
  theme_bw() +
  xlab("y") +
  ggtitle("distribution of 60 monthly max by station \ntruth data from Gumbel Distribution")
# 

write.csv(
  df,
  file=here::here("data/sim-gev-2.csv")
)