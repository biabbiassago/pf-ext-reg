" Simulate a Gaussian Spatial Process with Exponential covariance fnct "
set.seed(4649)
library(tidyverse)
library(MASS)
library(fields)

SIGMA2 <- 1
PHI <- 2
BETA1 <- 10
OBS <- 30

# define exponential covariance function
exp_cov <- function(sigma2,phi, distp){
    return(
      (sigma2*exp(-(distp/phi)))
    )
}

# generate locations on 10X10 grid
x_coords <- seq(0.1,1,by=0.1)
y_coords <- seq(0.1,1,by=0.1)
n_coords <- length(x_coords)

coords <- expand.grid(x_coords,y_coords)
N <- dim(coords)[1]
# plot(coords[,1],coords[,2])

distance_mat <- fields::rdist(coords)
# dim(distance_mat)

# todo: fix so you don't do double calcs...
Sigma_exp <- matrix(0,nrow=N,ncol=N)
for(i in 1:N){
  for(j in 1:N){
    Sigma_exp[i,j] <- exp_cov(SIGMA2,PHI,distance_mat[i,j])
  }
}

X <- mvrnorm(
  1,
  mu = rep(1,N),
  Sigma=Sigma_exp
)

eta <- mvrnorm(
  1,
  mu = rep(0,N),
  Sigma=Sigma_exp
)

## check?????

rep_data <- matrix(NA, nrow=N, ncol=OBS)

for(j in 1:OBS){
  gau_process <- mvrnorm(
    1,
    mu=BETA1*X+eta,
    Sigma=Sigma_exp
  )
  rep_data[,j] <- gau_process
}
names(rep_data)<- c(1:OBS)

df <- data.frame(
    x_coords = x_coords,
    y_coords=y_coords,
    rep_data,
    x = X
  ) %>%
  pivot_longer(
    cols=-c("x_coords","y_coords","x"),
    names_to = "measurements"
  ) %>%
  mutate(measurements = parse_number(measurements)) %>%
  mutate(station = rep(1:100,each=OBS))


par(mfrow=c(3,2))

for(i in c(1:6)){
  image.plot(
    x_coords,
    y_coords, 
    matrix(rep_data[,i], n_coords,n_coords),
    col=terrain.colors(100)
  ) + title(paste("True Data, Meas:",i))
}

head(df)


write.csv(
  df,
  file="data/sim1.csv"
)
