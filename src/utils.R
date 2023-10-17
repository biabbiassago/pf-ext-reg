PHI <- 0.6
SIGMA2X <- 0.2
BETA1 <- 1.2
OBS <- 720 
MONTHS <- 5

if(!exists("SIGMA2")){
  SIGMA2 <- 0.8
}

if(!exists("GAP_SIZE")){
  GAP_SIZE = 0.1
}

if(!exists("PHIX")){
  PHIX = 0.5
}


# define exponential covariance function
exp_cov <- function(sigma2,phi, distp){
  return(
    (sigma2*exp(-(distp/phi)))
  )
}

# generate locations on grid
x_coords <- seq(GAP_SIZE,1.0,by=GAP_SIZE)
y_coords <- seq(GAP_SIZE,1.0,by=GAP_SIZE)

n_coords <- length(x_coords)

true_coords <- expand.grid(x_coords,y_coords)
N <- dim(true_coords)[1] # number of stations
# plot(coords[,1],coords[,2])

distance_mat <- fields::rdist(true_coords)
# dim(distance_mat)

# fix so you don't do double calcs...
Sigma_x <- matrix(0,nrow=N,ncol=N)
for(i in 1:N){
  for(j in 1:N){
    Sigma_x[i,j] <- exp_cov(SIGMA2X,PHIX,distance_mat[i,j])
  }
}


Sigma_eta <- matrix(0,nrow=N,ncol=N)
for(i in 1:N){
  for(j in 1:N){
    Sigma_eta[i,j] <- exp_cov(SIGMA2,PHI,distance_mat[i,j])
  }
}

make_mean_by_station_df <- function(true_data){
  mean_by_station <- true_data %>%
    group_by(station) %>%
    summarize(
      mean_value = mean(value),
      x_coords = first(x_coords),
      y_coords = first(y_coords),
      x = first(x)
    )
  return(mean_by_station)
}

