
SIGMA2 <- 0.1
PHI <- 0.8
SIGMA2X <- 0.2
BETA1 <- 0.5
OBS <- 720 #
MONTHS <- 50


# define exponential covariance function
exp_cov <- function(sigma2,phi, distp){
  return(
    (sigma2*exp(-(distp/phi)))
  )
}

# generate locations on 20X20 grid
x_coords <- seq(0.05,1.0,by=0.05)
y_coords <- seq(0.05,1.0,by=0.05)
n_coords <- length(x_coords)

true_coords <- expand.grid(x_coords,y_coords)
N <- dim(true_coords)[1]
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

