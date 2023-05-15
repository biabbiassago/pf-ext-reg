
SIGMA2 <- 1
PHI <- 0.2
SIGMA2X <- 3
PHIX <- 0.5
BETA1 <- 10
OBS <- 30
MONTHS <- 50
SAMPLE_SIZE = 30


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
