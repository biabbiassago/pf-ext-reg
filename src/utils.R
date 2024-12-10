library(spatstat)

MCMC_ITERS <- 100000

SCALE <- 1.5
SHAPE <- 0.2
SIGMA2ETA <- 0.2
PHIETA <- 0.4

SIGMA2NU <- 1
PHINU <- 0.3

# SIGMA2X <- 0.1
# PHIX <- 0.2
ALPHA0 <- 1
ALPHA1 <- 1.2
RHO0 <- 0.5

PREF_VAR <- 1.2
PREF_SCALE <- 0.7



# define exponential covariance function
exp_cov <- function(sigma2,phi, distp){
  return(
    (sigma2*exp(-(distp/phi)))
  )
}

# if(TYPE == "GRID"){
#   # generate locations on grid
#   x_coords <- seq(GAP_SIZE,1.0,by=GAP_SIZE)
#   y_coords <- seq(GAP_SIZE,1.0,by=GAP_SIZE)
#   
#   n_coords <- length(x_coords)
#   i
#   true_coords <- expand.grid(x_coords,y_coords)
#   N <- dim(true_coords)[1] # number of stations
#   # plot(coords[,1],coords[,2])
# }

sim_random_locs <- function(N){
  dp <- runifpoint(
    N,
    win=owin(c(0,1),c(0,1))
  )
  df <- as.data.frame(dp)
  x_coords <- df$x
  y_coords <- df$y
  distance_mat <- fields::rdist(df)
  return(list("x_coords"=x_coords,"y_coords"=y_coords,"distance_mat"=distance_mat))
}

sim_grid <- function(GAP_SIZE=0.01){
    x_coords <- seq(GAP_SIZE,1.0,by=GAP_SIZE)
    y_coords <- seq(GAP_SIZE,1.0,by=GAP_SIZE)
    n_coords <- length(x_coords)
    true_coords <- expand.grid(x_coords,y_coords)
    names(true_coords) <- c("x_coords","y_coords")
    distance_mat <- fields::rdist(true_coords)
    return(list(
      "x_coords"=true_coords$x_coords,
      "y_coords"=true_coords$y_coords,
      "distance_mat"=distance_mat)
    )
}


# dim(distance_mat)

# fix so you don't do double calcs...
make_Sigma_x <- function(distance_mat,N){
  Sigma_x <- matrix(0,nrow=N,ncol=N)
  for(i in 1:N){
    for(j in 1:N){
      Sigma_x[i,j] <- exp_cov(SIGMA2X,PHIX,distance_mat[i,j])
    }
  }
  return(Sigma_x)
}

make_Sigma_eta <- function(distance_mat,N){
  Sigma_eta <- matrix(0,nrow=N,ncol=N)
  for(i in 1:N){
    for(j in 1:N){
      Sigma_eta[i,j] <- exp_cov(SIGMA2ETA,PHIETA,distance_mat[i,j])
    }  
  }
  return(Sigma_eta)  
}


make_Sigma_nu <- function(distance_mat,N){
  Sigma_nu <- matrix(0,nrow=N,ncol=N)
  for(i in 1:N){
    for(j in 1:N){
      Sigma_nu[i,j] <- exp_cov(SIGMA2NU,PHINU,distance_mat[i,j])
    }  
  }
  return(Sigma_nu)  
}



make_mean_by_station_df <- function(true_data,cov=FALSE){
  if(cov==TRUE){
    mean_by_station <- true_data %>%
      group_by(station) %>%
      summarize(
        mean_value = mean(value),
        x_coords = first(x_coords),
        y_coords = first(y_coords),
        x = first(x)
      )    
  }
  if(cov==FALSE){
    mean_by_station <- true_data %>%
      group_by(station) %>%
      summarize(
        mean_value = mean(value),
        x_coords = first(x_coords),
        y_coords = first(y_coords),
        #x = first(x)
      )
  }
  return(mean_by_station)
}

return_levels <- function(mu, sigma, xi, year){
  return(mu - (sigma/xi)*(1-(-log(1-1/year))^(-xi)))
}
