true_return_levels <- function(fd,nyear){
  mus <- fd$true_mus
  sigmas <- fd$true_sigmas
  xi <- fd$true_xi
  return(return_levels(mus,sigmas,xi,nyear))
}
return_levels <- function(mu, sigma, xi, year){
  return(mu - (sigma/xi)*(1-(-log(1-1/year))^(-xi)))
}