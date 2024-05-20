## Reparametrize GEV
library(SpatialExtremes)
elle <- function(a,xi){
  # auxilliary function 
  if(xi !=0){
    return((-log(a))^(-xi))
  }
  if(xi == 0){
    return(log(-log(a)))
  }
}

rep2gev <- function(q, s, xi, a=0.5, b=0.5){
  # transforms REP-GEV parameters to regular GEV parameters
  if(xi != 0){
    den <- elle((1-b/2),xi)-elle(b/2,xi)
    mu <- q - s*(elle(a,xi)-1)/den
    sigma <- xi*s/den
  
  }
  if(xi == 0){
    den <- elle(b/2,0) - elle(1-b/2,0)
    mu <- q+ s*elle(a,0)/den
    sigma <- s / den 
  }
  return(list("mu"=mu,"sigma"=sigma,"xi"=xi))
}

gev2rep <- function(mu, sigma, xi, a=0.5, b=0.5){
  # transforms regular GEV parameters to REP-GEV parameters
  if(xi != 0){
    den <- elle((1-b/2),xi)-elle(b/2,xi)
    q <- mu + sigma*(elle(a,xi)-1)/xi
    s <- sigma*den/xi
  }
  if(xi == 0){
    den <- elle(b/2,xi) - elle(1-b/2,x)
    q <- mu - sigma*(elle(a,xi))
    s <- sigma*den

  }
  return(list("q"=q,"s"=s,"xi"=xi))
}


rgevrep <- function(N,q,s,xi,a=0.5,b=0.5){
  # given q, s, xi REP-GEV parameters, simulate iid GEV
  gev_pars <- rep2gev(q,s,xi,a,b)
  return(
    SpatialExtremes::rgev(N, loc=gev_pars$mu,scale=gev_pars$sigma,shape=gev_pars$xi)
  )
}

qgevrep <- function(quant,q,s,xi,a=0.5,b=0.5){
  # given a quantile q, s, xi REP-GEV parameters, find quantile fnc 
  gev_pars <- rep2gev(q,s,xi,a,b)
  return(
    SpatialExtremes::qgev(quant, loc=gev_pars$mu,scale=gev_pars$sigma,shape=gev_pars$xi)
  )
}


## check
make_check_plot <- function(median,iqr,xi){
  N <- seq(0.1,0.99,0.05)
  new_pars <- list(median = median, iqr = iqr, xi = xi)
  gev_pars <- rep2gev(new_pars$median, new_pars$iqr, new_pars$xi)
  reg_gev <- SpatialExtremes::qgev(N, loc=gev_pars$mu,scale=gev_pars$sigma,shape=gev_pars$xi)
  rep_gev <- qgevrep(N, new_pars$median, new_pars$iqr, new_pars$xi)
  
  df_plot <- data.frame("RegularGEV"=reg_gev,"Reparametrized"=rep_gev)
  
  df_plot %>% 
    pivot_longer(everything(),names_to="Parametrization") %>%
    ggplot(aes(x = value, fill = Parametrization)) +
    geom_density(alpha = 0.3) + 
    theme_classic()
}