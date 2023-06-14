
# sampling probabilities

true_data <- read.csv("data/sim2.csv")
mean_by_station <- true_data %>%
  group_by(station) %>%
  summarize(
    mean_value = mean(value),
    x_coords = first(x_coords),
    y_coords = first(y_coords),
    x = first(x)
  )


generate_sampling_probs <- function(b,dat, type="power"){
  n <- length(dat$x)
  
  eps <- rnorm(n,0,1)
  if(b != 0){
    
    if(type=="linear"){
      pi <- abs(b*dat$mean_value+eps)
    }
    if(type=="quadratic"){
      pi <- (b*dat$mean_value+eps)^2
    }
    if(type=="power"){
      pi <- (abs(dat$mean_value)+eps)^(2*b)
    }
    if(type=="exponential"){
      pi <- exp(abs(b*dat$mean_value) + eps)
    }
    pi_std <- pi/sum(pi)
  }
  else if(b==0){
    pi_std <- rep(1/n,n)
  }
  return(list("pi"=pi_std,"b"=b))
}



plot(mean_by_station$mean_value,generate_sampling_probs(b=5,mean_by_station, type="quadratic")$pi,xlab="Mean Value by station",ylab="Probability",col="darkgreen")
points(mean_by_station$mean_value, generate_sampling_probs(b=3,mean_by_station, type="quadratic")$pi,col="red")
points(mean_by_station$mean_value,generate_sampling_probs(b=1,mean_by_station, type="quadratic")$pi,col="blue")
points(mean_by_station$mean_value,generate_sampling_probs(b=0,mean_by_station, type="quadratic")$pi,col="darkgray")
ylim = c(0,0.5)
legend("topleft",legend=c(5,3,1,0),col=c("darkgreen","red","blue","darkgray"),pch=1)



plot(mean_by_station$mean_value,generate_sampling_probs(b=5,mean_by_station, type="linear")$pi,xlab="Mean Value by station",ylab="Probability",col="darkgreen")
points(mean_by_station$mean_value, generate_sampling_probs(b=3,mean_by_station, type="linear")$pi,col="red")
points(mean_by_station$mean_value,generate_sampling_probs(b=1,mean_by_station, type="linear")$pi,col="blue")
points(mean_by_station$mean_value,generate_sampling_probs(b=0,mean_by_station, type="linear")$pi,col="darkgray")
ylim = c(0,0.5)
legend("topleft",legend=c(5,3,1,0),col=c("darkgreen","red","blue","darkgray"),pch=1)





plot(mean_by_station$mean_value,generate_sampling_probs(b=5,mean_by_station, type="power")$pi,xlab="Mean Value by station",ylab="Probability",col="darkgreen")
points(mean_by_station$mean_value, generate_sampling_probs(b=3,mean_by_station, type="power")$pi,col="red")
points(mean_by_station$mean_value,generate_sampling_probs(b=1,mean_by_station, type="power")$pi,col="blue")
points(mean_by_station$mean_value,generate_sampling_probs(b=0,mean_by_station, type="power")$pi,col="darkgray")
ylim = c(0,0.5)
legend("topleft",legend=c(5,3,1,0),col=c("darkgreen","red","blue","darkgray"),pch=1)




