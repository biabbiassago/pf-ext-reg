fit_gev <- '
  functions{
      real vector_comp_less(vector x,real thresh){
        int N = rows(x);
        real tmp;
        real comp;
        comp=0;
        for(k in 1:N){
          tmp = x[k] < thresh; 
          comp = comp + tmp;
        }
        // returns 
        return N != comp;
      }
      real vector_comp_more(vector x,real thresh){
        // returns TRUE if any elements is BELOW the threshold
        // sample CANNOT be from xi>0
        int N = rows(x);
        real tmp;
        real comp; // number of elements that are above the threshold
        comp=0;
        for(k in 1:N){
          tmp = x[k] > thresh; 
          comp = comp + tmp;
        }
        return N != comp;
      }
  
      real gev_lpdf(vector x, real mu, real sigma, real xi) {
        vector[rows(x)] t;
        vector[rows(x)] lp;
        int N;
        real thresh;
        thresh = (mu-sigma)/xi;
        N = rows(x);
        
        real comp_less;
        comp_less = vector_comp_less(x,thresh);

        real comp_more;
        comp_more = vector_comp_more(x,thresh);

        if(comp_less && xi<0){
          reject("sample cannot be generated from such xi ",xi);
          
          
              //for(n in 1:N){
                //lp[n] = 0
              //}
        }
        
        else if(comp_more && xi>0){
          reject("sample cannot be generated from such xi",xi);
          
          
             //for(n in 1:N){
                //lp[n] = 0
              //}
        }
        
        else{
          for(n in 1:N){
            t[n] = abs(xi) < 1e-10 ? exp((mu - x[n]) / sigma) : pow(1 + xi * ((x[n] - mu ) / sigma), -1/xi);
            lp[n] = -log(sigma) + (xi + 1) * log(t[n]) - t[n];
          
          }
          
        }
      
        return sum(lp);
      }
      
  }
  data{
   int N;
   vector[N] x;
  }
  parameters{
    real mu;
    real<lower=0> sigma;
    real xi;
  }
  model{
    
    mu ~ normal(0,2);
    sigma ~ inv_gamma(3,1);
    xi ~ normal(0,0.25);
    
    
    x ~ gev(mu,sigma,xi);
    
  }
'

N <- 1000
gev_vec <- rgev(N, location = 0, scale = 1.5,shape = 0.1)
gev_fit <- stan(
  model_code = fit_gev,
  data = list(N=1000, x = gev_vec),
  chains = 1,
  warmup = 100,
  iter = 1000,
  control = list(adapt_delta = 0.999, max_treedepth=13)
)


# if the prior on xi is really sparse it does really converge. 

