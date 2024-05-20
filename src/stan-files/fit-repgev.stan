// quello che hai usato per le simulazioni!!!!!1
// This Stan program defines the model
functions{
    int comp_above_size(vector x, real thresh){
        int N = rows(x);
        int comp;
        comp=0;
        for(k in 1:N){
          comp = comp + (x[k] > thresh);
        }
        return comp;
     }
    vector vector_above(vector x,real thresh, int comp){
        int N = rows(x);
        vector[comp] x_red;
        int pos = 1;
        for(j in 1:N){
          if(x[j]>thresh){
            x_red[pos] = x[j];
            pos = pos + 1;
          }
        }
        return x_red;
      }
      
     int comp_below_size(vector x, real thresh){
        int N = rows(x);
        int comp;
        comp=0;
        for(k in 1:N){
          comp = comp + (x[k] < thresh);
        }
        return comp;
     }
      
      vector vector_below(vector x,real thresh,int comp){
        int N = rows(x);
        vector[comp] x_red;
        int pos = 1;
        
        for(j in 1:N){
          if(x[j]<thresh){
            x_red[pos] = x[j];
            pos = pos + 1;
          }
        }
        return x_red;
      }
      
  
      real gev_lpdf(vector x, real mu, real sigma, real xi) {
        
        vector[rows(x)] t;
        vector[rows(x)] lp;
        real thresh;
        int N = rows(x);
        int comp = N;
        vector[comp] x_red;

        
        
        thresh = mu-(sigma)/xi;

        if(xi < 0 ){
          
          comp = comp_below_size(x, thresh);
          x_red = vector_below(x, thresh, comp);
          
        }
        else if (xi > 0){
          comp = comp_above_size(x, thresh);
          x_red = vector_above(x, thresh, comp);

        }
        else {
          comp = N;
          x_red = x;
        }
        
      
        for(n in 1:comp){
            t[n] = abs(xi) < 1e-10 ? exp((mu - x_red[n]) / sigma) : pow(1 + xi * ((x_red[n] - mu ) / sigma), -1/xi);
            lp[n] = -log(sigma) + (xi + 1) * log(t[n]) - t[n];
        }
        return sum(lp);
      }
  real elle(real quant, real xi){
    real elle_value;
    elle_value = pow(-log(quant), -xi);
    return elle_value;
  }
  real elle_star(real quant){
    real elle_star_value;
    elle_star_value = log(-log(quant));
    return elle_star_value;
  }
  //real loc_to_median(real mu, real sigma, real xi){
    //if(xi == 0){
      //q = mu - sigma*elle_star(0.05);
    //}
    //else{
      //q = mu + sigma*(elle(0.05,xi)-1)/xi
    //}
    //return q
  //}
  
  real median_to_loc(real q, real iqr, real xi){
    real mu;
    if(xi == 0){
      mu = q + iqr*(elle_star(0.5))/(elle_star(0.5/2)-elle_star(1-0.5/2));
    }
    else{
      mu = q - iqr*(elle(0.5,xi)-1)/(elle(1-0.5/2,xi) - elle(0.5/2,xi));
    }
    return mu;
  }
  
  real iqr_to_scale(real iqr, real xi){
    // reparametrization from iqr to scale par
    real sigma;
    if(xi==0){
      sigma = iqr / (elle_star(0.5/2)-elle_star(1-0.5/2));
    }
    else{
      sigma = xi*iqr / (elle(1-0.5/2,xi)-elle(0.5/2,xi));
    }
    return sigma;
  }
  
  
  // exponential covariance function def
  matrix exp_cov(matrix DMat, real sigma_sq, real scale) {
    int S = dims(DMat)[1];
    matrix[S, S] K;
    for (i in 1:(S-1)) {
      K[i, i] = sigma_sq;
      for (j in (i + 1):S) {
        K[i, j] = sigma_sq * exp(- DMat[i,j] / scale);
        K[j, i] = K[i, j];
      }
    }
    K[S, S] =sigma_sq;
    return K;
  }
  vector return_level(vector mu, vector sigma, real xi, real nyear){
    int s = rows(mu);
    
    vector[s] q;
    q = mu - (sigma/xi) * (1-pow(-log(1-1/nyear),-xi));
    return q;
  }
  
}
data {
  int s;
  int months;
  int N;
  vector[N] z;
  matrix[s, s] DMat; // Distance matrix
}
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  //real<lower=0> gamma0;
  real<lower=0> gamma1;
  
  vector[s] eta;
  vector[s] nu;
  
  
  real alpha0;
  real rho0;
  real xi;
}
transformed parameters{
  
  matrix[s,s] SIGMAETA;
  matrix[s,s] SIGMANU;
  matrix[s,s] L_SE;
  matrix[s,s] L_SN;
  
  vector[s] mu;
  vector[s] sigma;
  vector[s] q;
  vector[s] iqr;
    
  SIGMAETA = exp_cov(DMat, beta0, beta1);
  SIGMANU = exp_cov(DMat,1,gamma1);
  
  L_SE = cholesky_decompose(SIGMAETA);
  L_SN = cholesky_decompose(SIGMANU);
 
  // data model

  q = alpha0 + eta;
  iqr = exp(rho0 + nu);
  
  for(i in 1:s){
    mu[i] = median_to_loc(q[i], iqr[i],xi);
  }
  for(j in 1:s){
    sigma[j] = iqr_to_scale(iqr[j],xi);
  }
}
model {
  
  int pos; 
  
  // priors
  target += normal_lpdf(alpha0|0,10);
  target += normal_lpdf(rho0| 0,10);
  
  // priors on variance and range of the location par
  // true values are beta0 0.2 and gamma1 0.4
  // putting a very informative prior on the 
  // range parameter
  target += inv_gamma_lpdf(beta0 | 6, 1);
  target += gamma_lpdf(beta1 | 80, 200);
  

  //target += inv_gamma_lpdf(gamma0 | 6, 0.5);
  target += gamma_lpdf(gamma1 | 60, 200);
  
  
  target += normal_lpdf(xi | 0, 0.5);

  target += multi_normal_cholesky_lpdf(eta | rep_vector(0,s) , L_SE);
  target += multi_normal_cholesky_lpdf(nu |rep_vector(0,s), L_SN);
  
  
  // likelihood
  pos = 1;
  for (i in 1:s) {
    segment(z, pos, months) ~ gev(mu[i], sigma[i],xi);
    pos = pos + months;
  }
}
generated quantities{
  vector[s] ret_10;
  ret_10 = return_level(mu, sigma, xi, 10.0);
}


