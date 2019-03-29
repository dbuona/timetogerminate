//based on a poor understanding of https://groups.google.com/forum/#!topic/stan-users/OexQuNGnsko
  //other tries were with https://magesblog.com/post/2015-10-27-non-linear-growth-curves-with-stan/
  
  functions { 
    
    real loglogistic(real t, real t50, real beta) { 
      return 1/(1+(t/t50)^-beta); 
      
      // loglogistic_cdf_log(real t, real t50, real beta) 
      //return -log1p_exp(-beta * (log(t) - log(t50)); 
    } 
    
  }


data {
  int<lower=1> N;//I think this is the number of replicates
  real t[N]; //this is time 
  real Y[N]; //this is the response

} 

parameters {
  real d;
  real beta;
  real t50; 
  real<lower=0>  sigma; //sigma
  
  
} 

transformed parameters {
}

model {
  // priors
  t50 ~ uniform(0, 50); 
  beta ~ normal(0, 20); 
   d ~ normal(1,1); 

  sigma ~ normal(0, 1);
  
  
  // likelihood
  //Y ~ logistic(y_hat, sigma);
  
  for (i in 1:N)
    target += loglogistic(t[i],t50,beta); 
}
