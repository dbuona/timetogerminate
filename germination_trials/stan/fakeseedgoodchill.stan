data {
  int<lower=1> N; 
  real t[N]; //this is time 
  real Y[N]; //this is the response
  vector[N] chill;// chilling treatmet (0,1)
} 

parameters {
  real <lower=0,upper=20> d; 
   real<lower=0> a_beta;
   real b_beta;
   real<lower=0> a_t50;
   real b_t50;
 
  real<lower=0>  sigma;
} 

transformed parameters {
  real beta;
  real<lower=0> t50;
  real y_hat[N];
  
    for (i in 1:N) 
  beta=a_beta+b_beta*chill[i];
    for (i in 1:N) 
  t50=a_t50+b_t50*chill[i];
  
  for (i in 1:N) 
    y_hat[i] = d/(1+((t[i]/t50)^-beta));
    
 } 
 
 model {
  // priors
  //t50 ~ uniform(0, 100); 
  //beta ~ normal(0, 50); 
  d ~ logistic(0,1); 
   sigma ~ normal(0,1);
  // likelihood
  Y ~ normal(y_hat, sigma);
}
