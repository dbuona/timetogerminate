data {
  int<lower=1> N; 
  real t[N]; //this is time 
  vector<lower=0>[N] Y; //this is the response
  vector[N] chill;// chilling treatmet (0,1)
} 

parameters {
  real <lower=0,upper=1> d; 
  real<lower=0> beta;
  real<lower=0> a_t50;
  real b_t50;
  real<lower=0>  sigma;
} 

transformed parameters {
  vector<lower=0>[N] y_hat;
  
  for (i in 1:N)
     y_hat[i] = d/(1+exp(-beta * (log(t[i]) - log(a_t50+b_t50*chill[i]))));
 } 
 
 model {
  // priors
  a_t50 ~ normal(20, 10); // previous prior did not reach to 15 ...
  b_t50 ~ normal(0, 3);
  beta ~ normal(0, 50); 
  d ~ normal(0.5, 0.5); 
  sigma ~ normal(0,1);
  // likelihood
  Y ~ normal(y_hat, sigma);
}
