data {
  int<lower=1> N; 
  vector<lower=0,upper=24>[N] t; //this is time 
  vector<lower=0>[N] Y; //this is the response
  vector[N] chill;// chilling treatmet (0,1)
} 

parameters {
  //real <lower=0,upper=20> d; 
   real<lower=0> a_beta;
   real b_beta;
   real<lower=0> a_t50;
   real b_t50;
 
  real<lower=0>  sigma;
} 

transformed parameters {
  real<lower=0> beta;
  real<lower=0> t50;
  vector<lower=0>[N] y_hat;
  
    for (i in 1:N) 
  beta=a_beta+b_beta*chill[i];
    for (i in 1:N) 
  t50=a_t50+b_t50*chill[i];
  
  for (i in 1:N) 
    y_hat[i] = 1/(1+((t[i]/t50)^-beta));
  
 } 
 
 model {
  // priors
  //t50 ~ uniform(0, 100); 
  //beta ~ normal(0, 50); 
  //d ~ uniform(0,1); 
   sigma ~ normal(0,1);
  // likelihood
Y~ normal(y_hat, sigma);
}
