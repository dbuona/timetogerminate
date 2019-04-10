data {
  int<lower=1> N; 
  vector<lower=0,upper=24>[N] t; //this is time 
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
  real<lower=0> t50;  //
  vector<lower=0>[N] y_hat;
  

    // equiv to: t50 = a_t50 + b_50 * chill[N] // (overwriting)
    // 
    for (i in 1:N) 
      t50=a_t50+b_t50*chill[i];
  
  for (i in 1:N)
//y_hat[i] = d/(1+pow((t[i]/t50),-beta)); or it can be written
y_hat[i]= d /(1+exp(-beta * (log(t[i]) - log(t50)))); 
 } 
 
 model {
  // priors
  a_t50~ normal(0,100)T[0,]; # truncations add CDF, don't need if just adding constants (on right side)
  b_t50 ~ normal (0,100);
  beta ~ normal(0, 100)T[0,]; 
  d ~ normal(0.5,0.5) T[0,1]; # beta(2,2)
   sigma ~ normal(0,1);
   //likelihood  
   
Y~ normal(y_hat, sigma);
}
