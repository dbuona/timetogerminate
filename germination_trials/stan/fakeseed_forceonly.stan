data {
  int<lower=1> N;

  real t[N]; //this is time 
  real Y[N]; //this is the response
  vector[N] warm; // this is the amount of fforcing
} 

parameters {

  real b_warm_beta; //effect of chilling on sliope of S curve
  real b_warm_t50; //effect of chilling on t50.
  real b_warm_d; //effect of chill on max germination
  
  real a_beta;  // alpha for beta
  real a_t50; //alpha for t50
  real a_d;  //alpha for max germination
  real<lower=0>  tau; //sigma
  

} 

transformed parameters {
real y_hat[N];
real sigma; 

sigma= 1 / sqrt(tau);
   for (i in 1:N)
y_hat[i] =(b_warm_d*warm[i]+a_d)/(1+(((t[i])/(b_warm_t50*warm[i]+a_t50))^(b_warm_beta*warm[i]+a_beta))); //log logistic equation where paramenters beta, t50 and d are represent as chilling sub mod
     
     }
 

 model {
  // priors
  //t50 ~ uniform(0, 100); 
  //beta ~ normal(0, 50); 
  //d ~ uniform(0, 30); 
  a_beta~normal(0,10);
  a_t50~normal(0,10);
  a_d ~ normal(0,10);
  b_warm_beta ~normal(0,1);
  b_warm_t50 ~normal(0,1);
  b_warm_d ~normal(0,1);
  
 
 
  // likelihood
  Y ~ normal(y_hat, sigma);
}
