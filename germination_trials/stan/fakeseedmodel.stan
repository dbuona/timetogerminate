data {
  int<lower=1> N; 
  real t[N]; //this is time 
  real Y[N]; //this is the response
} 

parameters {
  real <lower=0,upper=20> d; 
  real beta;  
  real<lower=0> t50; 
  real<lower=0>  sigma;
} 

transformed parameters {
  real y_hat[N];
  for (i in 1:N) 
    y_hat[i] = d/(1+((t[i]/t50)^beta));
    
 } 
 
 model {
  // priors
  t50 ~ uniform(0, 100); 
  beta ~ normal(0, 50); 
  d ~ uniform(0, 30); 
   sigma ~ normal(0, 10);
  // likelihood
  Y ~ normal(y_hat, sigma);
}

//generated quantities{
  //real Y_mean[N]; 
 //real Y_pred[N];
  //for(i in 1:N){
    // Posterior parameter distribution of the mean
  //  Y_mean[i] = d/(1+((t[i]/t50)^beta));
    // Posterior predictive distribution
  // Y_pred[i] = normal_rng(Y_mean[i], sigma);   
//}
//}

