data {
  int<lower=1> N; 
  real t[N]; //this is time 
  real Y[N]; //this is the response
} 

parameters {
  real <lower=0,upper=1> d; 
  real <lower=0>beta;  
  real<lower=0> t50; 
  real<lower=0>  sigma;
} 

transformed parameters {
  real y_hat[N];
  for (i in 1:N) 
    //y_hat[i] = d/(1+((t[i]/t50)^-beta));
     y_hat[i] = d/(1+exp(-beta * (log(t[i]) - log(t50))));
 } 
 
 model {
  // priors
  t50 ~ normal (20,10 )T[0,100]; 
  beta ~ normal(0, 50); 
  d ~ normal(.5, .5)T[0,1]; 
   sigma ~ normal(0,1);
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

