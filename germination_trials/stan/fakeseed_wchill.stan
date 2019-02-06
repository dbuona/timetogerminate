data {
  int<lower=1> N;
  int<lower=0,M;
  real t[N]; //this is time 
  real Y[N]; //this is the response
  vector[N] chill; // this is whether something got chilling or not
  
} 

parameters {
  real <lower=0,upper=20> d; 
  real beta;  
  real<lower=0> t50; 
  real<lower=0>  sigma;
  real b_chill_beta;
  real b_chill_t50;
  real b_chill_d;
} 

transformed parameters {
  real y_hat[N];
  for (i in 1:N) 
    y_hat[i] = (d*b_chill_d)/(1+((t[i]/(t50*b_chill_t50))^(beta*b_chill_beta)));
    
 } 
 
 model {
  // priors
  t50 ~ uniform(0, 100); 
  beta ~ normal(0, 50); 
  d ~ uniform(0, 30); 
  b_chill_beta ~normal(0,10);
  b_chill_t50 ~normal(0,10);
  b_chill_d ~normal(0,10);
   sigma ~ normal(0, 10);
  // likelihood
  Y ~ normal(y_hat, sigma);
for (i in 1:M) {
    y[i,] ~ normal(y_hat,sigma[i]);
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

