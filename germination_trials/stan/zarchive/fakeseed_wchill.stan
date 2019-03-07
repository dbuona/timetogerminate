data {
  int<lower=1> N;

  real t[N]; //this is time 
  real Y[N]; //this is the response
  vector[N] chill; // this is the amount of chilling
} 

parameters {

  real b_chill_beta; //effect of chilling on sliope of S curve
  real b_chill_t50; //effect of chilling on t50.
  real b_chill_d; //effect of chill on max germination
  
  real a_beta;  // alpha for beta
  real a_t50; //alpha for t50
  real a_d;  //alpha for max germination
  real<lower=0>  sigma; //sigma
  
  real<lower=0>  sigma_b_beta; 
  real<lower=0>  sigma_b_d;
  real<lower=0>  sigma_b_t50;
  
  real<lower=0>  sigma_a_beta; 
  real<lower=0>  sigma_a_d;
  real<lower=0>  sigma_a_t50;

} 

transformed parameters {
real y_hat[N];
   for (i in 1:N)
y_hat[i] =(b_chill_d*chill[i]+a_d)/(1+(((t[i])/(b_chill_t50*chill[i]+a_t50))^(b_chill_beta*chill[i]+a_beta))); //log logistic equation where paramenters beta, t50 and d are represent as chilling sub mod
     
     }
 

 model {
  // priors
  //t50 ~ uniform(0, 100); 
  //beta ~ normal(0, 50); 
  //d ~ uniform(0, 30); 
  a_beta~normal(0,10);
  a_t50~normal(0,10);
  a_d ~ normal(0,10);
  b_chill_beta ~normal(0,1);
  b_chill_t50 ~normal(0,1);
  b_chill_d ~normal(0,1);
   sigma ~ normal(0, 10);
  
  sigma_b_beta ~ normal(0, 1);
  sigma_b_d ~ normal(0, 1);
  sigma_b_t50 ~ normal(0, 1);
  
   sigma_a_beta ~ normal(0, 10);
  sigma_a_d ~ normal(0, 10);
  sigma_a_t50 ~ normal(0, 10);
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

