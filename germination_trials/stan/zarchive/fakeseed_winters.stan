data {
  int<lower=1> N;

  real t[N]; //this is time 
  real Y[N]; //this is the response
  vector[N] chill; // this is the amount of chilling
  vector[N] inc; //this is the amound of incubation temperature
} 

transformed data {
    vector[N] inter_IC;
    inter_IC    = inc .* chill; 
}

parameters {

  real b_chill_beta; //effect of chilling on sliope of S curve
  real b_chill_t50; //effect of chilling on t50.
  real b_chill_d; //effect of chill on max germination
 
  real b_inc_beta;
  real b_inc_t50;
  real b_inc_d;
  
  real b_inter_beta;
  real b_inter_t50;
  real b_inter_d;
  
  real a_beta;  // alpha for beta
  real a_t50; //alpha for t50
  real a_d;  //alpha for max germination
  
  real<lower=0>  sigma_y; //sigma
  real<lower=0> sigma_a_beta;
   real<lower=0> sigma_a_t50;
   real<lower=0> sigma_a_d;
   
  real<lower=0> sigma_b_chill_beta; 
  real<lower=0> sigma_b_chill_t50; 
  real<lower=0> sigma_b_chill_d;  
  
  real<lower=0> sigma_b_inc_beta; 
  real<lower=0> sigma_b_inc_t50; 
  real<lower=0> sigma_b_inc_d;  
  
  real<lower=0> sigma_b_inter_beta; 
  real<lower=0> sigma_b_inter_t50; 
  real<lower=0> sigma_b_inter_d;  

} 

transformed parameters {
  real y_hat[N];
   for (i in 1:N)
y_hat[i] =(b_chill_d*chill[i]+b_inc_d*inc[i]+b_inter_d*inter_IC[i]+a_d)/(1+(((t[i])/(b_chill_t50*chill[i]+b_inc_t50*inc[i]+b_inter_t50*inter_IC[i]+a_t50))^(b_chill_beta*chill[i]+b_inc_beta*inc[i]+b_inter_beta*inter_IC[i]+a_beta))); //log logistic equation where paramenters beta, t50 and d are represent as chilling sub mod
     
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
  
  b_inc_beta ~normal(0,1);
  b_inc_t50 ~normal(0,1);
  b_inc_d ~normal(0,1);
  
  b_inter_beta ~normal(0,1);
  b_inter_t50 ~normal(0,1);
  b_inter_d ~normal(0,1);
   
   
   sigma_y ~ normal(0, 10);
    sigma_a_beta ~ normal(0, 1);
    sigma_a_t50 ~ normal(0, 1);
   sigma_a_d ~ normal(0, 1);
   
  sigma_b_chill_beta ~ normal(0, 1); 
  sigma_b_chill_t50 ~ normal(0, 1); 
  sigma_b_chill_d ~ normal(0, 1);  
  
  sigma_b_inc_beta ~ normal(0, 1); 
  sigma_b_inc_t50 ~ normal(0, 1); 
  sigma_b_inc_d ~ normal(0, 1);  
  
  sigma_b_inter_beta ~ normal(0, 1); 
  sigma_b_inter_t50 ~ normal(0, 1); 
  sigma_b_inter_d ~ normal(0, 1);  
   
  // likelihood
  Y ~ normal(y_hat, sigma_y);
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

