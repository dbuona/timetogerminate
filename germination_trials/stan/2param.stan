data {
  int<lower=1> N; 
  real t[N]; //this is time 
  vector<lower=0>[N] Y; //this is the response
  vector[N] chill;// chilling treatmet (0,1)
} 

parameters {
  real<lower=0> a_beta;
  real<lower=0> b_beta;
  real<lower=0> a_t50;
  real b_t50;
  real<lower=0, upper=1>  sigma;
} 

transformed parameters {
  vector<lower=0, upper=1>[N] y_hat;
  
  for (i in 1:N)
    y_hat[i] = (1)/(1+exp(-(a_beta+b_beta*chill[i]) * (log(t[i]) - log(a_t50+b_t50*chill[i]))));
} 

model {
  // priors
  a_t50 ~ normal(30, 15); // fake
  b_t50 ~ normal(-3, 3);
  a_beta ~ normal(1, 4); //fake data
  b_beta ~ normal (0,2);

  //sigma ~ normal(0,.1); //sigma for fake data
  sigma ~ normal(0, .3) ; //real plants
  // likelihood
  Y ~ normal(y_hat, sigma);
}

generated quantities {
  vector[N] Y_mean; 
  vector[N] Y_pred; 
  for(i in 1:N){
    // Posterior parameter distribution of the mean
    Y_mean[i] = (1)/(1+exp(-(a_beta+b_beta*chill[i]) * (log(t[i]) - log(a_t50+b_t50*chill[i]))));
    
    Y_pred[i] = normal_rng(Y_mean[i], sigma); // Posterior predictive distribution 
  }
}
