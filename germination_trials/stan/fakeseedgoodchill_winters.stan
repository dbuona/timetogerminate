//based on lizzie's fakeseedgoodchill_alt.2param.stan but adding interactions as of June 25
data {
  int<lower=1> N; 
  real t[N]; //this is time 
  vector<lower=0>[N] Y; //this is the response
  vector[N] chill;// chilling treatmet (0,1,2,4,5)
  vector[N] force;// forcing treatmet (0,1)
} 

transformed data {
  vector[N] inter_cxf;
  inter_cxf = chill .* force;
}

parameters {
  real <lower=0,upper=1> a_d;
  real <lower=0,upper=1> bc_d;
  real <lower=0,upper=1> bf_d;
  
  real <lower=0> a_beta;
  real <lower=0> bc_beta;
  real <lower=0> bf_beta;
  
  real <lower=0> a_t50;
  real bc_t50;
  real bf_t50;
  
  real inter_d;
  real inter_beta;
  real inter_t50; //interactions
  
  
  real <lower=0>  sigma;
} 

transformed parameters {
  vector<lower=0>[N] y_hat;
  
  for (i in 1:N)
     y_hat[i] = (a_d+bc_d*chill[i]+bf_d*force[i]+inter_d*inter_cxf[i])/(1+exp(-(a_beta+bc_beta*chill[i]+bf_beta*force[i]+inter_beta*inter_cxf[i]) * (log(t[i]) - log(a_t50+bc_t50*chill[i]+bf_t50*force[i]+inter_t50*inter_cxf[i]))));
 } 
 
 model {
  // priors
  a_t50 ~ normal(20, 10); // previous prior did not reach to 15 ...
  bc_t50 ~ normal(0, 10);
  bf_t50 ~ normal(0,10);
  a_beta ~ normal(0, 50); 
  bc_beta ~ normal (0,10);
  bf_beta ~ normal (0,10);
  a_d ~ normal(0.5, 0.5); 
  bc_d ~ normal(0.5, 0.5);
  bf_d ~ normal(0.5, 0.5);
  inter_t50 ~ normal(0,10);
  inter_beta ~ normal(0,10);
  inter_d ~ normal (0.5, 0.5); 
  sigma ~ normal(0,1);
  // likelihood
  Y ~ normal(y_hat, sigma);
}

generated quantities {
 vector[N] Y_mean; 
  vector[N] Y_pred; 
  for(i in 1:N){
    // Posterior parameter distribution of the mean
    Y_mean[i] = (a_d+bc_d*chill[i]+bf_d*force[i])/(1+exp(-(a_beta+bc_beta*chill[i]+bf_beta*force[i]) * (log(t[i]) - log(a_t50+bc_t50*chill[i]+bf_t50*force[i]))));
    
    Y_pred[i] = normal_rng(Y_mean[i], sigma); // Posterior predictive distribution 
}
}
