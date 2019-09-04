data {
  int<lower=1> N; 
  vector[N] Y; //this is the response
  vector[N] Den_c;// Desity of species c
  vector[N] Den_e;// Density of species e
  //vector[N] T_c; //average Time of germination for c
  //vector[N] T_e; //average Time of germination for e

}

parameters {
  real <lower=-0,upper=10> alpha;
  real <lower=0,upper=20> beta_c;
  real <lower=0,upper=10> beta_e;
  //real <lower=-1,upper=1> B;
  real <lower=0> sigma;
}

transformed parameters {
//real ratio[N];
//for (i in 1:N)
//ratio[i]=T_c[i]./T_e[i];
  
  vector<lower=0>[N] y_hat;
  
  for (i in 1:N)
     y_hat[i] = (alpha+beta_c*Den_c[i]+beta_e*Den_e[i]);
 }
 //*(pow((T_c[i]./T_e[i]),B)
  model {
  // priors
  alpha ~ normal(0, 5); 
  beta_c ~ normal(0, 10);
  beta_e ~ normal(0,10);
 // T_c ~ normal(0,5);
 // T_e ~ normal (0,5);
 // B ~ uniform(-1,1);

  sigma ~ normal(0,1);
  // likelihood
  Y ~ normal(y_hat, sigma);
}
