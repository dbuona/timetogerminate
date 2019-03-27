  //based on a poor understanding of https://groups.google.com/forum/#!topic/stan-users/OexQuNGnsko
  //other tries were with https://magesblog.com/post/2015-10-27-non-linear-growth-curves-with-stan/
  
  functions { 

    real loglogistic(real t, real t50, real beta) { 
      return 1/(1+(t/t50)^-beta); 
   
    // loglogistic_cdf_log(real t, real t50, real beta) 
    //return -log1p_exp(-beta * (log(t) - log(t50)); 
  } 
   
    }


 data {
  int<lower=1> N;//I think this is the number of replicates
  real t[N]; //this is time 
  real Y[N]; //this is the response
  vector[N] chill; // this is the amount of chilling
 // number of replicates
} 

parameters {

  real b_chill_beta; //effect of chilling on sliope of S curve
  real b_chill_t50; //effect of chilling on t50.
 // real b_chill_d; //effect of chill on max germination
  
  real a_beta;  // alpha for beta
  real a_t50; //alpha for t50
  //real a_d;  //alpha for max germination
  real<lower=0>  sigma; //sigma
  

} 

transformed parameters {
real t50; 
//real d;
real beta;
//real y_hat[N];
  // for (i in 1:N)
//y_hat[i]=(b_chill_d*chill[i]+a_d)/(1+(((t[i])/(b_chill_t50*chill[i]+a_t50))^(b_chill_beta*chill[i]+a_beta)));
//log logistic equation where paramenters beta, t50 and d are represent as chilling sub mod not in current use due to function block

for (i in 1:N)
t50=b_chill_t50*chill[i]+a_t50;

//for (i in 1:N)
//d=b_chill_d*chill[i]+a_d;

for (i in 1:N)
beta=b_chill_beta*chill[i]+a_beta;
     }
 

 model {
  // priors
  //t50 ~ uniform(0, 100); 
  //beta ~ normal(0, 50); 
 // d ~ uniform(0, 50); 
  a_beta~normal(0,10);
  a_t50~uniform(0,20);
  //a_d ~ uniform(0,10);
  b_chill_beta ~normal(0,10);
  b_chill_t50 ~normal(0,10);
  //b_chill_d ~normal(0,10);
  
  sigma ~ normal(0, 1);
 
 
  // likelihood
  //Y ~ logistic(y_hat, sigma);

  for (i in 1:N)
  target += loglogistic(t[i],t50,beta); 
}
