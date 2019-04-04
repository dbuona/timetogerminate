  //based on a poor understanding of https://groups.google.com/forum/#!topic/stan-users/OexQuNGnsko
  //other tries were with https://magesblog.com/post/2015-10-27-non-linear-growth-curves-with-stan/
  
  functions { 
    real loglogistic_lcdf(vector Y,real d, real t50, real beta){ 
    return d/ -log1p_exp(-beta * (log(t) - log(t50))); 
  } 
  }

 data {
  int<lower=1> N;//I think this is the number of replicates
  real <lower=0, upper=24> t; //this is time 
  vector[N] Y; //this is the response
  vector[N] chill; // this is the amount of chilling 0,1
 // number of replicates
} 

parameters {

  //real b_chill_beta; //effect of chilling on sliope of S curve
  real b_chill_t50; //effect of chilling on t50.
 // real b_chill_d; //effect of chill on max germination

  real<lower=0> a_beta;  // alpha for beta
  real<lower=0> a_t50; //alpha for t50
  real a_d;  //alpha for max germination
  real<lower=0>  sigma; //sigma
} 

transformed parameters {
//real<lower=0> scale;
//real shape;

real t50=b_chill_t50*chill[N]+a_t50;

//real beta=b_chill_beta*chill[N]+a_beta;

//scale=log(t50); 
//shape=(1/beta); 

 
}
 model{
  a_beta~normal(0,10) T[0,];
  a_t50~normal(0,20) T[0,];
 //b_chill_beta ~normal(0,5);
 b_chill_t50 ~normal(0,5);
sigma ~ normal(0, 1);
 
// likelihood
 loglogistic_lcdf(Y|,a_d,t50,a_beta); 

//Y~logistic(scale,shape); //wikipedia say this is an alternative parameterization
}
