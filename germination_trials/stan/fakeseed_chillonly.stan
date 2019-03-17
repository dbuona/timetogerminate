functions {
 real loglogistic_cdf_log(real t, real t50, real beta) {
    return -log1p_exp(-beta * (log(t) - log(t50))); 
 }

}
data {
  int<lower=1> N;

  int<lower = 0> t[N]; //this is time 
  real<lower = 0> Y[N]; //this is the response
  vector[N] chill; //Level of chilling
} 

parameters {

  real b_chill_beta; //effect of chilling on sliope of S curve
  real b_chill_t50; //effect of chilling on t50.
  //real b_chill_d; //effect of chill on max germination
  
  real a_beta;  // alpha for beta
  real a_t50; //alpha for t50
  //real a_d;  //alpha for max germination
  real<lower=0>  sigma; //sigma
  

} 
transformed parameters {
real t50[N];
real beta[N];
//real d[N];

for (i in 1:N)
t50[i]=a_t50+b_chill_t50*chill[i];
for (i in 1:N)
beta[i]=a_beta+b_chill_beta*chill[i];
//for (i in 1:N)
//d[i]=a_d+b_chill_d*chill[i];
}
model{
 a_beta~normal(0,10);
  a_t50~normal(0,10);
  //a_d ~ normal(0,10);
  b_chill_beta ~normal(0,10);
  b_chill_t50 ~normal(0,10);
 // b_chill_d ~normal(0,10);
  
 t ~ loglogistic(t50, beta); 
 
 
  // likelihood


}
