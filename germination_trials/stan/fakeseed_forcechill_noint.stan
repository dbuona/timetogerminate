data {
  int<lower=1> N;

  int<lower = 0> t[N]; //this is time 
  real<lower = 0> Y[N]; //this is the response
  int<lower = 0> warm[N]; // this is the amount of fforcing
  int<lower = 0> chill[N]; //Level of chilling
} 

parameters {

  real b_warm_beta; //effect of chilling on sliope of S curve
  real b_warm_t50; //effect of chilling on t50.
  real b_warm_d; //effect of chill on max germination
  
   real b_chill_beta; //effect of chilling on sliope of S curve
  real b_chill_t50; //effect of chilling on t50.
  real b_chill_d; //effect of chill on max germination
  
  real<upper=0> a_beta;  // alpha for beta
  real<lower=0> a_t50; //alpha for t50
  real<lower=0> a_d;  //alpha for max germination
  real<lower=0>  sigma; //sigma
  

} 

transformed parameters {
real y_hat[N];
   for (i in 1:N)
y_hat[i] =(b_warm_d*warm[i]+b_chill_d*chill[i]+a_d)/(1+(((t[i])/(b_warm_t50*warm[i]+b_chill_t50*chill[i]+a_t50))^(b_warm_beta*warm[i]+b_chill_beta*chill[i]+a_beta))); //log logistic equation where paramenters beta, t50 and d are represent as chilling sub mod
     
     }
 

 model {
  // priors
  //t50 ~ uniform(0, 100); 
  //beta ~ normal(0, 50); 
  //d ~ uniform(0, 30); 
  a_beta~normal(0,10);
  a_t50~uniform(0,20);
  a_d ~ uniform(0,10);
  b_warm_beta ~normal(0,5);
  b_warm_t50 ~normal(0,5);
  b_warm_d ~normal(0,5);
  
  b_chill_beta ~normal(0,1);
  b_chill_t50 ~normal(0,1);
  b_chill_d ~normal(0,1);
  
  sigma ~ normal(0, 10);
 
 
  // likelihood
  Y ~ normal(y_hat, sigma);//it seems like the model is having trouble with the normal distrubtion here. Should y actually be a log logtistic?
                            //maybe this would be helpful https://groups.google.com/forum/#!topic/stan-users/OexQuNGnsko
                            //or this https://mc-stan.org/docs/2_18/functions-reference/logistic-distribution.html
}
