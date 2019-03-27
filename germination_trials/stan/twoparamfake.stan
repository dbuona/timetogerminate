functions {
    real growth_factor_loglogistic(real t, real omega, real theta) {
        real pow_t_omega = t^omega;
        return pow_t_omega / (pow_t_omega + theta^omega);
    }
}

data{
 int n_data;
  int n_time;
  int chill;

vector<lower=0>[n_time] t_value;
}

parameters {
real beta_b;
real beta_a;

real t50_a;
real t50_b;
}

transformed parameters{
      real<lower=0> omega;
    real<lower=0> theta;

omega= beta_b*chill+beta_a;
theta=t50_b*chill+t50_a;
}

model{
  for(i in 1:n_time) 
target +=growth_factor_loglogistic(t_value[i], omega, theta);

}
