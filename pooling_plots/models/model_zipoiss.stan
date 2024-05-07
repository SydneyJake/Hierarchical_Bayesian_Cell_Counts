data{
  int                            N;
  array[N] int<lower=0>          y;
  real<lower=0>                  tau;
}
parameters{
  real                   theta;
  vector[N]              gamma_tilde;
  real<lower=0, upper=1> pi;
}
transformed parameters{
  vector[N] gamma;
  gamma = theta + gamma_tilde * tau;
}
model{
  pi ~ beta(1,5);
  gamma_tilde ~ std_normal();

  for(i in 1:N){
    if(y[i] == 0){
      target += log_sum_exp(log(pi),
                            log1m(pi) + poisson_log_lpmf(0 | gamma[i]));
    }
    else{
      target += log1m(pi) + poisson_log_lpmf(y[i] | gamma[i]);
    }
  }
}
