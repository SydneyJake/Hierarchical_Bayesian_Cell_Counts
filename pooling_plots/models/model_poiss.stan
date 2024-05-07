data{
  int                            N;
  array[N] int<lower=0>          y;
  real<lower=0>                  tau;
}
parameters{
  real      theta;
  vector[N] gamma_tilde;
}
transformed parameters{
  vector[N]          gamma;
  gamma = theta + gamma_tilde * tau;
}
model{
  gamma_tilde ~ std_normal();
  y ~ poisson_log(gamma);
}
