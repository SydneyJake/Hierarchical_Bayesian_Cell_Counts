data{
  int                            N;
  array[N] int<lower=0>          y;
  real<lower=0>                  tau;
}
parameters{
  real               theta;
  vector[N]          gamma_tilde;
  vector<lower=0>[N] lambda;
}
transformed parameters{
  vector[N] gamma;

  gamma = theta + lambda .* gamma_tilde * tau;
}
model{
  gamma_tilde  ~ std_normal();
  lambda       ~ normal(0, 1);
  y            ~ poisson_log(gamma);
}
