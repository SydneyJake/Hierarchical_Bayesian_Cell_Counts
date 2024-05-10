data{
  int R;             // Number of regions
  int N;
  int G;
  array[N] int<lower=1, upper=G> group_idx;
  array[N] int<lower=1, upper=R> region_idx;
  vector[N] E;
  array[N] int<lower=0> y; // Observations
}
parameters{
  array[G] vector<lower=0>[R] tau;
  array[G] vector[R] theta;
  vector[N]          gamma_raw;
}
transformed parameters{
  vector[N] gamma;

  for(i in 1:N){
    gamma[i] = theta[group_idx[i], region_idx[i]] + gamma_raw[i] * tau[group_idx[i], region_idx[i]];
  }
}
model{
  for(i in 1:G){
    tau[i]   ~ normal(0, log(1.05));
    theta[i] ~ normal(5,2);
  }

  gamma_raw ~ std_normal();

  y ~ poisson_log(E + gamma);
}
generated quantities{
  array[N] int y_rep;

  for(i in 1:N){
      y_rep[i] = poisson_log_rng(E[i] + gamma[i]); // posterior predictions
  }
}
