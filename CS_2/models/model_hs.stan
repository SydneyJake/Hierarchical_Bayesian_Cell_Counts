data{
  int                                A;          // Number of animals
  int                                R;          // Number of regions
  int                                G;          // Number of groups
  array[A]     int<lower=1, upper=G> group_idx;  // Group membership
  array[A, R]  int<lower=0>          y;          // Observations
}
parameters{
  array[G] vector<lower=0>[R] tau;
  array[G] vector[R]          theta;
  array[A] vector[R]          gamma_raw;
  array[A] vector<lower=0>[R] kappa;
}
transformed parameters{
  array[A] vector[R] gamma;

  for(i in 1:A){
    gamma[i]  = theta[group_idx[i]] + tau[group_idx[i]] .* gamma_raw[i] .* kappa[i];
  }
}
model{

  // Population parameters
  for(i in 1:G){
    tau[i]   ~ normal(0, log(1.05));
    theta[i] ~ normal(5, 2);
  }

  // Random effects
  for(i in 1:A){
    gamma_raw[i] ~ std_normal();
    kappa[i]    ~ std_normal();
  }

  // Observed Data
  for(i in 1:A){
    y[i] ~ poisson_log(gamma[i]);
  }
}
generated quantities{
  array[A, R] int y_rep;

  for(i in 1:A){
    for(j in 1:R){
      y_rep[i, j] = poisson_log_rng(gamma[i,j]); // Posterior predictions
    }
  }

}
