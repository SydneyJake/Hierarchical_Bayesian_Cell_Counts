functions{
  int num_zeros(array[] int y){
    int sum = 0;
    for(n in 1:size(y)){
      if(y[n] == 0){
        sum += 1;
      }
    }
    return sum;
  }
}
data{
  int                                A;          // Number of animals
  int                                R;          // Number of regions
  int                                G;          // Number of groups
  array[A]     int<lower=1, upper=G> group_idx;  // Group membership
  array[A, R]  int<lower=0>          y;          // Observations
}
transformed data{
  int<lower=0> N = A*R;
  array[N] int y_vec;
  int<lower=0> N_zero;

  // Flatten the array
  for(a in 1:A){
    for(r in 1:R){
      y_vec[(a - 1) * R + r] = y[a,r];
    }
  }

  N_zero = num_zeros(y_vec);
  array[N - N_zero] int<lower=1> y_nonzero_idx;
  array[N_zero]     int<lower=1> y_zero_idx;

  int zero_count    = 0;
  int nonzero_count = 0;

  for(n in 1:N){
    if(y_vec[n] == 0){
      zero_count += 1;
      y_zero_idx[zero_count] = n;
    }
    else{
      nonzero_count += 1;
      y_nonzero_idx[nonzero_count] = n;
    }
  }
}
parameters{
  array[G] vector<lower=0>[R] tau;
  array[G] vector[R]          theta;
  array[A] vector[R]          gamma_raw;
  real<lower=0, upper=1>      pi;
}
transformed parameters{
  array[A] vector[R] gamma;
  array[N] real      gamma_vec;

  for(a in 1:A){
    gamma[a] = theta[group_idx[a]] + tau[group_idx[a]] .* gamma_raw[a];
    for(r in 1:R){
      gamma_vec[(a - 1) * R + r] = gamma[a,r];
    }
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
  }

  // inflation
  pi ~ beta(1, 5);

  // observed data
  for(n in 1:size(y_zero_idx)){
    target += log_sum_exp(log(pi), log1m(pi) + poisson_log_lpmf(0 | gamma_vec[y_zero_idx[n]]));
  }
  target += nonzero_count * log1m(pi);
  target += poisson_log_lpmf(y_vec[y_nonzero_idx] | gamma_vec[y_nonzero_idx]);
}
generated quantities{
  array[A, R] int y_rep; // posterior predictions

  for(i in 1:A){
    for(j in 1:R){
      int ind = bernoulli_rng(pi); // flip a coin with probability P(Heads) = pi.
      y_rep[i,j] = (1-ind) * poisson_log_rng(gamma[i,j]); // If heads y_rep[i] = 0, else draw from a Poisson.
    }
  }
}
