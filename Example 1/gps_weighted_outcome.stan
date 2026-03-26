data {
  int<lower=1> N;                 // observations
  int<lower=1> N_subj;            // subjects
  int<lower=1> Kd;                // spline basis dimension
  int<lower=1,upper=N_subj> id[N];
  matrix[N, Kd] B_d;              // spline basis for dose
  vector[N] logy;                 // outcome
  vector[N] w;                    // stabilized GPS weights
}

parameters {
  real alpha;
  vector[Kd] beta_d;
  vector[N_subj] b_y;
  real<lower=0> sigma_y;
  real<lower=0> sigma_by;
}

model {
  // Priors
  alpha ~ normal(0, 100);
  beta_d ~ normal(0, 100);
  b_y ~ normal(0, sigma_by);
  sigma_y ~ normal(0, 10);
  sigma_by ~ normal(0, 10);

  // Weighted likelihood
  for (n in 1:N) {
    target += w[n] * normal_lpdf(
      logy[n] |
        alpha + B_d[n] * beta_d + b_y[id[n]],
      sigma_y
    );
  }
}
