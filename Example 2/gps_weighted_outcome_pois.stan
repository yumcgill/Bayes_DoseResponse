data {
  int<lower=1> N;                      // observations
  int<lower=1> N_subj;                 // subjects
  int<lower=1> Kd;                     // spline basis dimension
  int<lower=1, upper=N_subj> id[N];    // subject id
  matrix[N, Kd] B_d;                   // spline basis for dose
  int<lower=0> y[N];                   // outcome (counts)
  vector[N] w;                         // stabilized GPS weights
}

parameters {
  real alpha;
  vector[Kd] beta_d;
  vector[N_subj] b_y;                  // random intercepts
  real<lower=0> sigma_by;              // SD of random intercepts
}

model {
  // Priors
  alpha     ~ normal(0, 10);
  beta_d   ~ normal(0, 10);
  b_y      ~ normal(0, sigma_by);
  sigma_by ~ normal(0, 5);

  // Weighted Poisson log-likelihood (log link)
  for (n in 1:N) {
    target += w[n] * poisson_log_lpmf(
      y[n] |
      alpha + B_d[n] * beta_d + b_y[id[n]]
    );
  }
}


