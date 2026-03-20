msm_one_replication <- function(seed,
                                N_subj = 100,
                                T_obs  = 10,
                                doses  = c(6, 8, 10)) {
  start.time <- Sys.time()
  set.seed(seed)
  true_vals  = sapply(doses, true_mu_mc)
  ############################################################
  ## Simulation
  ############################################################
  N <- N_subj * T_obs
  id <- rep(1:N_subj, each = T_obs)

  x <- rnorm(N, 0.2, sqrt(0.1))
  u <- rnorm(N, 1.0, sqrt(0.6))
  z <- rep(rnorm(N_subj, 0, sqrt(0.1)), each = T_obs)

  d <- 5 + 4 * x + 2 * u + 1 * z + rnorm(N)
  y <- 20 * exp(d + x - 0.25 * u + 0.5 * z) * exp(rnorm(N))

  sim.data <- data.frame(id = id, x = x, u = u, z = z, d = d, logy = log(y))


  ##--------------------------------------------------
  ## GPS model (Bayesian parametric)
  ##--------------------------------------------------
  ## Bayesian GPS model
  gps_fit <- stan_glmer(d ~ x + u + (1 | id), data = sim.data, family = gaussian(), chains = 1, iter = 2000, warmup = 1000, refresh = 0)

  ## Posterior mean of conditional mean
  mu_d_hat <- posterior_linpred(gps_fit, transform = FALSE)
  mu_d_bar <- colMeans(mu_d_hat)

  ## Posterior mean of sigma
  sigma_d_hat <- mean(as.matrix(gps_fit, pars = "sigma"))

  ## Conditional GPS
  gps_cond <- dnorm(sim.data$d, mu_d_bar, sigma_d_hat)

  ## Marginal dose density
  d_dens <- npudens(~ d, data = sim.data)
  gps_marg <- predict(d_dens, newdata = data.frame(d = sim.data$d))

  ## Stabilized weights
  sim.data$sw <- gps_marg / gps_cond


  ## Weighted outcome GEE
  dose_B <- ns(sim.data$d, df = 3)
  stan_data_msm <- list(
    N = nrow(sim.data),
    N_subj = length(unique(sim.data$id)),
    Kd = ncol(dose_B),
    id = sim.data$id,
    B_d = dose_B,
    d = sim.data$d,
    logy = sim.data$logy,
    w = sim.data$sw
  )

  fit_msm <- stan(file = "gps_weighted_outcome.stan", data = stan_data_msm, chains = 1, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.95))



  ### dose-reponse curve


  d_vals <- doses

  ## MSM
  post_msm <- extract(fit_msm)

  B_grid <- ns(
    d_vals, df = 3,
    Boundary.knots = attr(dose_B, "Boundary.knots"),
    knots = attr(dose_B, "knots")
  )

  mu_msm <- sweep(post_msm$beta_d %*% t(B_grid),  # S x G
                  1,                              # sweep over rows
                  post_msm$alpha,                 # length S
                  "+"
  )

 
  ############################################################
  ## Coverage indicators
  ############################################################
  wor_ci  <- apply(mu_msm, 2, quantile, probs = c(0.025, 0.975))


  wor_cover <- as.numeric(true_vals >= wor_ci[1, ] & true_vals <= wor_ci[2, ])
  end.time <- Sys.time()
  tt <- end.time - start.time
  ############################################################
  ## Return posterior summaries
  ############################################################
  list(wor_msm_mean = colMeans(mu_msm), wor_cover  = wor_cover, runtime = tt)

}


