bb_one_pois <- function(N_subj = 100,T_obs  = 10, BB = 500, Nv = 200, doses  = c(6, 8, 10), seed = seed){
  set.seed(seed)
  true_vals  = sapply(doses, true_mu_mc)
  ############################################################
  ## Simulation parameters
  ############################################################
  N <- N_subj * T_obs
  mux <- 0.2; varx <- 0.1
  muu <- 1.0; varu <- 0.6
  muz <- 0.2; varz <- 0.1


  ############################################################
  ## Generate data
  ############################################################
  id   <- rep(1:N_subj, each = T_obs)
  time <- rep(1:T_obs, times = N_subj)

  z <- rep(rnorm(N_subj, muz, sqrt(varz)), each = T_obs)
  x <- rnorm(N, mux, sqrt(varx)) + 0.25 * z
  u <- rnorm(N, muu, sqrt(varu))
  d <- 5 + 4 * x + 2 * u + rnorm(N)

  y <- rpois(1000, exp(1 + 0.2 * d + 0.005 * x/100 -0.002 * u/100 + 0.1 * z))

  sim.data <- data.frame(id = id, time = time, x = x, u = u, z = z, d = d, y = y)


  cov_bb  <- matrix(NA, BB, length(doses))
  wor_bb  <- matrix(NA, BB, length(doses))

  for (b in 1:BB){

    w_subj <- as.numeric(rdirichlet(1, rep(1, N_subj)))
    w <- rep(w_subj, each = T_obs)

    ##########################################################
    ## GPS model
    ##########################################################
    ps_bb <- geeglm(d ~ x + u, id = id, data = sim.data, weights = w, family = gaussian, corstr = "independence")

    sigma_bb <- as.numeric(sqrt(summary(ps_bb)$geese$scale[1]))
    gps_bb <- dnorm(sim.data$d, ps_bb$fitted.values, sigma_bb)

    ##########################################################
    ## Outcome regression (GPS as covariate)
    ##########################################################
    or_bb <- geeglm(y ~ ns(gps_bb, 3) + d, id = id, data = sim.data, weights = w, family = poisson(link="log"), corstr = "independence")


    ##########################################################
    ## Stabilized-weight GPS (WOR)
    ##########################################################
    d_density <- npudens(~ d, data = sim.data)
    f_d <- predict(d_density, newdata = sim.data)

    sw <- f_d / gps_bb

    or_wbb <- geeglm(y ~ ns(d, 3), id = id, data = sim.data, weights = w * sw, family = poisson(link="log"), corstr = "independence")

    ## Store posterior means at fixed doses
    for (k in seq_along(doses)){
      cov_bb[b, k] <- mean(predict(or_bb, newdata = data.frame(d = rep(doses[k],nrow(sim.data))), type = "response"))
      wor_bb[b, k] <- mean(predict(or_wbb, newdata = data.frame(d = rep(doses[k],nrow(sim.data))), type = "response"))
    }
  }

  ############################################################
  ## Coverage indicators
  ############################################################
  cov_ci  <- apply(cov_bb, 2, quantile, probs = c(0.025, 0.975))
  wor_ci  <- apply(wor_bb, 2, quantile, probs = c(0.025, 0.975))

  cov_cover <- as.numeric(true_vals >= cov_ci[1, ] & true_vals <= cov_ci[2, ])

  wor_cover <- as.numeric(true_vals >= wor_ci[1, ] & true_vals <= wor_ci[2, ])

  ############################################################
  ## Return posterior summaries
  ############################################################
  list(cov_bb_mean = colMeans(cov_bb), wor_bb_mean = colMeans(wor_bb), cov_cover   = cov_cover, wor_cover  = wor_cover)
}
