bb_one <- function(N_subj = 100,T_obs  = 10, BB = 500, Nv = 200, alpha  = 5, doses  = c(6, 8, 10), seed = seed){
  set.seed(seed)
  true_vals  = sapply(doses, true_mu_mc)
  ############################################################
  ## Simulation parameters
  ############################################################
  N <- N_subj * T_obs
  mux <- 0.2; varx <- 0.1
  muu <- 1.0; varu <- 0.6
  muz <- 0; varz <- 0.1


  stick.breaking <- function(alpha, Nv){
    u <- rbeta(Nv, 1, alpha)
    w <- c(1, cumprod(1 - u[-Nv]))
    u * w
  }

  ############################################################
  ## Generate data
  ############################################################
  id   <- rep(1:N_subj, each = T_obs)
  time <- rep(1:T_obs, times = N_subj)

  x <- rnorm(N, mux, sqrt(varx))
  u <- rnorm(N, muu, sqrt(varu))
  z <- rep(rnorm(N_subj, muz, sqrt(varz)), each = T_obs)

  d <- 5 + 4 * x + 2 * u + 1 * z + rnorm(N)
  y <- 20 * exp(d + x - 0.25 * u + 0.5 * z) * exp(rnorm(N))

  sim.data <- data.frame(id = id, time = time, x = x, u = u, z = z, d = d, logy = log(y))


  cov_bb  <- matrix(NA, BB, length(doses))
  wor_bb  <- matrix(NA, BB, length(doses))

  for (b in 1:BB){

    w_subj <- as.numeric(rdirichlet(1, rep(1, N_subj)))
    w <- rep(w_subj, each = T_obs)

    ##########################################################
    ## GPS model
    ##########################################################
    ps_bb <- geeglm( d ~ x + u, id = id, data = sim.data, weights = w, family = gaussian, corstr = "independence")

    sigma_bb <- as.numeric(sqrt(summary(ps_bb)$geese$scale[1]))
    gps_bb <- dnorm(sim.data$d, ps_bb$fitted.values, sigma_bb)

    ##########################################################
    ## Outcome regression (GPS as covariate)
    ##########################################################
    or_bb <- geeglm(logy ~ ns(gps_bb, 3) + d, id = id, data = sim.data, weights = w, family = gaussian, corstr = "independence")


    ##########################################################
    ## Stabilized-weight GPS (WOR)
    ##########################################################
    d_density <- npudens(~ d, data = sim.data)
    f_d <- predict(d_density, newdata = sim.data)

    sw <- f_d / gps_bb

    or_wbb <- geeglm(logy ~ ns(d, 3), id = id, data = sim.data, weights = w * sw, family = gaussian, corstr = "independence")

    ## Store posterior means at fixed doses
    for (k in seq_along(doses)){
      cov_bb[b, k] <- mean(predict(or_bb, newdata = data.frame(d = rep(doses[k],nrow(sim.data)))))
      wor_bb[b, k] <- mean(predict(or_wbb, newdata = data.frame(d = rep(doses[k],nrow(sim.data)))))
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

