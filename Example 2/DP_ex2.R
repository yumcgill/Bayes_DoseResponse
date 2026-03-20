dp_one_pois <- function(N_subj = 100,T_obs  = 10, BB = 1000, Nv = 100, alpha  = 5, doses  = c(6, 8, 10) , seed = seed){
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
  y <- rpois(1000, exp(1 + 0.2 * d + 0.005 * x/100 -0.002 * u/100 + 0.1 * z))
  
  sim.data <- data.frame(id = id, time = time, x = x, u = u, z = z, d = d, y = y)
  
  ############################################################
  ## Initial GPS + outcome model (base measure G0)
  ############################################################
  ps_gee <- geeglm(d ~ x + u, id = id, data = sim.data, family = gaussian, corstr = "independence")
  
  sigma_hat <- as.numeric(sqrt(summary(ps_gee)$geese$scale[1]))
  
  gps_hat   <- dnorm(sim.data$d, ps_gee$fitted.values, sigma_hat)
  
  or_init <- geeglm(y ~ ns(gps_hat, 3) + d, id = id, data = sim.data, family = poisson(link="log"), corstr = "independence")
  
  sim.data$pred <- predict(or_init)
  
  ############################################################
  ## DP posterior sampling
  ############################################################
  ids <- unique(sim.data$id)
  split_data <- split(sim.data, sim.data$id)
  
  cov_dp  <- matrix(NA, BB, length(doses))
  wor_dp  <- matrix(NA, BB, length(doses))
  
  for (b in 1:BB){
    
    subj_idx <- sample(ids, Nv, replace = TRUE)
    w_dp <- stick.breaking(alpha + length(ids), Nv)
    
    newdata <- NULL
    for (j in seq_along(subj_idx)){
      tmp <- split_data[[as.character(subj_idx[j])]]
      tmp$w_dp <- w_dp[j]
      
      if (runif(1) < alpha / (alpha + length(ids))){
        tmp$y <- rpois(nrow(tmp), tmp$pred)
      }
      newdata <- rbind(newdata, tmp)
    }
    
    ## GPS re-estimation
    ps_bb <- geeglm(d ~ x + u, id = id, data = newdata, weights = newdata$w_dp, family = gaussian, corstr = "independence")
    
    sigma_bb <- as.numeric(sqrt(summary(ps_bb)$geese$scale[1]))
    gps_bb <- dnorm(newdata$d, ps_bb$fitted.values, sigma_bb)
    
    ## COV-DP
    or_bb <- geeglm(y ~ ns(gps_bb, 3) + d, id = id, data = newdata, weights = w_dp, family = poisson(link="log"), corstr = "independence")
    
    ## WOR-DP
    d_density <- npudens(~ d, data = newdata)
    f_d <- predict(d_density, newdata = newdata)
    sw <- f_d / gps_bb
    
    or_wbb <- geeglm(y ~ ns(d, 3), id = id, data = newdata, weights = w_dp * sw, poisson(link="log"), corstr = "independence")
    
    
    ## Store posterior means at fixed doses
    for (k in seq_along(doses)){
      cov_dp[b, k] <- mean(predict(or_bb, newdata = data.frame(d = rep(doses[k],nrow(newdata)), gps_bb=dnorm(doses[k], ps_bb$fitted.values, sigma_bb)),type = "response"))
      wor_dp[b, k] <- mean(predict(or_wbb, newdata = data.frame(d = doses[k]), type = "response"))
    }
  }
  
  ############################################################
  ## Coverage indicators
  ############################################################
  cov_ci  <- apply(cov_dp, 2, quantile, probs = c(0.025, 0.975))
  wor_ci  <- apply(wor_dp, 2, quantile, probs = c(0.025, 0.975))
  
  cov_cover <- as.numeric(true_vals >= cov_ci[1, ] & true_vals <= cov_ci[2, ])
  
  wor_cover <- as.numeric(true_vals >= wor_ci[1, ] & true_vals <= wor_ci[2, ])
  
  ############################################################
  ## Return posterior summaries
  ############################################################
  list(cov_dp_mean = colMeans(cov_dp), wor_dp_mean = colMeans(wor_dp), cov_cover   = cov_cover, wor_cover  = wor_cover)
}
