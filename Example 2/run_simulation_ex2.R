### Load packages
library(splines)
library(np)
library(MCMCpack)
library(geepack)
library(rstan)
library(rstanarm)
rstan_options(auto_write = TRUE)

############################################################
## Load source functions
## Each file implements one simulation replication
############################################################
source("Bayes_MSM_ex2.R")
source("BB_ex2.R")
source("DP_ex2.R")

############################################################
## True marginal dose–response curve (MC)
############################################################
set.seed(1213)
n_mc <- 100000

## Draw covariates from the true DGP
mux <- 0.2; varx <- 0.1
muu <- 1.0; varu <- 0.6
muz <- 0.2; varz <- 0.1

z_mc <- rep(rnorm(n_mc / 10, muz, sqrt(varz)), each = 10)
x_mc <- rnorm(n_mc, mux , sqrt(varx)) + 0.25 * z_mc
u_mc <- rnorm(n_mc, muu, sqrt(varu))
## True dose–response function via Monte Carlo
true_mu_mc <- function(d){
  mean(exp(1 + 0.2 * d + 0.005 * x_mc/100 -0.002 * u_mc/100 + 0.1 * z_mc))
}

############################################################
## Run 1000 simulation replications (parallelized)
############################################################

library(parallel)
options(mc.cores = 20)
set.seed(213)
seeds <- sample(c(1:100000000),1000)


BB_ex2_res<-mclapply(seeds,function(x) bb_one_pois(seed = x))
DP_ex2_res<-mclapply(seeds,function(x) dp_one_pois(seed = x))
MSM_ex2_res<-mclapply(seeds,function(x) msm_one_replication_pois(seed = x))
