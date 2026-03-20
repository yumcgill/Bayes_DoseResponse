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
source("Bayes_MSM_ex1.R")
source("BB_ex1.R")
source("DP_ex1.R")

############################################################
## True marginal dose–response curve (MC)
############################################################
set.seed(1213)
n_mc <- 100000

## Draw covariates from the true DGP
mux <- 0.2; varx <- 0.1
muu <- 1.0; varu <- 0.6
muz <- 0; varz <- 0.1
x_mc <- rnorm(n_mc, mux, sqrt(varx))
u_mc <- rnorm(n_mc, muu, sqrt(varu))
z_mc <- rep(rnorm(n_mc / 10, muz, sqrt(varz)), each = 10) 

## True dose–response function via Monte Carlo
true_mu_mc <- function(d){
  mean(log(20)  + (d + x_mc - 0.25 * u_mc + 0.5 * z_mc))
}

############################################################
## Run 1000 simulation replications (parallelized)
############################################################

library(parallel)
options(mc.cores = 20)
set.seed(213)
seeds <- sample(c(1:100000000),1000)


BB_ex1_res<-mclapply(seeds,function(x) bb_one(seed = x))
DP_ex1_res<-mclapply(seeds,function(x) dp_one(seed = x))
MSM_ex1_res<-mclapply(seeds,function(x) msm_one_replication(seed = x))
