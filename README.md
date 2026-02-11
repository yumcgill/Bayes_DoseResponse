## Overview
This repository contains the R code supporting the paper by Yu Luo, Kuan Liu, Ramandeep Singh and Daniel J. Graham (2025) ["A longitudinal Bayesian framework for estimating causal dose-response relationships"](https://arxiv.org/abs/2505.20893). 

A nonparametric Bayesian framework is developed for estimating marginal longitudinal causal dose–response functions with repeated outcome measurements. This approach targets the average potential outcome at any fixed dose level and accommodates time-varying confounding through the generalized propensity score.  The proposed framework combines generalized estimating equations (GEE) with the Bayesian bootstrap (BB) and Dirichlet process (DP) priors to provide uncertainty quantification without specifying a full likelihood.

## Repository Structure

```
├── Example1/
│   ├── BB_ex1.R                        # Bayesian bootstrap methods for Example 1 in the paper
│   ├── DP_ex1.R                        # Dirichlet process methods for Example 1 in the paper
│   ├── Bayes_MSM_ex1.R                 # Parametric Bayesian MSM (Stan) for Example 1 in the paper
│   ├── gps_weighted_outcome.stan       # Stan code for a Gaussian Bayesian marginal structural model 
│                                       # with inverse GPS–weighted likelihood
│   ├── BB_ex1.R                        # Bayesian bootstrap methods for Example 1 in the paper
│   └── run_simulation_ex1.R            # Master simulation script for Example 1 in the paper
│
├── Example2/
│   ├── BB_ex2.R                        # Bayesian bootstrap methods for Example 2 in the paper
│   ├── DP_ex2.R                        # Dirichlet process methods for Example 2 in the paper
│   ├── gps_weighted_outcome_pois.stan  # Stan code for a Poisson Bayesian marginal structural model 
│                                       # with inverse GPS–weighted likelihood
│   ├── Bayes_MSM_ex2.R                 # Parametric Bayesian MSM (Stan) for Example 2 in the paper
│   └── run_simulation_ex2.R            # Master simulation script for Example 2 in the paper

```

## Running the Simulation Study

To reproduce the simulation results:
```
source("run_simulation_ex1.R")
source("run_simulation_ex2.R")
```

Parallel computation is enabled via:
```
options(mc.cores = 20)
```
Adjust the number of cores based on your system.

