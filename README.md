# Fitting baysian models to simulated data
Final project for WILD 6900 Applied Bayesian Analysis of Ecological Data

## Description of files--

### `scripts/functions.R`

Contains the functions used to build the metropolis MCMC sampler. Also code to set up the JAGS models.

### `scripts/jags_simulations.R`

Simulate 500 data sets, and fit correct hierarchical model and incorrect simple linear model. 

### `scripts/HoldregeMartin_final_proj.rmd`  

Class project report that pulls in results from `jags_simulations.R`. Knitted pdf of this doc can also be found here: `scripts/HoldregeMartin_final_proj.pdf`. 
