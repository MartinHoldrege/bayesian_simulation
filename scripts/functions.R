# Martin Holdrege

# Script started 4/23/20

# functions used in WILD 6900 final proj


library(dplyr)


# simulate data -----------------------------------------------------------

sim_dat <- function(K = 3, J = 6, mu_beta1 = 0.2, sigma_beta1 = 0.3,
                    sigma_eps = 3) {
  # args:
  #     K replicates per treatment
  #     J  number of years
  #    mu_beta1--mean slope (hyperparameter)
  #   sigma_beta1  sd of slope (hyperperameter)
  #   sigma_eps sd of error term
  # returns:
  #   dataframe


  # treatment levels
  trmts <- c(5, 10, 15, 20)

  mu_beta0 <- 4

  sigma_beta0 <- 0.4


  beta0 <- rnorm(J, mu_beta0, sigma_beta0)
  beta1 <- rnorm(J, mu_beta1, sigma_beta1)

  df1 <- expand.grid(year = 1:J,
                     replicate = 1:J,
                     trmt = trmts)

  out <- df1 %>%
    mutate(y = beta0[year] + beta1[year]*trmt + rnorm(nrow(df1), 0, sigma_eps),
           year = as.factor(year)) %>%
    as_tibble()

  out
}

sim_dat()


# priors ------------------------------------------------------------------

prior_beta <- function(x){
  stopifnot(length(x) == 1) # i only want it to return scalar prob
  prob <- dnorm(x, 0, 10, log = TRUE)
  return(prob)
}

prior_sigma <- function(x) {
  stopifnot(length(x) == 1) # i only want it to return scalar prob
  prob <- dunif(x, min = 0, max = 20, log = TRUE)
  return(prob) # log of probability
}



