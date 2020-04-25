# Martin Holdrege

# Script started 4/23/20

# functions used in WILD 6900 final proj


library(dplyr)


# simulate data -----------------------------------------------------------

sim_dat <- function(K = 3, I = 6, mu_beta1 = 0.2, sigma_beta1 = 0.3,
                    sigma_eps = 3) {
  # args:
  #     K replicates per treatment
  #     i  number of years
  #    mu_beta1--mean slope (hyperparameter)
  #   sigma_beta1  sd of slope (hyperperameter)
  #   sigma_eps sd of error term
  # returns:
  #   dataframe


  # treatment levels
  trmts <- c(5, 10, 15, 20)

  mu_beta0 <- 4

  sigma_beta0 <- 0.4


  beta0 <- rnorm(I, mu_beta0, sigma_beta0)
  beta1 <- rnorm(I, mu_beta1, sigma_beta1)

  df1 <- expand.grid(year = 1:I,
                     replicate = 1:I,
                     trmt = trmts)

  out <- df1 %>%
    mutate(y = beta0[year] + beta1[year]*trmt + rnorm(nrow(df1), 0, sigma_eps),
           year = as.factor(year)) %>%
    as_tibble()

  out
}

# sim_dat()


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


# likelihood --------------------------------------------------------------


calc_like <- function(y, year, x, beta0, beta1, mu_beta0, mu_beta1, sigma_beta0,
                      sigma_beta1, sigma_eps) {
  # args:
  #   y--respone
  #   year--same length as y
  #   x--same length as y
  #   beta0 and beta1--same length as unique(year)
  # mu_beta0, mu_beta1, sigma_beta0, sigma_beta1--hyper parameters
  # returns:
  #   log likelihood

  mu_y <- beta0[year] + beta1[year]*x

  ll1 <- sum(dnorm(y, mu_y, sigma_eps, log = TRUE))

  # likelihood of intercepts given hyper parameter
  ll2 <- sum(dnorm(beta0, mu_beta0, sigma_beta0, log = TRUE))



  # likelihood of slopes given hyper parameter
  ll3 <- sum(dnorm(beta1, mu_beta1, sigma_beta1, log = TRUE))


  total_ll <- ll1 + ll2 + ll3

  total_ll
}


# likelihood*prior --------------------------------------------------------

calc_joint <- function(y, year, x, beta0, beta1, mu_beta0, mu_beta1, sigma_beta0,
                       sigma_beta1, sigma_eps) {

  ll <- calc_like(y = y, year = year , x = x, beta0 = beta0, beta1 = beta1,
                  mu_beta0 = mu_beta0, mu_beta1 = mu_beta1,
                  sigma_beta0 = sigma_beta0,
                  sigma_beta1 = sigma_beta1, sigma_eps = sigma_eps)

  joint <- ll + prior_beta(mu_beta0) + prior_beta(mu_beta1) +
    prior_sigma(sigma_beta0) + prior_sigma(sigma_beta1) + prior_sigma(sigma_eps)
  joint
}


# initial values ----------------------------------------------------------

jags_inits <- function(){list(mu_beta0 = rnorm(1),
                              mu_beta1 = rnorm(1),
                              sigma_beta0 = runif(1),
                              sigma_beta1 = runif(1),
                              beta0 = rnorm(6),
                              beta1 = rnorm(6),
                              sigma_eps = runif(1))}

# df <- sim_dat()
# inits <- jags_inits()

# create df for mcmc ---------------------------------------------------------------

# creates empty (except first row) df for mcmc
create_mcmc_df <- function(inits, niter) {
  # inits--list of initial values
  # niter--number of iterations

  n_groups <- length(inits$beta0)

  mcmc_df <- tibble(iteration = 1:niter,
                    mu_beta0 = NA_real_,
                    mu_beta1 = NA_real_,
                    sigma_beta0 = NA_real_,
                    sigma_beta1 = NA_real_,
                    sigma_eps = NA_real_,
                    mu_beta0_accept = 0,
                    mu_beta1_accept = 0,
                    sigma_beta0_accept = 0,
                    sigma_beta1_accept = 0,
                    sigma_eps_accept = 0)

  # add beta columns
  for (i in 1:n_groups) {
    mcmc_df[[paste0("beta0", i)]] <- NA_real_
    mcmc_df[[paste0("beta1", i)]] <- NA_real_
    mcmc_df[[paste0("beta0", i, "_accept")]] <- 0
    mcmc_df[[paste0("beta1", i, "_accept")]] <- 0
    # initial values
    mcmc_df[[paste0("beta0", i)]][1] <- inits$beta0[i]
    mcmc_df[[paste0("beta1", i)]][1] <- inits$beta1[i]
  }


  # reorder cols so easier to read
  mcmc_df <- mcmc_df %>%
    select(iteration, matches("^mu.+\\d$"), matches("sigma"),
           matches("^beta.+\\d$"), matches("accept$"), everything())

  # set initial values

  # initial values for the rest of the columns
  single_inits <- names(inits)[!names(inits) %in% c("beta0", "beta1")]

  for (col in single_inits) {
    mcmc_df[1, ][[col]] <- inits[[col]]
  }

  mcmc_df
}

# mcmc_df <- create_mcmc_df(inits = jags_inits(), niter = 100)


# proposal values ---------------------------------------------------------

# propose a new values
proposal <- function(current, tune, is_sigma = FALSE) {

  if (current < 0 & is_sigma == TRUE) stop("can't provide negative current value")

  proposal <- rnorm(1, current, sd = tune)

  # recursion to get new value that is positive--when don't want to propose
  # negative sd's
  while (is_sigma & proposal < 0) {
    proposal <- rnorm(1, current, sd = tune)
  }

  proposal
}

if (FALSE){
  for (i in 1:10) {
    print(proposal(0, 1))
  }

  x <- vector("numeric", 10000)
  for (i in 1:10000) {
    x[i] <- proposal(0, 1, is_sigma = TRUE)
  }
  hist(x)
}

# propose new values and keep/reject--------------------------------------

# last not na value in a vector
last_value <- function(x) {
  x[max(which(!is.na(x)))]
}

# accept or reject proposal values
propose_check <- function(df, mcmc_df, i, col_name, tune, likelihood = NULL,
                          is_sigma) {
  # df--df of data
  # mcmc_df--df from create_mc_df()
  # i--row number
  # col_name--name in mcmc_df--only differet from target_name for beta vectors
  # tuning parameter
  # likelihood --current likelihood (if known), ie prior to proposal value

  # important columns in mcmc_df
  cols <- c('mu_beta0', 'mu_beta1', 'sigma_beta0', 'sigma_beta1', 'sigma_eps',
            'beta01', 'beta11', 'beta02', 'beta12', 'beta03', 'beta13', 'beta04',
            'beta14', 'beta05', 'beta15', 'beta06', 'beta16')
  stopifnot(
    cols %in% names(mcmc_df),
    # check that previous step mcmc filled in a value
    max(which(!is.na(mcmc_df[[col_name]]))) == i - 1
  )

  cols <- sort(cols)

  # current row (take last value because some will be from row i other i-1)
  l <- mcmc_df[cols] %>%
    summarise_all(.funs = last_value) %>%
    as.numeric()

  names(l) <- cols

  beta0_names <- paste0("beta0", 1:6)
  beta1_names <- paste0("beta1", 1:6)



  if (is.null(likelihood)) {
    joint_old <- calc_joint(
      y = df$y, year = df$year , x = df$trmt, beta0 = l[beta0_names],
      beta1 = l[beta1_names], mu_beta0 = l["mu_beta0"], mu_beta1 = l["mu_beta1"],
      sigma_beta0 = l["sigma_beta0"], sigma_beta1 = l["sigma_beta1"],
      sigma_eps = l["sigma_eps"])
  } else {
    joint_old <- likelihood
  }


  # update with proposal values
  l[[col_name]] <- proposal(l[[col_name]], tune, is_sigma = is_sigma)

  # candidate likelihood
  joint_cand <- calc_joint(
    y = df$y, year = df$year , x = df$trmt, beta0 = l[beta0_names],
    beta1 = l[beta1_names], mu_beta0 = l["mu_beta0"], mu_beta1 = l["mu_beta1"],
    sigma_beta0 = l["sigma_beta0"], sigma_beta1 = l["sigma_beta1"],
    sigma_eps = l["sigma_eps"])

  ## 1d: Acceptance probability
  r <- exp(joint_cand - joint_old) # joint likelihoods are on log scale
  R <- min(1, r)
  ## 1e: Decide whether to accept or not

  if(!is.na(R) && runif(1) <= R) {   # if accepted
    mcmc_df[[col_name]][i] <- l[[col_name]]
    mcmc_df[[paste0(col_name, "_accept")]][i] <-1
    likelihood <- joint_cand
  } else {
    mcmc_df[[col_name]][i] <-mcmc_df[[col_name]][i - 1]
    likelihood <- joint_old
  }

  out <- list(mcmc_df = mcmc_df, likelihood = likelihood)
  out
}


# mcmc loop ---------------------------------------------------------------


mcmc_loop <- function(df, inits, tune, niter = 100) {

  mcmc_df <- create_mcmc_df(inits, niter = niter)
  # at this point the likelihood is unknown
  mcmc_l <- list(mcmc_df = mcmc_df, likelihood = NULL)

  cols <- c('mu_beta0', 'mu_beta1', 'sigma_beta0', 'sigma_beta1', 'sigma_eps',
            'beta01', 'beta11', 'beta02', 'beta12', 'beta03', 'beta13', 'beta04',
            'beta14', 'beta05', 'beta15', 'beta06', 'beta16')

  for(i in 2:niter) {
    for (col in cols) {
      is_sigma <- stringr::str_detect(col, "sigma")

      mcmc_l <- propose_check(df = df, mcmc_df = mcmc_l$mcmc_df,
                              i = i, col_name = col, tune = tune,
                              likelihood = mcmc_l$likelihood,
                              is_sigma = is_sigma)
    }

  }

  out <- mcmc_l$mcmc_df
  out
}

# mcmc_loop(df = sim_dat(),
#           inits = jags_inits(),
#           tune = 0.5,
#           niter = 100)


# determine tuning parameter ----------------------------------------------


if (FALSE) {
  mcmc_df <- mcmc_loop(df = sim_dat(),
            inits = jags_inits(),
            tune = 0.5,
            niter = 1000)
  mcmc_df %>%
    select(matches("accept")) %>%
    unlist() %>%
    mean()

  # note beta0's have lower variability so acceptance rate is
  # quite low, they would need a separate tuning parameter.
  mcmc_df %>%
    select(matches("accept")) %>%
    summarise_all(mean) %>%
    unlist()
}

