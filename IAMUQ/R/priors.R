#########################################################
# priors.R                                              #
#                                                       #
# This file contains a function for setting up the      #
# list of priors for the MCMC calibration.              #
#########################################################

#########################################################
# prior_list(parnames):                                 #
#   This function creates the list of prior             #
#     distribution information, including density       #
#     functions and parameters for those functions,     #
#     based on a dataframe that is passed in which      #
#     contains the name of parameters and required      #
#     information about the distribution for each       #
#     parameter. It converts that information to the    #
#     prior distribution function parameters.           #
#                                                       #
#   If the prior is uniform, the lower and upper data   #
#    frame parameters become the min and max            #
#    parameters. If the prior is normal, they're the    #
#    bounds of the 95% confidence interval. If the      #
#    prior is half-normal, they're irrelevant, as we    #
#    just use the standard half-Cauchy.                 #
#########################################################
create_prior_list <- function(prior_df) {
  parnames <- prior_df[, 'name']
  priors <- vector('list', length(parnames))
  names(priors) <- parnames
  
  for (i in 1:nrow(prior_df)) {
    name <- prior_df[i, 'name']
    priors[[name]] <- list(type=prior_df[i, 'type'])
    if (priors[[name]][['type']] == 'uniform') {
      priors[[name]][['dens.fun']] <- 'dunif'
      priors[[name]][['quant.fun']] <- 'qunif'
      priors[[name]][['rand.fun']] <- 'runif'
      priors[[name]][['min']] <- prior_df[i, 'lower']
      priors[[name]][['max']] <- prior_df[i, 'upper']
    } else if (priors[[name]][['type']] == 'normal') {
      priors[[name]][['dens.fun']] <- 'dnorm'
      priors[[name]][['quant.fun']] <- 'qnorm'
      priors[[name]][['rand.fun']] <- 'rnorm'
      priors[[name]][['mean']] <- mean(c(prior_df[i, 'lower'], prior_df[i, 'upper']))
      priors[[name]][['sd']] <- (prior_df[i, 'upper'] - prior_df[i, 'lower'])/(qnorm(0.975) - qnorm(0.025))
    } else if (priors[[name]][['type']] == 'log-normal') {
      priors[[name]][['dens.fun']] <- 'dlnorm'
      priors[[name]][['quant.fun']] <- 'qlnorm'
      priors[[name]][['rand.fun']] <- 'rlnorm'
      priors[[name]][['meanlog']] <- -1
      priors[[name]][['sdlog']] <- 1
    } else if (priors[[name]][['type']] == 'truncnorm') {
      priors[[name]][['dens.fun']] <- 'dtnorm'
      priors[[name]][['quant.fun']] <- 'qtnorm'
      priors[[name]][['rand.fun']] <- 'rtnorm'
      priors[[name]][['lower']] <- prior_df[i, 'lower']
      priors[[name]][['upper']] <- prior_df[i, 'upper']
      priors[[name]][['mean']] <- mean(c(prior_df[i, 'lower'], prior_df[i, 'upper']))
      priors[[name]][['sd']] <- (prior_df[i, 'upper'] - prior_df[i, 'lower'])/(qnorm(0.975) - qnorm(0.025))
    } else if (priors[[name]][['type']] == 'beta') {
      priors[[name]][['dens.fun']] <- 'dbeta'
      priors[[name]][['quant.fun']] <- 'qbeta'
      priors[[name]][['rand.fun']] <- 'rbeta'
      mu <- mean(prior_df[i, 'upper'], prior_df[i, 'lower'])
      sigma.sq <- (prior_df[i, 'upper'] - prior_df[i, 'lower'])^2 / 16
      priors[[name]][['shape1']] <- ((1 - mu) / sigma.sq - 1 / mu) * mu^2
      priors[[name]][['shape2']] <- priors[[name]][['shape1']] * (1 / mu - 1)
    }
  }
  
  priors
}
