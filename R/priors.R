#########################################################
# priors.R                                              #
#                                                       #
# This file contains functions for setting up the       #
# list of priors for the MCMC calibration.              #
#########################################################

#########################################################
# set_prior_params(parnames):                           #
#   This function sets up a data frame of prior         #
#   parameters with more "intuitive" parameters such as #
#   confidence interval bounds, then we convert those   #
#   parameters to those required by the density         #
#   functions when necessary.                           #
#########################################################

set_prior_params <- function(parnames) {

  # initialize prior parameter data frame
  prior_df <- data.frame(name=parnames, type=character(length(parnames)), lower=numeric(length(parnames)), upper=numeric(length(parnames)), stringsAsFactors=FALSE)
  
  # loop over data frame rows and fill in values
  for (i in 1:nrow(prior_df)) {
    name <- prior_df[i, 'name']
    if (name == 'psi1') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0.0001
      prior_df[i, 'upper'] <- 0.15
    } else if (name == 'psi2') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0
      prior_df[i, 'upper'] <- 100
    } else if (name == 'psi3') {
      prior_df[i, 'type'] <- 'normal'
      prior_df[i, 'lower'] <- 6.9
      prior_df[i, 'upper'] <- 14.4
    } else if (name == 'P0') {
      prior_df[i, 'type'] <- 'normal'
      prior_df[i, 'lower'] <- 0.3
      prior_df[i, 'upper'] <- 0.9
    } else if (name == 'lambda') {
      prior_df[i, 'type'] <- 'normal'
      prior_df[i, 'lower'] <- 0.6
      prior_df[i, 'upper'] <- 0.8
    } else if (name == 's') {
      prior_df[i, 'type'] <- 'normal'
      prior_df[i, 'lower'] <- 0.18
      prior_df[i, 'upper'] <- 0.26
    } else if (name == 'delta') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0.01
      prior_df[i, 'upper'] <- 0.14
    } else if (name == 'alpha') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0.0007
      prior_df[i, 'upper'] <- 0.0212
    } else if (name == 'As') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 5.3
      prior_df[i, 'upper'] <- 16.11
    } else if (name == 'pi') {
      prior_df[i, 'type'] <- 'normal'
      prior_df[i, 'lower'] <- 0.51
      prior_df[i, 'upper'] <- 0.67
    } else if (name == 'A0') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0
      prior_df[i, 'upper'] <- 3
    } else if (name == 'rho2') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0
      prior_df[i, 'upper'] <- 0.5
    } else if (name == 'rho3') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0
      prior_df[i, 'upper'] <- 0.5
    } else if (name == 'tau2') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 1700
      prior_df[i, 'upper'] <- 2100
    } else if (name == 'tau3') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 1700
      prior_df[i, 'upper'] <- 2100
    } else if (name == 'tau4') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 2010
      prior_df[i, 'upper'] <- 2500
    } else if (name == 'kappa') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0.005
      prior_df[i, 'upper'] <- 0.2
    } else if (grepl('rho', name)) {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0.4
      prior_df[i, 'upper'] <- 1
    } else if (grepl('sigma', name)) {
      prior_df[i, 'type'] <- 'log-normal'
      prior_df[i, 'lower'] <- 0
      prior_df[i, 'upper'] <- Inf
    } else if (grepl('eps', name)) {
      prior_df[i, 'type'] <- 'log-normal'
      prior_df[i, 'lower'] <- 0
      prior_df[i, 'upper'] <- Inf
    }
    
  }
  
  prior_df
}

#########################################################
# prior_list(parnames):                                 #
#   This function creates the list of prior             #
#     distribution information, including density       #
#     functions and parameters for those functions.     #
#     First it calls set_prior_params() to construct a  #
#     data frame of more intuitive prior information,   #
#     and then it converts that information to the      #
#     prior distribution function parameters.           #
#                                                       #
#   If the prior is uniform, the lower and upper data   #
#    frame parameters become the min and max            #
#    parameters. If the prior is normal, they're the    #
#    bounds of the 95% confidence interval. If the      #
#    prior is half-normal, they're irrelevant, as we    #
#    just use the standard half-Cauchy.                 #
#########################################################
create_prior_list <- function(parnames) {
  prior_df <- set_prior_params(parnames)
  priors <- vector('list', length(parnames))
  names(priors) <- parnames
  
  for (i in 1:nrow(prior_df)) {
    name <- prior_df[i, 'name']
    priors[[name]] <- list(type=prior_df[i, 'type'])
    if (priors[[name]][['type']] == 'uniform') {
      priors[[name]][['dens.fun']] <- 'dunif'
      priors[[name]][['quant.fun']] <- 'qunif'
      priors[[name]][['min']] <- prior_df[i, 'lower']
      priors[[name]][['max']] <- prior_df[i, 'upper']
    } else if (priors[[name]][['type']] == 'normal') {
      priors[[name]][['dens.fun']] <- 'dnorm'
      priors[[name]][['quant.fun']] <- 'qnorm'
      priors[[name]][['mean']] <- mean(c(prior_df[i, 'lower',], prior_df[i, 'upper']))
      priors[[name]][['sd']] <- (prior_df[i, 'upper'] - prior_df[i, 'lower'])/(qnorm(0.975) - qnorm(0.025))
    } else if (priors[[name]][['type']] == 'log-normal') {
      priors[[name]][['dens.fun']] <- 'dlnorm'
      priors[[name]][['quant.fun']] <- 'qlnorm'
      priors[[name]][['meanlog']] <- -1
      priors[[name]][['sdlog']] <- 1
    }
  }
  
  priors
}
