#' Run Markov Chain Monte Carlo (MCMC) to calibrate the IAM.
#'
#' \code{run_mcmc} uses adaptive MCMC (with the \code{adaptMCMC} package) to
#'  calibrate the IAM.
#'
#' When run in parallel, \code{run_mcmc} uses all of the available cores (or
#'  as many as are necessary to run the number of chains). The acceptance
#'  rate is based on the number of parameters, between 44% for one parameter
#'  and 23.4% for many parameters (from Roberts, Gelman, and Gilks (1997)).
#'
#' The covariance jump matrix is initialized using 10% of the prior
#'  distribution standard deviations. Adaptation begins after the maximum of
#'  500 parameters and 1% of the number of iterations, and continues until
#'  5e5 parameters (intended as burn-in).
#'
#' @param post Function call (not string) for the log-posterior density
#'  function.
#' @param parnames Character vector of parameter names (see the documentation
#'  for the desired likelihood function for which parameters are required).
#' @param residtype String for the type of residual structure (by default
#'  "iid" or "var"). The log-likelihood function associated with the
#'  structure should have the name log_lik_residtype.
#' @param prior_df Data frame of prior information (such as that produced by
#'  \code{\link{create_prior_list}}).
#' @param data_yrs Numeric vector of the years for which the data will be
#'  assimilated.
#' @param init Numeric vector of initial values for the MCMC chain.
#' @param n_iter Number of iterations for \code{\link{adaptMCMC}}.
#' @param n_chain Number of MCMC chains to be run (more chains facilitates
#'  covergence diagnostics).
#' @param parallel Boolean: Should multiple MCMC chains be run in parallel?
#' @param ... Additional parameters to be passed to the log-posterior
#'  function.
#' @return List of \code{\link{adaptMCMC}} output.
#' @export
run_mcmc <- function(post, parnames, residtype, prior_df, data_yrs, init, n_iter=1e6, n_chain=4, parallel=TRUE, ...) {

  ## prune data to specified years
  dat <- lapply(iamdata, function(l) {l[l$year %in% data_yrs,]})
  ## obtain parameter names based on type of residual strucure

  ## set up prior distributions
  priors <- create_prior_list(prior_df)
  
  ## set up MCMC run parameters
  rate_accept_many <- 0.234 # optimal acceptance rate for multiple parameters
  rate_accept_one <- 0.44 # optimal acceptance rate for one parameter
  rate_accept <- rate_accept_many + (rate_accept_one - rate_accept_many)/length(parnames)
  
  # set initial step size to be 10% of the standard deviation of the prior distributions when defined, or else some initial guess otherwise
  stepsize <- numeric(length(parnames))
  for (name in parnames) {
    if (priors[[name]][['type']] == 'uniform') {
      stepsize[match(name, parnames)] <- 0.1*(priors[[name]][['max']] - priors[[name]][['min']])/sqrt(12)
    } else if (priors[[name]][['type']] == 'normal') {
      stepsize[match(name, parnames)] <- 0.1*priors[[name]][['sd']]
    } else if (priors[[name]][['type']] == 'log-normal') {
      stepsize[match(name, parnames)] <- 0.1*sqrt((exp(priors[[name]][['sdlog']]^2)-1)*exp(2*priors[[name]][['meanlog']]+priors[[name]][['sdlog']]^2))
    }
  }

  adapt_start <- max(500, round(0.01*n_iter)) # when to start stepsize adaptation


  # if MCMC is run in parallel, set up cluster and run chains
  # otherwise just run chains
  if (parallel) {
    ncores <- parallel::detectCores()
    mcmc_out <- MCMC.parallel(post, n_iter, init=init, n.chain=n_chain, n.cpu=ncores, acc.rate=rate_accept, n.start = adapt_start, gamma=0.51, scale=stepsize, adapt=5e5, list=TRUE, packages=c('IAMUQ'), parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', residtype), ...)
  } else {
      mcmc_out <- MCMC(post, n_iter, init=init, n.chain=n_chain, acc.rate=rate_accept, n.start = adapt_start, gamma=0.51, scale=stepsize, adapt=5e5, list=TRUE, parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', residtype), ...)
  }

  # return MCMC output
  mcmc_out
}