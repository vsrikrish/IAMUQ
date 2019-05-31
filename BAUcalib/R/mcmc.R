run_mcmc <- function(post, parnames, residtype, prior_df, data_yrs, init, n_iter=1e6, n_chain=4, parallel=TRUE, ...) {

  ## prune data to specified years
  dat <- lapply(baudata, function(l) {l[l$year %in% data_yrs,]})
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
    mcmc_out <- MCMC.parallel(post, n_iter, init=init, n.chain=n_chain, n.cpu=ncores, acc.rate=rate_accept, n.start = adapt_start, gamma=0.51, scale=stepsize, adapt=5e5, list=TRUE, packages=c('BAUcalib'), parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', residtype), ...)
  } else {
      mcmc_out <- MCMC(post, n_iter, init=init, n.chain=n_chain, acc.rate=rate_accept, n.start = adapt_start, gamma=0.51, scale=stepsize, adapt=5e5, list=TRUE, parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', residtype), ...)
  }

  # return MCMC output
  mcmc_out
}