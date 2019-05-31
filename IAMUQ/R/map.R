find_map <- function(post, parnames, residtype, prior_df, data_yrs, n_iter=5e3, NP_scale=25, parallel=TRUE, trace=FALSE, ...) {

  ## prune data to specified years
  dat <- lapply(baudata, function(l) {l[l$year %in% data_yrs,]})
  ## obtain parameter names based on type of residual structure
  # list all possible parameter names for DE bounds
  all_parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis', 'a_11', 'a_22', 'a_33', 'a_21', 'a_31', 'a_12', 'a_23', 'a_13', 'a_32', 'eps_pop', 'eps_prod', 'eps_emis')
  
  ## set up prior distributions
  priors <- create_prior_list(prior_df)
  
  ## set DEoptim parameters for MAP estimate
  # set upper and lower bounds for each variable
  lbound <- c(0.0001, 0, 6.9, 0.3, 0.5, 0.1, 0.01, 0.0007, 5.3, 0.4, 0, 0, 0, 1700, 1700, 2010, 0.005, 0, 0, 0, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0, 0, 0)
  ubound <- c(0.15, 100, 15, 1, 0.9,  0.3, 0.14, 0.0212, 16.11, 0.7, 3, 0.5, 0.5, 2100, 2100, 2500, 0.2, 5, 5, 5, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 5, 5, 5)


  # if DEoptim is run in parallel, set up cluster and evaluate
  # otherwise just evaluate
  if (parallel) {
    map_out <- DEoptim(post, lbound[match(parnames, all_parnames)], ubound[match(parnames, all_parnames)], control=list(parallelType=1, trace=trace, itermax=n_iter, NP=NP_scale*length(parnames)), parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', residtype), ...)
  } else {
      map_out <- DEoptim(post, lbound[match(parnames, all_parnames)], ubound[match(parnames, all_parnames)], control=list(NP=NP_scale*length(parnames), itermax=n_iter, trace=trace), parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', residtype), ...)
  }
  
  # set best estimate parameter names
  names(map_out$optim$bestmem) <- parnames

  # return MAP estimate
  map_out
}