#' Find the maximum a posteriori (MAP) estimate for the IAM.
#'
#' \code{find_map} finds a MAP estimate using a differential evolution
#'  function (in this case, \code{DEoptim}).
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
#' @param n_iter Number of iterations for \code{\link{DEoptim}}.
#' @param NP_scale Numeric factor for the number of particles (relative to
#'  the number of parameters).
#' @param parallel Boolean: Should \code{\link{DEoptim}} be run in parallel?
#' @param trace Boolean: Should \code{\link{DEoptim}} print its trace?
#' @param ... Other parameters passed to the log-posterior function.
#' @return Output of \code{\link{DEoptim}} call. The best parameter vector is
#'  named using parnames.
#' @export
find_map <- function(post, parnames, residtype, prior_df, data_yrs, n_iter=5e3, NP_scale=25, parallel=TRUE, trace=FALSE, ...) {

  ## prune data to specified years
  dat <- lapply(iamdata, function(l) {l[l$year %in% data_yrs,]})
  ## obtain parameter names based on type of residual structure
  # list all possible parameter names for DE bounds
  all_parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis', 'a_11', 'a_22', 'a_33', 'a_21', 'a_31', 'a_12', 'a_23', 'a_13', 'a_32', 'eps_pop', 'eps_prod', 'eps_emis')
  
  ## set up prior distributions
  priors <- create_prior_list(prior_df)
  
  ## set DEoptim parameters for MAP estimate
  # set upper and lower bounds for each variable
  lbound <- c(0.0001, 0, 6.9, 0.3, 0.5, 0.1, 0.01, 0.0007, 5.3, 0.4, 0, 0, 0, 1850, 1900, 2030, 0.005, 0, 0, 0, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0, 0, 0)
  ubound <- c(0.15, 100, 15, 1, 0.9,  0.3, 0.14, 0.0212, 16.11, 0.7, 3, 0.5, 0.5, 1950, 2100, 2100, 0.2, 5, 5, 5, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 5, 5, 5)


  # if DEoptim is run in parallel, set up cluster and evaluate
  # otherwise just evaluate
  if (parallel) {
    map_out <- DEoptim(post, lbound[match(parnames, all_parnames)], ubound[match(parnames, all_parnames)], control=list(parallelType=1, trace=trace, itermax=n_iter, NP=NP_scale*length(parnames), packages=c('truncnorm')), parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', residtype), ...)
  } else {
      map_out <- DEoptim(post, lbound[match(parnames, all_parnames)], ubound[match(parnames, all_parnames)], control=list(NP=NP_scale*length(parnames), itermax=n_iter, trace=trace), parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', residtype), ...)
  }
  
  # set best estimate parameter names
  names(map_out$optim$bestmem) <- parnames

  # return MAP estimate
  map_out
}
