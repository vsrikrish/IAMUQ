##########################################################
# likelihood.R                                           #
#   This file contains the code for the likelihood,      #
#     prior, and posterior functions                     #
##########################################################

R_path <- 'R'
source(file.path(R_path, 'model.R'))

## evaluate log-prior densities
log_pri <- function(pars, parnames, priors) {

   # this function evaluates the log-prior density for a given parameter
   log_dens <- function(name) {
      val <- pars[match(name, parnames)]
      do.call(match.fun(priors[[name]][['dens.fun']]),
                  c(list(x=val, log=TRUE),
                    priors[[name]][-which(names(priors[[name]]) %in% c('type', 'dens.fun', 'quant.fun'))])
      )
    }


  # evaluate log-prior densities for each parameter
  lp <- vapply(parnames, log_dens, numeric(1))

  # return sum of log-priors
  sum(lp)
}

## calculate the model residuals
residuals <- function(pars, parnames, dat, start = 1700, end = 2017) {
  # run model
  model_out <- mod(pars, parnames, start, end)

  # compute residuals for each output
  r <- list()
  pop <- merge(dat[['pop']], model_out[, c('year', 'P')], by='year')
  r[['pop']] <- log(pop$pop) - log(pop$P)
  prod <- merge(dat[['prod']], model_out[, c('year', 'Q')], by='year')
  r[['prod']] <- log(prod$prod) - log(prod$Q)
  emis <- merge(dat[['emissions']], model_out[, c('year', 'C')], by='year')
  r[['emissions']] <- log(emis$emissions) - log(emis$C)

  # return residuals
  r
}

## calculate AR1 residuals
ar_residuals <- function(pars, parnames, dat, start=1700, end = 2017) {
  r <- residuals(pars, parnames, dat, start, end)

  rho <- pars[match(c('rho_pop', 'rho_prod', 'rho_emis'), parnames)] # AR model error coefficient
  names(rho) <- c('pop', 'prod', 'emissions')
  lapply(names(rho), function(n) {r[[n]][-1] - rho[n]*r[[n]][-length(r[[n]])]})
}

## evaluates the log-likelihood under the assumption of iid residuals
log_lik_iid <- function(pars, parnames, dat, start = 1700, end = 2017) {

  r <- residuals(pars, parnames, dat, start, end)

  # extract likelihood parameters
  sigma <- pars[match(c('sigma_pop', 'sigma_prod', 'sigma_emis'), parnames)]
  names(sigma) <- c('pop', 'prod', 'emissions')
  # compute log-likelihoods
  log_lik <- function(datname) {
    sum(dnorm(r[[datname]], mean = 0, sd = sigma[datname], log = TRUE))
  }
  # return sum of log-likelihood across data
  sum(vapply(names(sigma), log_lik, numeric(1)))
}

## evaluates the log-likelihood with AR(1) model errors and iid observation errors
log_lik_ar <- function(pars, parnames, dat, start = 1700, end = 2017) {
  # compute residuals
  r <- residuals(pars, parnames, dat, start, end)

  # extract likelihood parameters
  rho <- pars[match(c('rho_pop', 'rho_prod', 'rho_emis'), parnames)] # AR model error coefficient
  names(rho) <- c('pop', 'prod', 'emissions')
  sigma <- pars[match(c('sigma_pop', 'sigma_prod', 'sigma_emis'), parnames)] # AR model error sd
  names(sigma) <- c('pop', 'prod', 'emissions')
  eps <- pars[match(c('eps1_pop', 'eps1_prod', 'eps1_emis'), parnames)] # observation error coefficient
  names(eps) <- c('pop', 'prod', 'emissions')
#  eps2 <- pars[match(c('eps2_pop', 'eps2_prod', 'eps2_emis'), parnames)] # observation error coefficient

  # compute log-likelihoods
  log_lik <- function(datname) {
    # compute covariance matrix for likelihood
    n <- length(r[[datname]])
    H <- abs(outer(1:n, 1:n, '-'))

    V <- sigma[datname]^2 / (1-rho[datname]^2) * rho[datname]^H # variance for model errors
    epv <- c(rep(eps[datname]^2, (sum(dat[[datname]]$year < 1950))), rep(0, sum(dat[[datname]]$year >= 1950)))
    Sigma <- V + diag(epv)
    dmvnorm(r[[datname]], sigma = V, log = TRUE) # compute log-likelihood
  }
  # return sum of log-likelihood across data
  sum(vapply(names(rho), log_lik, numeric(1)))
}

## check the prior constraints on parameter values
check_constraints <- function(pars, parnames) {
   # check for parameter constraints; if not satisfied, return -Inf
   delta <- pars[match('delta', parnames)]
   s <- pars[match('s', parnames)]
   rho <- pars[match(c('rho2', 'rho3'), parnames)]
   tau <- pars[match(c('tau2', 'tau3', 'tau4'), parnames)]
   if (any(grepl('eps2', parnames))) {
    eps1 <- pars[match(c('eps1_pop', 'eps1_prod', 'eps1_emis'), parnames)]
     eps2 <- pars[match(c('eps2_pop', 'eps2_prod', 'eps2_emis'), parnames)]
     eps <- rbind(eps1, eps2)
     return((delta >= s) || (rho[1] < rho[2]) || (any(tau != cummax(tau))) || (any(apply(eps, 2, function(x) x[1] < x[2]))))
    }

   (delta < s) && (rho[1] >= rho[2]) && (all(tau == cummax(tau)))

}

## compute the negative log-likelihood (for use with DEoptim)
neg_log_lik <- function(pars, parnames, dat, lik_fun) {
   # check for parameter constraints; if not satisfied, return Inf
   if (!check_constraints(pars, parnames)) {
     return(Inf)
   }

  # evaluate log-likelihood
  ll <- match.fun(lik_fun)(pars, parnames, dat)
  # return negative log-likelihood
  -1*ll
}


## compute the log-posterior density
log_post <- function(pars, parnames, priors, dat, lik_fun) {
  # check for parameter constraints and return -Inf if not satisfied
   if (!check_constraints(pars, parnames))  {
     return(-Inf)
   }

  # compute log-priors
  lp <- log_pri(pars, parnames, priors)
  # if log-prior density is -Inf, no need to evaluate likelihood
  if (lp == -Inf) {
    return(-Inf)
  }
  # evaluate likelihood
  ll <- match.fun(lik_fun)(pars, parnames, dat)
  # return sum of log-prior and log-likelihood
  lp + ll
}
