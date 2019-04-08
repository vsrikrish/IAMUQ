##########################################################
# likelihood.R                                           #
#   This file contains the code for the likelihood,      #
#     prior, and posterior functions                     #
##########################################################

#### evaluate log-expert assessment likelihood for probabilistic inversion ####

## GWP expert assessment
# we use the trimmed-mean expert assessment distribution from
# Christensen et al (2018) for global per-capita growth rate from 2010 to 2100
log_exp_gwp <- function(model_out) {
  # set parameter values for expert assessment distribution
  mu=2.54 # mean on percentage scale
  sigma=1.07 # sd on percentage scale

  # compute per-capita growth rates from model output
  gwp_rt <- avg_gwp_rate(model_out, start=2010, end=2100)
  
  # compute and return log-likelihood of growth rate given expert assessment
  # distribution parameters are on a percentage scale, so multiply rate by 100
  dnorm(gwp_rt * 100, mean=mu, sd=sigma, log=TRUE)
}
  
## CO2 emissions expert assessment
# we use the median of medians and IQRs of CO2 emissions in 2100 from Ho et al (in prep)
# we assume a normal distribution for mapping median and IQR to a parametric distribution
log_exp_co2 <- function(model_out) {
  # set parameter values from Ho et al
  med <- 54.5 # median of individual expert medians
  iqr <- 58.2 # median of individual expert IQRs
  # compute parameters of normal approximation
  mu <- med # mean is equal to median
  sigma <- iqr / (qnorm(0.75) - qnorm(0.5)) / 2 # simplified equation using the 75% and 50% percentiles
  
  # compute and return log-likelihood of CO2 emissions given expert assessment
  dnorm(model_out$C[model_out$year == 2100], mean=mu, sd=sigma, log=TRUE)
}
  
## evaluate log-prior densities
log_pri <- function(pars, parnames, priors) {
  
  # this function evaluates the log-prior density for a given parameter
  log_dens <- function(name) {
    val <- pars[match(name, parnames)]
    do.call(match.fun(priors[[name]][['dens.fun']]),
            c(list(x=val, log=TRUE),
              priors[[name]][-which(names(priors[[name]]) %in% c('type', 'dens.fun', 'quant.fun', 'rand.fun'))])
    )
  }
  # evaluate log-prior densities for each parameter
  lp <- vapply(parnames, log_dens, numeric(1))
  
  # return sum of log-priors
  sum(lp)
}


## check the prior constraints on parameter values
check_param_constraints <- function(pars, parnames) {
  # check for parameter constraints; if not satisfied, return -Inf
  delta <- pars[match('delta', parnames)]
  s <- pars[match('s', parnames)]
  rho <- pars[match(c('rho2', 'rho3'), parnames)]
  tau <- pars[match(c('tau2', 'tau3', 'tau4'), parnames)]
  # check if the eigenvalues of the VAR coefficient matrix are within the unit circle (for stability of the VAR process)
  eig_flag <- TRUE
  if ('a_21' %in% parnames) {
    a <- pars[match(c('a_11', 'a_21', 'a_31', 'a_12', 'a_22', 'a_23', 'a_13', 'a_32', 'a_33'), parnames)]
    ev <- eigen(matrix(a, ncol=sqrt(length(a)), nrow=sqrt(length(a))))$values
    eig_flag <- all(vapply(ev, abs, numeric(1)) <= 1)
  }
  
  (delta < s) && (rho[1] >= rho[2]) && (all(tau == cummax(tau))) && eig_flag
  
}

## check if a fossil fuel emissions constraint is satisfied
check_fossil_constraint <- function(model_out, start=1700, end=2500, thresh) {
  fossil_tot <- cum_co2(model_out, start, end)
  if (is.na(fossil_tot)) {
    FALSE
  } else {
    fossil_tot <= thresh
  }
}


## calculate the model residuals
residuals <- function(pars, parnames, model_out, dat) {
  # compute residuals for each output
  r <- list()
  pop <- merge(dat[['pop']], model_out[, c('year', 'P')], by='year')
  r[['pop']] <- pop$pop - pop$P
  prod <- merge(dat[['prod']], model_out[, c('year', 'Q')], by='year')
  r[['prod']] <- prod$prod - prod$Q
  emis <- merge(dat[['emissions']], model_out[, c('year', 'C')], by='year')
  r[['emissions']] <- emis$emissions - emis$C
  
  # return residuals
  r
}

## evaluates the log-likelihood under the assumption of iid residuals
log_lik_iid <- function(pars, parnames, model_out, dat) {
  # check for fossil fuel constraint
  if (!check_fossil_constraint(model_out, start=1700, end=2500, thresh=6000)) {
    return(-Inf)
  }
  
  r <- residuals(pars, parnames, model_out, dat)
  
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
log_lik_ar <- function(pars, parnames, model_out, dat) {
  # check fossil fuel constraint
  if (!check_fossil_constraint(model_out, start=1700, end=2500, thresh=6000)) {
    return(-Inf)
  }
  
  # compute residuals
  r <- residuals(pars, parnames, model_out, dat)
  
  # extract likelihood parameters
  a <- pars[match(c('a_pop', 'a_prod', 'a_emis'), parnames)] # AR model error coefficient
  names(a) <- c('pop', 'prod', 'emissions')
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
    
    V <- sigma[datname] / (1-a[datname]^2) * a[datname]^H # variance for model errors
    epv <- c(rep(eps[datname], (sum(dat[[datname]]$year < 1950))), rep(0, sum(dat[[datname]]$year >= 1950)))
    Sigma <- V + diag(epv)
    dmvnorm(r[[datname]], sigma = V, log = TRUE) # compute log-likelihood
  }
  # return sum of log-likelihood across data
  sum(vapply(names(a), log_lik, numeric(1)))
}

## evaluates the log-likelihood with AR(1) model errors and iid observation errors
log_lik_var <- function(pars, parnames, model_out, dat) {
  # check fossil fuel constraint
  if (!check_fossil_constraint(model_out, start=1700, end=2500, thresh=6000)) {
    return(-Inf)
  }
  
  # compute residuals
  r <- residuals(pars, parnames, model_out, dat)
  
  # create vectorized residual vector
  r_vec <- matrix(do.call(rbind,r), nrow=1)
  
  # extract likelihood parameters
  a <- pars[match(c('a_11', 'a_21', 'a_31', 'a_12', 'a_22', 'a_32', 'a_13', 'a_23', 'a_33'), parnames)] # VAR model error coefficient
  sigma <- pars[match(c('sigma_pop', 'sigma_prod', 'sigma_emis'), parnames)] # AR model error sd
  names(sigma) <- c('pop', 'prod', 'emissions')
  eps <- pars[match(c('eps_pop', 'eps_prod', 'eps_emis'), parnames)] # observation error coefficient
  names(eps) <- c('pop', 'prod', 'emissions')
  #  eps2 <- pars[match(c('eps2_pop', 'eps2_prod', 'eps2_emis'), parnames)] # observation error coefficient
  
  # construct VAR coefficient matrix
  A <- matrix(a, nrow=length(dat), ncol=length(dat))
   
  # compute covariance matrix of the residuals for each time
  W <- diag(sigma)  # construct covariance matrix of the innovations
  Sigma_x_vec <- solve(diag(1, nrow(A)^2) - kronecker(A, A)) %*% as.numeric(W)
  Sigma_x <- matrix(Sigma_x_vec, nrow=nrow(A), ncol=ncol(A))
  
  # compute powers of A for autocovariance matrix computation
  n <- length(r[[1]])

  H <- abs(outer(1:n, 1:n, '-')) # matrix of indices for blocks as they should be combined
  D <- diag(eps) # matrix of observation error variances
  
  Sigma <- cov_mat(A, Sigma_x, D, H) # generate
  
  # compute log-likelihood of residuals with respect to zero mean
  dmvnorm(r_vec, sigma=Sigma, log=TRUE)
}

## compute the log-posterior density
log_post <- function(pars, parnames, priors, dat, lik_fun, exp_co2=FALSE, exp_gwp=FALSE) {
  # check for parameter constraints and return -Inf if not satisfied
  if (!check_param_constraints(pars, parnames))  {
    return(-Inf)
  }
  
  # compute log-priors
  lpri <- log_pri(pars, parnames, priors)
  # if log-prior density is -Inf, no need to evaluate likelihood
  if (lpri == -Inf) {
    return(-Inf)
  }
  # run model to end date, which is 2500
  yr <- 1700:2500
  model_out <- run_model(pars, parnames, start=1700, end=2500)
  ll <- match.fun(lik_fun)(pars, parnames, model_out, dat) # evaluate likelihood
  
  lpost <- lpri + ll # store sum of log-likelihood and log-prior

  # if expert assessment inversion is included, add log-expert assessment density
  if (exp_gwp) {
    lpost <- lpost + log_exp_gwp(model_out)
  }
  if (exp_co2) {
    lpost <- lpost + log_exp_co2(model_out)
  }

  lpost # return log-posterior vlaue
}

## compute the negative log-posterior (for use with DEoptim)
neg_log_post <- function(pars, parnames, priors, dat, lik_fun, exp_co2=FALSE, exp_gwp=FALSE) {

  # evaluate log-likelihood
  lp <- log_post(pars,parnames, priors, dat, lik_fun, exp_co2, exp_gwp)
  # return negative log-likelihood
  -1*lp
}

