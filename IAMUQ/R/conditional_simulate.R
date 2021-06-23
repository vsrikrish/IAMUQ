#' Conditional simulation of out-of-sample data.
#'
#' \code{cond_sim_model} simulates out-of-sample data (held-out or
#'  projections) conditional on the observations.
#'
#' The joint likelihood function for a VAR(1) model is a multivariate
#'  Gaussian. \code{cond_sim_model} partitions the data into observations,
#'  which are treated as fixed, held-out data, and future projections. The
#'  model discrepancies of the simulation follows the desired joint
#'  likelihood. The variance of the held-out data takes into account
#'  observation error variances, while that is not included for projections.
#'
#' At least one of projyrs (vector of years for projections) and hoyrs
#'  (vector of years for held-out data) must be passed. It is ok for both to
#'  be passed.
#'
#' @param pars Numeric vector of parameter values. These parameters must
#'  include all of the model parameters listed in Table S1 and the
#'  statistical parameters listed in Table S2 of Srikrishnan &
#'  Keller (2019), with Greek letters spelled out: \itemize{
#'  \item psi1, the population growth rate;
#'  \item psi2, the population half-saturation constant;
#'  \item psi3, the population carrying capacity;
#'  \item P0, the initial population in year \code{start};
#'  \item lambda, the elasticity of production with respect to labor (must
#'  be less than 1);
#'  \item s, the savings rate;
#'  \item delta, the capital depreciation rate (must be less than s);
#'  \item alpha, the rate of technological progress for total factor
#'  productivity;
#'  \item As, the saturation level of total factor productivity;
#'  \item pi, the labor participation rate (must be less than 1);
#'  \item A0, the initial total factor productivity in year \code{start};
#'  \item rho2, the carbon emissions intensity of technology 2;
#'  \item rho3, the carbon emissions intensity of technology 3;
#'  \item tau2, the half-saturation year of technology 2;
#'  \item tau3, the half-saturation year of technology 3;
#'  \item tau4, the half-saturation year of technology 4;
#'  \item kappa, the rate of technological penetration;
#'  \item aij, for $i,j=1, 2, 3$, the elements of the VAR cefficient matrix
#'  (only when using a VAR model for the likelihood structure);
#'  \item sigma_pop, the variance of the VAR or normal innovations for
#'  population;
#'  \item sigma_prod, the variance of the VAR or normal innovations for
#'  economic output;
#'  \item sigma_emis, the variance of the VAR or normal innovations for
#'  emissions;
#'  \item eps_pop, the variance of the observation errors for population;
#'  \item eps_prod, the variance of the observation errors for economic
#'  output;
#'  \item eps_emis, the variance of the observation errors for emissions;
#'  }
#' @param parnames Character vector of parameter names. These names should
#'  align with the values in \code{pars}, but they don't need to be in any
#'  particular order otherwise.
#' @param dat List of data frames of data. List should have names 'pop'
#'  (population), 'prod' (production), and 'emissions,' and each data frame
#'  should have two columns, 'year' and 'value'.
#' @param projyrs Vector of years for projections. Should start in 2015 and
#'  be contiguous.
#' @param hoyrs Vector of years to hold out. These can be any years contained
#'  in the data.
#' @return Data frame of model simulation, with columns 'year, 'P'
#'  (population), 'Q' (economic output), and 'C' (CO2 emissions). Only the
#'  held-out and projected simulation years are returned.
#' @export
cond_sim_model <- function(pars, parnames, dat, projyrs=NULL, hoyrs=NULL) {
  # issue error if neither projyrs nor hoyrs is defined
  if (is.null(projyrs) & is.null(hoyrs)) {
    stop('Must specify either projyrs or hoyrs.')
  }

  # run model
  model_out <- run_model(pars, parnames, start=1700, end=max(c(dat[[1]]$year, projyrs)))
  
  # remove held-out data
  if (!is.null(hoyrs)) {
    d <- lapply(dat, function(l) l[-which(l$year %in% hoyrs), ])
  } else {
    d <- dat
  }
  # compute residuals
  r <- residuals(model_out, d)
  
  # create vectorized residual vector
  r_vec <- as.numeric(matrix(do.call(rbind,r), nrow=1))
  
  # extract parameters
  a <- pars[match(c('a_11', 'a_21', 'a_31', 'a_12', 'a_22', 'a_32', 'a_13', 'a_23', 'a_33'), parnames)] # VAR model error coefficient
  sigma <- pars[match(c('sigma_pop', 'sigma_prod', 'sigma_emis'), parnames)] # AR model error sd
  names(sigma) <- c('pop', 'prod', 'emissions')
  eps <- pars[match(c('eps_pop', 'eps_prod', 'eps_emis'), parnames)] # observation error coefficient
  names(eps) <- c('pop', 'prod', 'emissions')
  
  # construct VAR coefficient matrix
  A <- matrix(a, nrow=length(d), ncol=length(d))
   
  # compute covariance matrix of the residuals for each time
  W <- diag(sigma)  # construct covariance matrix of the innovations
  Sigma_x_vec <- solve(diag(1, nrow(A)^2) - kronecker(A, A)) %*% as.numeric(W)
  Sigma_x <- matrix(Sigma_x_vec, nrow=nrow(A), ncol=ncol(A))
  
  # compute powers of A for autocovariance matrix computation
  m <- length(c(dat[[1]]$year, projyrs)) # length of complete time series

  H <- abs(outer(1:m, 1:m, '-')) # matrix of indices for blocks as they should be combined
  D <- diag(eps) # matrix of observation error variances
  
  Sigma <- cov_mat(A, Sigma_x, D, H) # generate covariate matrix
  
  # partition covariance matrices into blocks based on data
  # define indices
  n <- nrow(dat[[1]]) # length of series given by data
  q <- length(dat) * n # length of vectorized series given by data
  N <- length(dat) * m # length of vectorized total series
  hoyridx <- which(dat[[1]]$year %in% hoyrs) # if any data points are held out, find indices with respect to years
  # set indices for data and projections, assuming no held out values
  datidx <- 1:q
  if (N > q) {
    projidx <- (q+1):N
  } else {
    projidx <- integer(0)
  }
  # if there are held out years, find indices with respect to the vectorized series, remove them from the data index vector, and add them to the projection index vector
  hoidx <- numeric(0)
  if (length(hoyridx) > 0) {
    for (i in 1:length(dat)) {
      hoidx <- c(hoidx, length(dat)*(hoyridx-1)+i)
    }
    hoidx <- sort(hoidx) # sort to group same year data together
    datidx <- datidx[-hoidx]
  }
  hoprojidx <- c(hoidx, projidx)
  Sigma_11 <- Sigma[datidx, datidx]
  Sigma_12 <- Sigma[datidx, hoprojidx]
  Sigma_21 <- Sigma[hoprojidx, datidx]
  Sigma_22 <- Sigma[hoprojidx, hoprojidx]
  # remove observation errors from projections, but not from held out data
  homat <- matrix(0, nrow=length(hoidx), ncol=length(hoidx))
  projmat <- kronecker(diag(1, length(projidx)/length(d)), D)
  Sigma_22 <- Sigma_22 - magic::adiag(homat, projmat)
  
  # compute conditional means and covariance
  mu_factor <- Sigma_21 %*% solve(Sigma_11)
  mu_cs <- mu_factor %*% r_vec
  Sigma_cs <- Sigma_22 - (mu_factor %*% Sigma_12)
  L_cs <- t(chol(Sigma_cs)) # lower-diagonal Cholesky of Sigma_cs
  # sample
  discrepancy <- matrix(mu_cs + L_cs %*% rnorm(N-length(datidx)), ncol=length(dat), byrow=TRUE) # compute random variates and re-form them into columns
  sim_out <- exp(log(model_out[model_out$year %in% c(hoyrs, projyrs), c('P', 'Q', 'C')]) + discrepancy) # add discrepancy terms to model output
  data.frame(year=sort(c(hoyrs, projyrs)), sim_out, model_out[model_out$year %in% c(hoyrs, projyrs), c('Frac_PreIndustrial', 'Frac_FossilHi', 'Frac_FossilLo', 'Frac_NonFossil')])
  
}