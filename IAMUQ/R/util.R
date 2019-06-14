#' Simulate white noise.
#'
#' \code{iid.sim} is a wrapper for simulating white noise. The samples are
#'  drawn from a univariate normal distribution with mean zero.
#'
#' @param N Length of white noise sequence.
#' @param sigma Standard deviation for the white noise.
#' @return Vector of generated white noise samples.
iid.sim <- function(N, sigma) {
  rnorm(N, mean=0, sd=sigma)
}

#' Simulate VAR(1) noise.
#'
#' \code{var.sim} simulates vector autoregressive of order 1 (VAR(1)) noise.
#'
#' @param N Length of VAR(1) noise sequence.
#' @param A VAR coefficient (square) matrix.
#' @param W VAR innovation (square) matrix.
#' @return Matrix of generated noise samples. Each row corresponds to a
#'  different variable.
var.sim <- function(N, A, W) {
  x <- matrix(NA, nrow=ncol(A), ncol=N)
  
  Sigma_x_vec <- solve(diag(1, nrow(A)^2) - kronecker(A, A)) %*% as.numeric(W)
  Sigma_x <- matrix(Sigma_x_vec, nrow=nrow(A), ncol=ncol(A))
  
  x[,1] <- mvtnorm::rmvnorm(n=1, sigma=Sigma_x)
  for (i in 2:N) {
    x[,i] <- A %*% x[,i-1] + t(mvtnorm::rmvnorm(1, sigma=W))
  }
  x
}

#' Simulate IAM based on provided parameters.
#'
#' \code{sim_model} simulates the IAM, including statistical noise, over the
#'  desired period of time. If initial parameters P0, Q0, and C0 are passed,
#'  the model will be initialized with them; if they are NULL, the model will
#'  be initialized as in Srikrishnan and Keller (2019), with an uncertain
#'  initial population.
#'
#' @param pars Numeric vector of parameter values. These parameters must
#'  include all of the model parameters listed in Table S1 of Srikrishnan &
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
#'  \item eps_pop, the variance of the observation errors for population
#'  (only when using a VAR model for the likelihood structure);
#'  \item eps_prod, the variance of the observation errors for economic
#'  output (only when using a VAR model for the likelihood structure);
#'  \item eps_emis, the variance of the observation errors for emissions.
#'  (only when using a VAR model for the likelihood structure);
#'  }
#' @param parnames Character vector of parameter names. These names should
#'  align with the values in \code{pars}, but they don't need to be in any
#'  particular order otherwise.
#' @param type String with the type of residual structure used for the
#'  statistical noise. This can be "iid," "var," or "none" (to just return
#'  model output with no noise).
#' @param start Starting year for the model (by default, 1700).
#' @param end Ending year for the model.
#' @param P0 Optional initial value for global population (in billions).
#' @param Q0 Optional initial value for gross world product (in trillions
#'  USD2011)
#' @param C0 Optional initial value for carbon emissions (in Gt C/yr).
#' @return Data frame with columns "year" (model year), "P" (global
#'  population in billions), "Q" (gross world product in trillions USD2011),
#'  and "C" (carbon emissions in Gt C/yr).
sim_model <- function(pars, parnames, type, start=1700, end=2100, P0=NULL, Q0=NULL, C0=NULL) {
  # evaluate model
  model_out <- run_model(pars, parnames, start, end, P0, Q0, C0)[, c('year', 'P', 'Q', 'C')]
  # add simulated model discrepancy terms
  if (type == 'iid') {
    sigma <- pars[match(c('sigma_pop', 'sigma_prod', 'sigma_emis'), parnames)] # model error sd
    names(sigma) <- c('P', 'Q', 'C')
    discrepancy <- vapply(sigma, function(s) iid.sim(length(start:end), s), numeric(length(start:end)))
    sim_out <- model_out[, c('P', 'Q', 'C')] + discrepancy
  } else if (type == 'var') {
    # extract VAR process parameters
    a <- pars[match(c('a_11', 'a_21', 'a_31', 'a_12', 'a_22', 'a_32', 'a_13', 'a_23', 'a_33'), parnames)] # VAR model error coefficient
    sigma <- pars[match(c('sigma_pop', 'sigma_prod', 'sigma_emis'), parnames)] # AR model error sd

    # construct VAR coefficient matrix
    A <- matrix(a, nrow=length(sigma), ncol=length(sigma))
  
    W <- diag(sigma)  # construct covariance matrix of the innovations
    discrepancy <- t(var.sim(length(start:end), A, W))
    sim_out <- exp(log(model_out[, c('P', 'Q', 'C')]) + discrepancy)
  } else if (type == 'none') {
    # if "none", just return the model output
    sim_out <- model_out[, c('P', 'Q', 'C')]
  }
  # return data frame
  data.frame(year=model_out[,'year'], sim_out)
}

#' Compute the average per-capita gross world product growth rate over the
#'  specified years.
#'
#' \code{avg_gwp_rate} computes the average per-capita gross world product
#'  (GWP) growth rate for a model simulation over the specified years.
#'
#' @param pop Numeric vector of global population.
#' @param gwp Numeric vector of gross world products.
#' @param yrs Numeric vector of years corresponding to the values in pop and
#'  gwp.
#' @param start Starting year for the average per-capita GWP growth rate
#'  calculation.
#' @param end Ending year for the average per-capita GWP growth rate
#'  calculation.
#' @return Average per-capita GWP growth rate value (on a decimal scale, not
#'  a percentage scale.)
avg_gwp_rate <- function(pop, gwp, yrs, start=2010, end=2100) {
  bdry_idx <- which(yrs %in% c(start, end)) # we only need values from start and end years
  gwp <- as.matrix(gwp)
  pop <- as.matrix(pop)
  gwp_pc <- gwp[bdry_idx, , drop=FALSE] /  pop[bdry_idx, , drop=FALSE]# compute per capita GWP
  exp((log(gwp_pc[2,]) - log(gwp_pc[1,]))/(end-start)) - 1 # compute average growth rate
}

#' Compute the cumulative CO2 emissions over the specified years.
#'
#' \code{cum_co2} computes the cumulative CO2 emissions for a model
#'  simulation over the specified years.
#'
#' @param emis Numeric vector of carbon emissions.
#' @param yrs Numeric vector of years corresponding to the values in pop and
#'  gwp.
#' @param start Starting year for the average per-capita GWP growth rate
#'  calculation.
#' @param end Ending year for the average per-capita GWP growth rate
#'  calculation.
#' @return Cumulative CO2 emissions (in Gt C) over the specified period.
cum_co2 <- function(emis, yrs, start=2000, end=2100) {
  yr_idx <- which(yrs %in% start:end)
  colSums(as.matrix(emis)[yr_idx, , drop=FALSE])
}
