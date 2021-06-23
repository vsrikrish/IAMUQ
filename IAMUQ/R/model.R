#' Run the integrated assessment model from Srikrishnan & Keller (2019).
#'
#' \code{run_model} runs the integrated assessment model over the specified
#'    period (from start to end).
#'
#' This function runs the model using the specified parameters. It can be
#'  started with arbitrary initial conditions, though these are detected
#'  based on whether or not \code{P0} (as a function argument) is
#'  NULL.
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
#'  \item kappa, the rate of technological penetration.
#'  }
#' @param parnames Character vector of parameter names. These names should
#'  align with the values in \code{pars}, but they don't need to be in any
#'  particular order otherwise.
#' @param start Numeric starting year for the simulation.
#' @param end Numeric ending year for the simulation.
#' @param P0 Optional initial population value (if not treated as a
#'  parameter).
#' @param Q0 Optional initial economic output value. If not specified, the
#'  initial value is treated initialized based on passed parameter values.
#' @param C0 Optional initial carbon emissions value. If not specified,
#'  the initial value is calculated based on the economic output and
#'  emissions parameters.
#' @return Data frame with columns \itemize{
#'  \item \code{year} (model years, between \code{start} and \code{end});
#'  \item \code{P} (population, in billions);
#'  \item \code{Q} (gross world product, or economic output, in trillions
#'  USD$2011);
#'  \item \code{C} (CO2 emissions, in Gt C/yr).}
#' @export
run_model <- function(pars, parnames, start=1700, end=2017, P0=NULL, Q0=NULL, C0=NULL) {
  # create time vector
  yr <- seq(start, end)
  n_yr <- length(yr)
  
  # extract parameter values
  if (!is.null(P0)) {
    init <- c(Q0, C0)
  } else {
    init <- NULL
    P0 <- pars[match('P0', parnames)]
  }
  psi <- pars[match(c('psi1', 'psi2', 'psi3'), parnames)]
  alpha <- pars[match('alpha', parnames)]
  A0 <- pars[match('A0', parnames)]
  As <- pars[match('As', parnames)]
  s <- pars[match('s', parnames)]
  lambda <- pars[match('lambda', parnames)]
  delta <- pars[match('delta', parnames)]
  pi <- pars[match('pi', parnames)]
  kappa <- pars[match('kappa', parnames)]
  tau <- pars[match(c('tau2', 'tau3', 'tau4'), parnames)]
  rho <- pars[match(c('rho2', 'rho3'), parnames)]
  
  
  # compute technology penetration values for technologies with non-zero emissions
  gamma <- cbind(
    1 - (1 / (1 + exp(-kappa * (yr - tau[1])))),
    1 / (1 + exp(-kappa * (yr - tau[1]))) - 1 / (1 + exp(-kappa * (yr - tau[2]))),
    1 / (1 + exp(-kappa * (yr - tau[2]))) - 1 / (1 + exp(-kappa * (yr - tau[3]))),
    1 / (1 + exp(-kappa * (yr - tau[3])))
  )
  
  model_out <- model_run(yr, P0, psi, alpha, A0, As, s, lambda, delta, pi, kappa, gamma[, 2:3], rho, init)

  colnames(gamma) <- c('Frac_PreIndustrial', 'Frac_FossilHi', 'Frac_FossilLo', 'Frac_NonFossil')
  
  data.frame(model_out, gamma)  
}
