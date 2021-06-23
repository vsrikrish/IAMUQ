#' Create list of prior density information based on input data frame.
#'
#' \code{create_prior_list} creates a list for use in \code{\link{log_pri}}.
#'  The parameters for each prior density function are estimated from the
#'  confidence interval bounds in the data frame. Available (by default)
#'  density functions are normal, uniform, log-normal, and truncated normals (or truncnorm).
#'
#' @param prior_df Data frame of information about the prior density
#'  functions. Each row (corresponding to a single parameter) should have the
#'  following entries columns:\itemize{
#'  \item name, the parameter name;
#'  \item type, the type of prior density function. Options are "normal,"
#'  "uniform," "log-normal," and "truncnorm." Truncated normals are by default truncated
#'  below at 2020;
#'  \item lower, numeric value defining the lower value. If the type is
#'  "normal," "log-normal," or "truncnorm," this should be the 0.05 quantile (or the lower
#'  bound of the 90% credible interval). If the type is "uniform," this
#'  should be the lower bound of the distribution;
#'  \item upper, numeric value defining the upper value. If the type is
#'  "normal," "log-normal," or "truncnorm," this should be the 0.95 quantile (or the upper
#'  bound of the 90% credible interval). If the type is "uniform," this
#'  should be the upper bound of the distribution. }
#' @return List with each entry corresponding to an individual parameter.
#'  Each of these lists has the following entries: \itemize{
#'  \item type, the parameter type;
#'  \item dens.fun, string giving the name of the function which evaluates
#'  the density;
#'  \item quant.fun, string giving the name of the function which computes
#'  quantiles;
#'  \item rand.fun, string giving the name of the function which generates
#'  random samples;
#'  \item other parameters required for the various function calls, which
#'  depend on the type of distribution.
#'  }
#' @export
create_prior_list <- function(prior_df) {
  parnames <- prior_df[, 'name']
  priors <- vector('list', length(parnames))
  names(priors) <- parnames
  
  for (i in 1:nrow(prior_df)) {
    name <- prior_df[i, 'name']
    priors[[name]] <- list(type=prior_df[i, 'type'])
    if (priors[[name]][['type']] == 'uniform') {
      priors[[name]][['dens.fun']] <- 'dunif'
      priors[[name]][['quant.fun']] <- 'qunif'
      priors[[name]][['rand.fun']] <- 'runif'
      priors[[name]][['min']] <- prior_df[i, 'lower']
      priors[[name]][['max']] <- prior_df[i, 'upper']
    } else if (priors[[name]][['type']] == 'normal') {
      priors[[name]][['dens.fun']] <- 'dnorm'
      priors[[name]][['quant.fun']] <- 'qnorm'
      priors[[name]][['rand.fun']] <- 'rnorm'
      priors[[name]][['mean']] <- mean(c(prior_df[i, 'lower'], prior_df[i, 'upper']))
      priors[[name]][['sd']] <- (prior_df[i, 'upper'] - prior_df[i, 'lower'])/(qnorm(0.975) - qnorm(0.025))
    } else if (priors[[name]][['type']] == 'log-normal') {
      priors[[name]][['dens.fun']] <- 'dlnorm'
      priors[[name]][['quant.fun']] <- 'qlnorm'
      priors[[name]][['rand.fun']] <- 'rlnorm'
      priors[[name]][['meanlog']] <- mean(c(log(prior_df[i, 'lower']), log(prior_df[i, 'upper'])))
      priors[[name]][['sdlog']] <- (log(prior_df[i, 'upper']) - log(prior_df[i, 'lower'])) / (qnorm(0.95) - qnorm(0.05))
    } else if (priors[[name]][['type']] == 'truncnorm') {
      priors[[name]][['dens.fun']] <- 'dtruncnorm'
      priors[[name]][['quant.fun']] <- 'qtruncnorm'
      priors[[name]][['rand.fun']] <- 'rtruncnorm'
      priors[[name]][['a']] <- 2020
      priors[[name]][['mean']] <- mean(c(prior_df[i, 'lower'], prior_df[i, 'upper']))
      priors[[name]][['sd']] <- (prior_df[i, 'upper'] - prior_df[i, 'lower'])/(qnorm(0.975) - qnorm(0.025))
    }
  }
  
  priors
}
