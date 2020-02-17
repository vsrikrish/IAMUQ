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
#' @export
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
#' @export
cum_co2 <- function(emis, yrs, start=2000, end=2100) {
  yr_idx <- which(yrs %in% start:end)
  colSums(as.matrix(emis)[yr_idx, , drop=FALSE])
}
