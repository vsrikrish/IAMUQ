#' \code{IAMUQ} package
#'
#' Code to run, calibrate, and simulate from the integrated assessment model
#'  in Srikrishnan, Guan, Tol, and Keller (2019).
#'
#' @doctype package
#' @name IAMUQ
#' @import Rcpp
#' @importFrom stats rnorm dnorm qnorm rlnorm dlnorm qlnorm
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom parallel detectCores
#' @importFrom adaptMCMC MCMC
#' @importFrom adaptMCMC MCMC.parallel
#' @importFrom DEoptim DEoptim
#' @importFrom magic adiag
#' @importFrom truncnorm qtruncnorm dtruncnorm rtruncnorm
#' @useDynLib IAMUQ, .registration = TRUE
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("iamdata"))
