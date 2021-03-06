% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood.R
\name{check_param_constraints}
\alias{check_param_constraints}
\title{Check if parameter values satisfy prior constraints}
\usage{
check_param_constraints(pars, parnames)
}
\arguments{
\item{pars}{Numeric vector of parameter values. These parameter values are the same that are passed to \code{\link{log_post}}.}

\item{parnames}{Character vector of parameter names. These names should
align with the values in \code{pars}, but they don't need to be in any
particular order otherwise.}
}
\value{
Boolean value corresponding to if the parameter values satisfied
 the given constraints.
}
\description{
\code{check_param_constraints} checks if the parameter values specify
 constraints which are specified ahead of time.
}
\details{
The specified constraints are (for the parameters passed to
 \code{\link{log_post}}): \itemize{
 \item delta < s;
 \item rho1 >= rho2;
 \item tau2 <= tau3 <= <= tau4;
 \item eigenvalues of the VAR coefficient matrix $A$ are <= 1 (for VAR
 residual model).
 }
}
