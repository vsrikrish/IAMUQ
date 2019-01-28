##########################################################
#   This file contains the code for the coupled          #
# model.R                                                #
#     population-economic-emissions model.               #
##########################################################

update_pop <- function(prev_vals, psi) {
  # compute previous time step per-capita consumption
  y <- prev_vals$Q / prev_vals$P
  # update population time step
  prev_vals$P * (1 + (psi[1] * y / (psi[2] + y)) * ((psi[3] - prev_vals$P) / psi[3]))
}

update_prod <- function(prev_vals, P, alpha, A0, As, s, lambda, delta, pi) {
  # update total factor productivity
  A <- prev_vals$A + alpha * prev_vals$A * (1 - (prev_vals$A / As))
  # compute labor input
  L <- pi * P
  # compute capital stock
  K <- (1 - delta) * prev_vals$K + s * prev_vals$Q
  # compute updated total world production Q
  Q <- A * (L ^ lambda) * (K ^ (1-lambda))
  # return list of A, L, K, and Q
  list(A=A, L=L, K=K, Q=Q)
}

update_emis <- function(Q, gamma, rho) {
  # compute carbon intensity
  phi <- sum(gamma*rho)
  # return carbon emissions
  Q*phi
}

mod <- function(pars, parnames, start=1700, end=2017) {
  # create time vector
  yr <- seq(start, end)
  # create data frame for storage of model values
  df <- data.frame(year=yr, P=numeric(length(yr)), Q=numeric(length(yr)), A=numeric(length(yr)), K=numeric(length(yr)), L=numeric(length(yr)), C=numeric(length(yr)))
  
  # extract parameter values
  P0 <- pars[match('P0', parnames)]
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
    1 / (1 + exp(-kappa * (yr - tau[1]))) - 1 / (1 + exp(-kappa * (yr - tau[2]))),
    1 / (1 + exp(-kappa * (yr - tau[2]))) - 1 / (1 + exp(-kappa * (yr - tau[3])))
  )
  
  # set initial values
  L0 <- pi * P0
  K0 <- L0 * (s * A0 / delta) ^ (1/lambda)
  Q0 <- A0 * (L0 ^ lambda) * (K0 ^ (1-lambda))
  C0 <- sum(gamma[1, ] * rho) * Q0
  df[1, -1] <- c(P0, Q0, A0, K0, L0, C0)
  
  # run model
  for (t in 2:length(yr)) {
    # compute updated population
    df[t, 'P'] <- update_pop(df[t-1, c('P', 'Q')], psi)
    # compute updated production
    prod_out <- update_prod(df[t-1, c('P', 'Q', 'A', 'K')], df[t, 'P'], alpha, A0, As, s, lambda, delta, pi)
    df[t, 'Q'] <- prod_out$Q
    df[t, 'A'] <- prod_out$A
    df[t, 'K'] <- prod_out$K
    df[t, 'L'] <- prod_out$L
    # compute updated emissions
    df[t, 'C'] <- update_emis(df[t, 'Q'], gamma[t, ], rho)
  }
  
  # return model output
  df
}