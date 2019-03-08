## add ar1 noise
# simulate stationary AR(1) process (approximate - faster, better convergence, and
# results not sensitive to use of this as opposed to exact AR1)
ar.sim <- function(N,rho1,sigma) {
    x <- rep(NA,N)
    if(length(sigma) > 1) {
        x[1] = rnorm(n=1,sd=sigma[1]/sqrt(1-rho1^2))
        for (i in 2:N) {
            x[i] <- rho1*x[i-1] + rnorm(1,sd=sigma[i])
        }
    } else {
        x[1] = rnorm(n=1,sd=sigma/sqrt(1-rho1^2))
        for (i in 2:N) {
            x[i] <- rho1*x[i-1] + rnorm(1,sd=sigma)
        }
    }
    x
}

## add iid normal noise
iid.sim <- function(N, sigma) {
  rnorm(N, mean=0, sd=sigma)
}

## simulate model output
sim_model <- function(pars, parnames, type, log=F) {
  # evaluate model
  model_out <- mod(pars, parnames, start=1700, end=2100)[, c('year', 'P', 'Q', 'C')]
  # add simulated model discrepancy terms
  if (type == 'ar') {
    # extract AR process parameters
    rho <- pars[match(c('rho_pop', 'rho_prod', 'rho_emis'), parnames)] # AR model error coefficient
    names(rho) <- c('P', 'Q', 'C')
    sigma <- pars[match(c('sigma_pop', 'sigma_prod', 'sigma_emis'), parnames)] # AR model error sd
    names(sigma) <- c('P', 'Q', 'C')
    discrepancy <- vapply(names(rho), function(n) match.fun(paste0(type, '.sim'))(length(1700:2100), rho[n], sigma[n]), numeric(length(1700:2100)))
    if (log) {
      sim_out <- exp(log(model_out[, c('P', 'Q', 'C')]) + discrepancy)
    } else {
      sim_out <- model_out[, c('P', 'Q', 'C')] + discrepancy
      sim_out[sim_out < 0] <- 0
    }
  }
  data.frame(year=model_out[,'year'], sim_out)
}      

## evaluate posterior predictive statistics

## compute annual average gdp per capita growth rate over some period
avg_gwp_rate <- function(model_out, start=2010, end=2100) {
  end_yrs <- model_out[model_out$year %in% c(start, end),] # we only need values from start and end years
  gwp_pc <- end_yrs$Q / end_yrs$P # compute per capita GWP 
  exp((log(gwp_pc[2]) - log(gwp_pc[1]))/(end-start)) - 1 # compute average growth rate
}

## compute cumulative CO2 emissions over some period
cum_co2 <- function(model_out, start=2000, end=2100) {
  sum(model_out$C[model_out$year %in% start:end])
}