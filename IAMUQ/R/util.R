## add iid normal noise
iid.sim <- function(N, sigma) {
  rnorm(N, mean=0, sd=sigma)
}

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

## add VAR1 noise
# simulate stationary VAR(1) process
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

## simulate model output
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
    sim_out <- model_out[, c('P', 'Q', 'C')]
  }
  data.frame(year=model_out[,'year'], sim_out)
}

## compute annual average gdp per capita growth rate over some period
avg_gwp_rate <- function(pop, gwp, yrs, start=2010, end=2100) {
  bdry_idx <- which(yrs %in% c(start, end)) # we only need values from start and end years
  gwp <- as.matrix(gwp)
  pop <- as.matrix(pop)
  gwp_pc <- gwp[bdry_idx, , drop=FALSE] /  pop[bdry_idx, , drop=FALSE]# compute per capita GWP
  exp((log(gwp_pc[2,]) - log(gwp_pc[1,]))/(end-start)) - 1 # compute average growth rate
}

## compute cumulative CO2 emissions over some period
cum_co2 <- function(emis, yrs, start=2000, end=2100) {
  yr_idx <- which(yrs %in% start:end)
  colSums(as.matrix(emis)[yr_idx, , drop=FALSE])
}
