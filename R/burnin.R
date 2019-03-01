library(coda)

# function to compute the burn-in length based on the Gelman-Rubin diagnostic of multiple parallel chains

burnin <- function(mcmc_out, niter_seq) {

  # convert the chains to MCMC objects
  mcmc_samp <- lapply(mcmc_out, function(l) l$samples)
  # compute G-R diagostic values
  
  vapply(niter_seq, function(n) {as.numeric(gelman.diag(lapply(mcmc_samp, function(l) as.mcmc(l[1:n,])))[2])}, numeric(1))
}
