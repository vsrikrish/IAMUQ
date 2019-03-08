library(adaptMCMC)
library(parallel)

source('R/priors.R')
source('R/likelihood.R')

# read in PBS job array index to specify type
aid <- Sys.getenv('PBS_ARRAYID')
# if PBS_ARRAYID doesn't exist, this should be passed as a command line argument or set interactively
if (aid == '') {
  if (!exists('type')) {
    args <- commandArgs(trailingOnly=TRUE)
    type <- args[1]
    expert <- as.logical(args[2])
  }
} else {
  types <- c('iid', 'ar')
  expert <- c(TRUE, FALSE)
  cases <- expand.grid(types=types, expert=expert)
  case <- cases[as.numeric(aid),]
  type <- case[, 'types']
  pi_flag <- case[, 'expert']
}

# read in results of preliminary MCMC run to get starting point and covariance jump matrix
if (pi_flag) {
  filename_prelim <- paste0('mcmc_prelim-pi-', type, '.rds') 
} else {
  filename_prelim <- paste0('mcmc_prelim-nopi-', type, '.rds') 
}

pmcmc_out <- readRDS(file.path('output', filename_prelim))
p0 <- pmcmc_out$samples[nrow(pmcmc_out$samples), ]
cov.jump <- pmcmc_out$cov.jump
parnames <- pmcmc_out$parnames
dat <- pmcmc_out$dat

# set up prior distributions
priors <- create_prior_list(parnames)

# set MCMC parameters
n_iter <- 1e6  # length of MCMC chains
n_chain <- 4 # number of parallel MCMC chains
n_cpu <- detectCores() # size of cluster if necessary

# set up cluster if necessary and run MCMC chain(s)
if (n_cpu > 1) {
  cl <- makeCluster(n_cpu)
  mcmc_out <- MCMC.parallel(log_post, n_iter, init=p0, n.chain=n_chain, n.cpu=n_cpu, scale=cov.jump, adapt=FALSE, list=TRUE, packages=c('mvtnorm'), parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', type), expert=TRUE)
  stopCluster(cl)
} else {
  # run MCMC chain(s)
  mcmc_out <- MCMC(log_post, n_iter, init=p0, scale=cov.jump, adapt=FALSE, list=TRUE, parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', type), expert=TRUE)
}

mcmc_out$dat <- dat
mcmc_out$parnames <- parnames
mcmc_out$type <- type

if (pi_flag) {
  filename_save <- paste0('output/mcmc-pi-', type, '.rds') 
} else {
  filename_save <- paste0('output/mcmc-nopi-', type, '.rds') 
}

saveRDS(mcmc_out, filename_save)
