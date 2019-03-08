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

# read in MLE estimate to start MCMC chain as well as data
mle_out <- readRDS(paste0('output/mle-', type, '.rds'))
mle <- mle_out$mle$optim$bestmem
parnames <- mle_out$parnames
dat <- mle_out$dat

# set up prior distributions
priors <- create_prior_list(parnames)

# set MCMC parameters for initial run
n_iter <- 1e5  # length of MCMC chains
n_chain <- 1 # number of parallel MCMC chains
n_cpu <- 1 # size of cluster if necessary
rate_accept_many <- 0.234 # optimal acceptance rate for multiple parameters
rate_accept_one <- 0.44 # optimal acceptance rate for one parameter
#rate_accept <- rate_accept_many + (rate_accept_one - rate_accept_many)/length(parnames)
rate_accept <- 0.3

# set initial step size to be 0.5% of the standard deviation of the prior distributions when defined, or else some initial guess otherwise
stepsize <- numeric(length(parnames))
for (name in parnames) {
  if (priors[[name]][['type']] == 'uniform') {
    stepsize[match(name, parnames)] <- 0.01*(priors[[name]][['max']] - priors[[name]][['min']])/sqrt(12)
  } else if (priors[[name]][['type']] == 'normal') {
    stepsize[match(name, parnames)] <- 0.01*priors[[name]][['sd']]
  } else if (priors[[name]][['type']] == 'log-normal') {
    stepsize[match(name, parnames)] <- 0.01*sqrt((exp(priors[[name]][['sdlog']]^2)-1)*exp(2*priors[[name]][['meanlog']]+priors[[name]][['sdlog']]^2))
  }
}

adapt_start <- max(1000, round(0.05*n_iter)) # when to start stepsize adaptation

# set up cluster if necessary and run MCMC chain(s)
if (n_cpu > 1) {
  mcmc_out <- MCMC.parallel(log_post, n_iter, init=mle, n.chain=n_chain, n.cpu=n_cpu, scale=stepsize, gamma=0.51, list=TRUE, n.start=adapt.start, acc.rate=rate_accept, packages=c('mvtnorm'), parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', type), expert=pi_flag)
} else {
  # run MCMC chain(s)
  mcmc_out <- MCMC(log_post, n_iter, init=mle, scale=stepsize, gamma=0.51, list=TRUE, n.start=adapt_start, acc.rate=rate_accept, parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', type), expert=pi_flag)
}

mcmc_out$dat <- dat
mcmc_out$parnames <- parnames
mcmc_out$mle <- mle
mcmc_out$type <- type

if (pi_flag) {
  filename_save <- paste0('output/log/mcmc_prelim-pi-', type, '.rds') 
} else {
  filename_save <- paste0('output/log/mcmc_prelim-nopi-', type, '.rds') 
}

saveRDS(mcmc_out, filename_save)
