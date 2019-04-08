library(adaptMCMC)
library(parallel)
library(BAUcalib)

rm(list=ls()) # clean up environment

## set case for this run
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
  types <- c('iid', 'var')
  expert <- c(TRUE, FALSE)
  cases <- expand.grid(types=types, expert=expert)
  case <- cases[as.numeric(aid),]
  type <- case[, 'types']
  pi_flag <- case[, 'expert']
}

## set up data
# need to filter the data to the common 1820-2014 range if the type is mvar
dat <- lapply(baudata, function(l) {l[l$year %in% 1820:2014,]})

## set parameter names for appropriate case
all_parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis', 'a_11', 'a_22', 'a_33', 'a_21', 'a_31', 'a_12', 'a_23', 'a_13', 'a_32', 'eps1_pop', 'eps1_prod', 'eps1_emis')

# set parameter name subset
if (type == 'ar') {
  parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis', 'a_pop', 'a_prod', 'a_emis', , 'eps1_pop', 'eps1_prod', 'eps1_emis')
  #, 'eps2_pop', 'eps2_prod', 'eps2_emis')
} else if (type == 'iid') {
  parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis')
} else if (type == 'var') {
  parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis', 'a_11', 'a_22', 'a_33', 'a_21', 'a_31', 'a_12', 'a_23', 'a_13', 'a_32', 'eps1_pop', 'eps1_prod', 'eps1_emis')
}

## set up prior distributions
priors <- create_prior_list(parnames)

########## set up cluster ########################################
# start cluster
print('Starting cluster...')
ncores <- detectCores()
cl <- makeCluster(ncores)

######### estimate MAP as starting point for MCMC #################
print('Estimating MAP...')

## set DEoptim parameters for MAP estimate
# set upper and lower bounds for each variable
lbound <- c(0.0001, 0, 6.9, 0.3, 0.5, 0.1, 0.01, 0.0007, 5.3, 0.4, 0, 0, 0, 1700, 1700, 2010, 0.005, 0, 0, 0, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0, 0, 0)
ubound <- c(0.15, 100, 15, 1, 0.9,  0.3, 0.14, 0.0212, 16.11, 0.7, 3, 0.5, 0.5, 2100, 2100, 2500, 0.2, 5, 5, 5, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 5, 5, 5)
# if this is run outside of batch mode, print DEoptim progress
if (aid == '') {
  trace_flag <- TRUE
} else {
  trace_flag <- FALSE
}

map <- DEoptim(neg_log_post, lbound[match(parnames, all_parnames)], ubound[match(parnames, all_parnames)], control=list(NP=10*length(parnames), itermax=3000, parallelType=1, cluster=cl, packages=c('BAUcalib', 'mvtnorm'), parVar=c('parnames', 'dat', 'priors', paste0('log_lik_', type))), parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', type), exp_gwp=pi_flag)$optim$bestmem
names(map) <- parnames
print(map)
print('Done.')

############# run MCMC chains ######################################
print('Running MCMC...')
# set MCMC parameters
n_iter <- 5e6  # length of MCMC chains
n_chain <- 4 # number of parallel MCMC chains

# set initial step size to be 10% of the standard deviation of the prior distributions when defined, or else some initial guess otherwise
stepsize <- numeric(length(parnames))
for (name in parnames) {
  if (priors[[name]][['type']] == 'uniform') {
    stepsize[match(name, parnames)] <- 0.1*(priors[[name]][['max']] - priors[[name]][['min']])/sqrt(12)
  } else if (priors[[name]][['type']] == 'normal') {
    stepsize[match(name, parnames)] <- 0.1*priors[[name]][['sd']]
  } else if (priors[[name]][['type']] == 'log-normal') {
    stepsize[match(name, parnames)] <- 0.1*sqrt((exp(priors[[name]][['sdlog']]^2)-1)*exp(2*priors[[name]][['meanlog']]+priors[[name]][['sdlog']]^2))
  }
}

adapt_start <- max(500, round(0.01*n_iter)) # when to start stepsize adaptation

# run MCMC chain(s)
mcmc_out <- MCMC.parallel(log_post, n_iter, init=map, n.chain=n_chain, n.cpu=ncores, scale=stepsize, adapt=2e6, list=TRUE, packages=c('mvtnorm', 'BAUcalib'), parnames=parnames, priors=priors, dat=dat, lik_fun=paste0('log_lik_', type), exp_gwp=pi_flag)

print('Done.')

########### save output ########################
print('Saving output...')

if (pi_flag) {
  filename_save <- paste0('output/mcmc-pi-', type, '.rds') 
} else {
  filename_save <- paste0('output/mcmc-nopi-', type, '.rds') 
}
saveRDS(mcmc_out, filename_save)
print('Done.')

######## stop cluster
print('Stopping cluster...')
stopCluster(cl)
print('Done.')