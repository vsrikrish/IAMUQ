rm(list=ls()) # clean up environment

library(IAMUQ)

source('R/calib_priors.R')

## set case for this run
# read in PBS job array index to specify type
aid <- Sys.getenv('PBS_ARRAYID')
# if PBS_ARRAYID doesn't exist, this should be passed as a command line argument or set interactively
if (aid == '') {
  if (!exists('type')) {
    args <- commandArgs(trailingOnly=TRUE)
    scenario <- args[1]
    exp_assess <- args[2]
  }
} else {
  scenarios <- c('iid', 'base', 'short', 'low', 'high', 'del_zc')
  exp_assess <- c('none', 'gwp', 'co2', 'both')
  cases <- expand.grid(scenarios=scenarios, exp=exp_assess)
  id <- as.numeric(aid)
  scenario <- cases[id, 'scenarios']
  exp_assess <- cases[id, 'exp'] 
}

exp_gwp <- FALSE
exp_co2 <- FALSE

if (exp_assess == 'gwp') {
  exp_gwp <- TRUE
} else if (exp_assess == 'co2') {
  exp_co2 <- TRUE
} else if (exp_assess == 'both') {
  exp_gwp <- TRUE
  exp_co2 <- TRUE
}

## set model run parameters associated with each scenario
# the base scenario corresponds to the model defaults, but we set it up here for consistency's sake
if (scenario == 'base') {
  data_yrs <- 1820:2014 # years for observational constraints
  ff_thresh <- 6000 # fossil fuel constraint in GtC
  ff_const_yrs <- 1700:2500 # years over which fossil fuel constraint is evaluated
  residtype <- 'var'
} else if (scenario == 'short') {
  data_yrs <- 1950:2014 # years for observational constraints
  ff_thresh <- 6000 # fossil fuel constraint in GtC
  ff_const_yrs <- 1700:2500 # years over which fossil fuel constraint is evaluated
  residtype <- 'var' # residual structure type
} else if (scenario == 'iid') {
  data_yrs <- 1820:2014 # years for observational constraints
  ff_thresh <- 6000 # fossil fuel constraint in GtC
  ff_const_yrs <- 1700:2500 # years over which fossil fuel constraint is evaluated
  residtype <- 'iid'
} else if (scenario == 'low') {
  data_yrs <- 1820:2014 # years for observational constraints
  ff_thresh <- 3000 # fossil fuel constraint in GtC
  ff_const_yrs <- 2015:2500 # years over which fossil fuel constraint is evaluated
  residtype <- 'var' # residual structure type
} else if (scenario == 'high') {
  data_yrs <- 1820:2014 # years for observational constraints
  ff_thresh <- 10000 # fossil fuel constraint in GtC
  ff_const_yrs <- 2015:2500 # years over which fossil fuel constraint is evaluated
  residtype <- 'var' # residual structure type
} else if (scenario == 'del_zc') {
  data_yrs <- 1820:2014 # years for observational constraints
  ff_thresh <- 6000 # fossil fuel constraint in GtC
  ff_const_yrs <- 1700:2500 # years over which fossil fuel constraint is evaluated
  residtype <- 'var' # residual structure type
}

if (residtype == 'ar') {
  parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis', 'a_pop', 'a_prod', 'a_emis', , 'eps1_pop', 'eps1_prod', 'eps1_emis')
} else if (residtype == 'iid') {
  parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis')
} else if (residtype == 'var') {
  parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis', 'a_11', 'a_22', 'a_33', 'a_21', 'a_31', 'a_12', 'a_23', 'a_13', 'a_32', 'eps_pop', 'eps_prod', 'eps_emis')
}

## set up prior dataframe
prior_df <- set_prior_params(parnames)
# if scenario involves changing prior distributions, do so here
if (scenario == 'del_zc') {
  # modify zero-carbon half-saturation year (tau_4) prior distribution
  zc_idx <- match('tau4', prior_df[, 'name'])
  prior_df[zc_idx, 'type'] <- 'truncnorm'
  prior_df[zc_idx, 'lower'] <- 2100
  prior_df[zc_idx, 'upper'] <- 2400
}

## read in MAP estimate for that scenario as the initial value as the
## save estimate
appendix <- ''
if (exp_gwp) {
  appendix <- paste0(appendix, '-gwp')
}
if (exp_co2) {
  appendix <- paste0(appendix, '-co2')
}
map <- readRDS(paste0('output/map_', scenario, appendix, '.rds'))

## run MCMC chains (using parallel default)
mcmc_out <- run_mcmc(log_post, parnames=parnames, residtype=residtype, prior_df=prior_df, data_yrs=data_yrs, init=map$optim$bestmem, n_iter=2e6, thresh=ff_thresh, ff_const_yrs=ff_const_yrs, exp_gwp=exp_gwp, exp_co2=exp_co2)

## save estimate
saveRDS(mcmc_out, paste0('output/mcmc_', scenario, appendix, '.rds'))
