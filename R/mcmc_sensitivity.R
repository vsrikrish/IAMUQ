rm(list=ls()) # clean up environment

library(IAMUQ)

source('R/calib_priors.R')
source('R/compute_fossil_threshold.R')

## set case for this run
# read in PBS job array index to specify type
aid <- Sys.getenv('PBS_ARRAYID')
# if PBS_ARRAYID doesn't exist, this should be passed as a command line argument or set interactively
if (aid == '') {
  if (!exists('type')) {
    args <- commandArgs(trailingOnly=TRUE)
    scenario <- args[1]
  }
} else {
  scenarios <- c('alt_s', 'alt_lambda', 'alt_pi', 'alt_As', 'alt_multi', 'alt_zc')
  id <- as.numeric(aid)
  scenario <- scenarios[id]
}

## set model run parameters
data_yrs <- 1820:2019 # years for observational constraints
ff_thresh <- compute_fossil_threshold('base') # fossil fuel constraint in GtC
ff_const_yrs <- 2012:2500 # years over which fossil fuel constraint is evaluated
exp_gwp <- TRUE # do we do probabilistic inversion for average GWP per capita?
exp_co2 <- TRUE
exp_pop <- TRUE
residtype <- 'var'
## set fossil fuel penetration windows for coal and renewable penetration
## based on data from BP Statistical Review of World Energy
ff_pen_yr <- 2019
ff_pen_window <- list(cbind(c(20, 30), c(NA, NA), c(10, 20))
  
if (residtype == 'ar') {
  parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis', 'a_pop', 'a_prod', 'a_emis', , 'eps1_pop', 'eps1_prod', 'eps1_emis')
} else if (residtype == 'iid') {
  parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis')
} else if (residtype == 'var') {
  parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis', 'a_11', 'a_22', 'a_33', 'a_21', 'a_31', 'a_12', 'a_23', 'a_13', 'a_32', 'eps_pop', 'eps_prod', 'eps_emis')
}

## set up prior dataframe
prior_df <- set_prior_params(parnames)

## change prior distribution based on scenario
if (scenario == 'alt_s') {
  idx <- match('s', prior_df[, 'name'])
  prior_df[idx, 'type'] <- 'log-normal'
} else if (scenario == 'alt_lambda') {
  idx <- match('lambda', prior_df[, 'name'])
  prior_df[idx, 'type'] <- 'log-normal'
} else if (scenario == 'alt_pi') {
  idx <- match('pi', prior_df[, 'name'])
  prior_df[idx, 'type'] <- 'log-normal'
} else if (scenario == 'alt_As') {
  idx <- match('As', prior_df[, 'name'])
  prior_df[idx, 'type'] <- 'normal'
} else if (scenario == 'alt_multi') {
  idx <- match('lambda', prior_df[, 'name'])
  prior_df[idx, 'type'] <- 'log-normal'
  idx <- match('s', prior_df[, 'name'])
  prior_df[idx, 'type'] <- 'log-normal'
  idx <- match('pi', prior_df[, 'name'])
  prior_df[idx, 'type'] <- 'log-normal'
  idx <- match('As', prior_df[, 'name'])
  prior_df[idx, 'type'] <- 'normal'
} else if (scenario == 'alt_zc') {
  idx <- match('tau4', prior_df[, 'name'])
  prior_df[idx, 'upper'] <- 2250
}



## read in MAP estimate for that scenario as the initial value as the
map <- readRDS(paste0('output/map_', scenario, '-gwp-co2-pop.rds'))

## run MCMC chains (using parallel default)
mcmc_out <- run_mcmc(log_post, parnames=parnames, residtype=residtype, prior_df=prior_df, data_yrs=data_yrs, init=map$optim$bestmem, n_iter=2e6, ff_thresh=ff_thresh, ff_const_yrs=ff_const_yrs, ff_pen_windows=ff_pen_window, ff_pen_yrs=ff_pen_yr, exp_gwp=exp_gwp, exp_co2=exp_co2, exp_pop=exp_pop)

## save estimate
saveRDS(mcmc_out, paste0('output/mcmc_', scenario, '-gwp-co2-pop.rds'))
