library(IAMUQ)
library(parallel)

source('R/calib_priors.R')

# use PBS array id to set the seed for partitioning the data
aid <- Sys.getenv('PBS_ARRAYID')
set.seed(as.numeric(aid))

# load data
dat <- lapply(iamdata, function(l) l[l$year %in% 1820:2014, ])

# partition data
hoyrs <- sample(dat[[1]]$year, size=floor(length(dat[[1]]$year)/5), replace=FALSE)

## run MCMC if no MCMC file exists
mcmc_fname <- paste0('~/scratch/crossval/mcmc-', aid, '.rds')
if (file.exists(mcmc_fname)) {
  mcmc_out <- readRDS(mcmc_fname)
} else {
  # read in full data MAP estimate as the MCMC initial value
  map <- readRDS('output/map_base-gwp-co2.rds')$optim$bestmem
  parnames <- names(map)
  
  # get prior dataframe
  prior_df <- set_prior_params(parnames)
  
  # run MCMC
  mcmc_out <- run_mcmc(log_post, parnames=parnames, residtype='var', prior_df=prior_df, data_yrs=1820:2014, init=map, n_iter=2e6, exp_gwp=TRUE, hoyrs=hoyrs)
  
  saveRDS(mcmc_out, paste0('~/scratch/crossval/mcmc-', aid, '.rds'))
}

# combine MCMC chains and thin MCMC chain for conditional simulation
nsamp <- 1e5
mcmc_length <- nrow(mcmc_out[[1]]$samples)
burnin <- 5e5
post <- do.call(rbind, lapply(mcmc_out[1:4], function(l) l$samples[(burnin+1):mcmc_length,]))
parnames <- colnames(post)
# obtain ensemble of posterior samples
idx <- sample(1:nrow(post), nsamp, replace=TRUE)
samps <- post[idx, ]

# conditionally simulate
# start cluster and simulate
cl <- makeCluster(detectCores())
clusterEvalQ(cl, library(IAMUQ))
clusterEvalQ(cl, library(magic))
clusterExport(cl, c('parnames', 'dat'))
sim_out <- parApply(cl, samps, 1, cond_sim_model, parnames=parnames, dat=dat, hoyrs=hoyrs)
stopCluster(cl)

# save conditional simulation output
saveRDS(sim_out, paste0('~/scratch/crossval/sim-', aid, '.rds'))
