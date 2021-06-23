library(parallel)
library(IAMUQ)

set.seed(1024)

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

appendix <- ''
if (exp_gwp) {
  appendix <- paste0(appendix, '-gwp')
}
if (exp_co2) {
  appendix <- paste0(appendix, '-co2')
}

nsamp <- 1e5 # number of simulations
yrs <- 2015:2200 # years which the conditional simulation should project

# read in data
dat <- lapply(iamdata, function(l) {l[l$year %in% 1820:2014,]})

# start cluster
cl <- makeCluster(detectCores())
# export packages, functions, and objects
clusterEvalQ(cl, library(IAMUQ))
clusterExport(cl, c('dat', 'yrs'))

# conditionally simulate to 2100
# get posterior samples
mcmc_out <- readRDS(paste0('output/mcmc_', scenario, appendix, '.rds'))
mcmc_length <- nrow(mcmc_out[[1]]$samples)
burnin <- 5e5
post <- do.call(rbind, lapply(mcmc_out[1:4], function(l) l$samples[(burnin+1):mcmc_length,]))
parnames <- colnames(post)
# obtain ensemble of posterior samples
idx <- sample(1:nrow(post), nsamp, replace=TRUE)
samps <- post[idx, ]
clusterExport(cl, c('parnames'))
sim_out <- parApply(cl, samps, 1, cond_sim_model, parnames=parnames, dat=dat, projyrs=yrs)
sim_out <- lapply(1:length(sim_out), function(i) list(out=sim_out[[i]], pars=samps[i,]))

# save simulation output
## save estimate
saveRDS(sim_out, paste0('output/sim_', scenario, appendix, '.rds'))

stopCluster(cl)

