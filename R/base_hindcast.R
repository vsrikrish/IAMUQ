library(parallel)
library(IAMUQ)

exp_gwp <- TRUE
exp_co2 <- TRUE
exp_pop <- TRUE

appendix <- ''
if (exp_gwp) {
  appendix <- paste0(appendix, '-gwp')
}
if (exp_co2) {
  appendix <- paste0(appendix, '-co2')
}
if (exp_pop) {
  appendix <- paste0(appendix, '-pop')
}

nsamp <- 1e5 # number of simulations
yrs <- 1900:2019 # years which the conditional simulation should project

# read in data
dat <- lapply(iamdata, function(l) {l[l$year %in% 1820:2019,]})

# start cluster
cl <- makeCluster(detectCores())
# export packages, functions, and objects
clusterEvalQ(cl, library(IAMUQ))
clusterExport(cl, c('dat', 'yrs'))

# conditionally simulate to 2100
# get posterior samples
mcmc_out <- readRDS(paste0('output/mcmc_base',appendix, '.rds'))
mcmc_length <- nrow(mcmc_out[[1]]$samples)
burnin <- 5e5
post <- do.call(rbind, lapply(mcmc_out[1:4], function(l) l$samples[(burnin+1):mcmc_length,]))
parnames <- colnames(post)
# obtain ensemble of posterior samples
idx <- sample(1:nrow(post), nsamp, replace=TRUE)
samps <- post[idx, ]
clusterExport(cl, c('parnames'))
sim_out <- parApply(cl, samps, 1, cond_sim_model, parnames=parnames, dat=dat, hoyrs=yrs)
sim_out <- lapply(1:length(sim_out), function(i) list(out=sim_out[[i]], pars=samps[i,]))

# save simulation output
## save estimate
saveRDS(sim_out, paste0('output/hindcast_base', appendix, '.rds'))

stopCluster(cl)

