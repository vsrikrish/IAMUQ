library(DEoptim)
library(parallel)

R_path <- 'R'
data_path <- 'data'
output_path <- 'output'

# read in PBS job array index to specify type
types <- c('iid', 'ar', 'mvar', 'ar2')
aid <- Sys.getenv('PBS_ARRAYID')
# if PBS_ARRAYID doesn't exist, this should be passed as a command line argument or set interactively
if (aid == '') {
   if (!exists('type')) {
    args <- commandArgs(trailingOnly=TRUE)
    type <- args[1]
   }
} else {
  type <- types[as.numeric(aid)]
}

# read in calibration data. if this has been done already and saved,
# read the processed data; otherwise, process the raw data.
if (file.exists(file.path(data_path, 'calib_data.rds'))) {
  dat <- readRDS(file.path(data_path,'calib_data.rds'))
} else {
  source(file.path(R_path, 'read_data.R'))
}

# need to filter the data to the common 1820-2014 range if the type is mvar
dat <- lapply(dat, function(l) {l[l$year %in% 1820:2014,]})

all_parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis', 'rho_pop', 'rho_prod', 'rho_emis',  'eps1_pop', 'eps1_prod', 'eps1_emis', 'eps2_pop', 'eps2_prod', 'eps2_emis')

# set parameter name subset
if (type == 'ar') {
  parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'rho_pop', 'rho_prod', 'rho_emis', 'sigma_pop', 'sigma_prod', 'sigma_emis', 'eps1_pop', 'eps1_prod', 'eps1_emis')
  #, 'eps2_pop', 'eps2_prod', 'eps2_emis')
} else if (type == 'iid') {
  parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'sigma_pop', 'sigma_prod', 'sigma_emis')
} else if (type == 'mvar') {
    parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'rho_pop', 'rho_prod', 'rho_emis')
} else if (type == 'ar2') {
  parnames <- c('psi1', 'psi2', 'psi3', 'P0', 'lambda', 's', 'delta', 'alpha', 'As', 'pi', 'A0', 'rho2', 'rho3', 'tau2', 'tau3', 'tau4', 'kappa', 'rho_pop', 'rho_prod', 'rho_emis', 'sigma_pop', 'sigma_prod', 'sigma_emis')
}


# estimate maximum-likelihood parameter values
# start cluster
print('Starting cluster...')
ncores <- detectCores()
print(ncores)
cl <- makeCluster(ncores)

# source model and likelihood function files
source(file.path(R_path, 'model.R'))
source(file.path(R_path, 'likelihood.R'))

lbound <- c(0.0001, 0, 6.9, 0.3, 0.5, 0.1, 0.01, 0.0007, 5.3, 0.4, 0, 0, 0, 1700, 1700, 2010, 0.005, 0, 0, 0, 0.4, 0.4, 0.4, 0, 0, 0, 0, 0, 0)
ubound <- c(0.15, 100, 15, 1, 0.9,  0.3, 0.14, 0.0212, 16.11, 0.7, 3, 0.5, 0.5, 2100, 2100, 2500, 0.2, 5, 5, 5, 2, 2, 2, 5, 5, 5, 5, 5, 5)
if (aid == '') {
   trace_flag <- TRUE
} else {
  trace_flag <- FALSE
}
mle <- DEoptim(neg_log_lik, lbound[match(parnames, all_parnames)], ubound[match(parnames, all_parnames)], control=list(NP=10*length(parnames), itermax=5000, F=0.6, CR=0.42, trace=trace_flag, parallelType=1, cluster=cl, packages=c('mvtnorm'), parVar=c('parnames', 'dat', 'update_pop', 'update_prod', 'update_emis', 'mod', 'check_constraints', 'residuals', 'ar_residuals',paste0('log_lik_', type), 'dmvnorm')), parnames=parnames, dat=dat, lik_fun=paste0('log_lik_', type))

model_out <- list()
model_out[['bestfit']] <- mod(mle$optim$bestmem, parnames)
model_out[['mle']] <- mle
model_out[['dat']] <- dat
model_out[['parnames']] <- parnames
model_out[['type']] <- type

saveRDS(model_out, file.path(output_path, 'log', paste('mle-', type, '.rds', sep='')))

# if (type == 'iid') {
#   plot_type <- 'iid'
# } else if (grepl('ar', type)) {
#   plot_type <- 'ar'
# }
# source(file.path(R_path, paste('plot_', plot_type, '_residuals.R', sep='')))

stopCluster(cl)
print('Stopping cluster...')
