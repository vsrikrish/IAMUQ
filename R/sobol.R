library(IAMUQ)
library(sensitivity)
library(readr)

source('R/sobol_functions.R')

n_samp <- 1e5 # define size of ensemble
n_boot <- 1e4 # number of bootstrap samples

## load and combine MCMC output
mcmc_out <- readRDS('output/mcmc_base.rds')
mcmc_length <- nrow(mcmc_out[[1]]$samples)
burnin <- 5e5
post <- do.call(rbind, lapply(mcmc_out[1:4], function(l) l$samples[(burnin+1):mcmc_length,]))
parnames <- colnames(post)[1:17]
names(parnames) <- parnames

## set up Sobol sample ensembles
# get max and min values for each posterior distribution's samples
bds <- lapply(parnames, function(p) {c(min(post[, p]), max(post[, p]))})
# fit KDEs to each marginal posterior
dens <- lapply(parnames, function(p) density(post[, p], from=bds[[p]][1], to=bds[[p]][2], kernel='gaussian'))
# sample MCMC indices
sob_samp1_idx <- sample(1:nrow(post), n_samp, replace=TRUE)
sob_samp2_idx <- sample(1:nrow(post), n_samp, replace=TRUE)
# sample values for each parameter
sob_samp1 <- lapply(parnames, function(p) sample_value(n=n_samp, parvals=post[sob_samp1_idx, p], bw=dens[[p]]$bw, bds=bds[[p]]))
sob_samp2 <- lapply(parnames, function(p) sample_value(n=n_samp, parvals=post[sob_samp2_idx, p], bw=dens[[p]]$bw, bds=bds[[p]]))
# convert sample lists to data frame for sobolSalt function
sob_samp1 <- as.data.frame(do.call(cbind, sob_samp1))
sob_samp2 <- as.data.frame(do.call(cbind, sob_samp2))

# sample TCRE
# we fit a log-normal distribution to the 90% range from Gillet et al (2013)
low <- 0.7
high <- 2
tcre_mean <- mean(c(log(low), log(high)))
tcre_sd <- (log(high) - log(low)) / (qnorm(0.95) - qnorm(0.05))
bds[['TCRE']] <- c(0, 6)
sob_samp1[, 'TCRE'] <- map_range(rlnorm(n_samp, meanlog=tcre_mean, sdlog=tcre_sd), bdin=c(0, 6), bdout=c(0, 1))
sob_samp2[, 'TCRE'] <- map_range(rlnorm(n_samp, meanlog=tcre_mean, sdlog=tcre_sd), bdin=c(0, 6), bdout=c(0, 1))
parnames['TCRE'] <- 'TCRE'


## calculate temperature anomalies using 1861-1880 baseline from HadCRUT4
# read in HadCRUT4 data file
temp_dat <- read_table('data/HadCRUT4-gl.dat', col_names=FALSE)
temp_dat <- temp_dat[seq(1, nrow(temp_dat), by=2), c(1, ncol(temp_dat))]
colnames(temp_dat) <- c('year', 'temp')
# compute 1861-1880 mean
temp_mean <- mean(temp_dat[temp_dat$year %in% 1861:1880,]$temp)
# compute anomaly as of 2014
temp_base <- temp_dat[temp_dat$year == 2014,]$temp - temp_mean

## run the Sobol analysis
sens_out <- sobolSalt(model=temp_eval_par, X1=sob_samp1, X2=sob_samp2, scheme='B', nboot=n_boot, parnames=parnames, par_bds=bds, baseline=temp_base)

# write output file as in Tony and Perry's analysis codes
sobolout1 <- paste0('output/Sobol-1-tot_temp.txt')
sobolout2 <- paste0('output/Sobol-2-tot_temp.txt')

headers.1st.tot <- matrix(c('Parameter', 'S1', 'S1_conf_low', 'S1_conf_high',
                            'ST', 'ST_conf_low', 'ST_conf_high'), nrow=1)
output.1st.tot  <- data.frame(cbind( parnames,
                                     sens_out$S[,1],
                                     sens_out$S[,4],
                                     sens_out$S[,5],
                                     sens_out$T[,1],
                                     sens_out$T[,4],
                                     sens_out$T[,5]))
write.table(headers.1st.tot, file=sobolout1, append=FALSE, sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
write.table(output.1st.tot , file=sobolout1, append=TRUE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
headers.2nd     <- matrix(c('Parameter_1', 'Parameter_2', 'S2', 'S2_conf_low',         'S2_conf_high'), nrow=1)
output2.indices <- sens_out$S2[,1]
output2.conf1   <- sens_out$S2[,4]
output2.conf2   <- sens_out$S2[,5]
# 2nd order index names ordered as: (assuming 39 parameters)
# 1. parnames[1]-parnames[2]
# 2. parnames[1]-parnames[3]
# 3. parnames[1]-parnames[4]
# ... etc ...
# 38. parnames[1]-parnames[39] << N=2:39 => p1-p[N]
# 39. parnames[2]-parnames[3]
# 40. parnames[2]-parnames[4]
# 38+37. parnames[2]-parnames[39] << N=3:39 => p2-p[N]
# ... etc ...
names2  <- rownames(sens_out$S2)
names2a <- rep(NA, length(names2))
names2b <- rep(NA, length(names2))
cnt <- 1
for (i in seq(from=1, to=(length(parnames)-1), by=1)) {           # i = index of first name
    for (j in seq(from=(i+1), to=(length(parnames)), by=1)) {   # j = index of second name
        names2a[cnt] <- parnames[i]
        names2b[cnt] <- parnames[j]
        cnt <- cnt+1
    }
}

output.2nd <- data.frame(cbind( names2a,
                                names2b,
                                output2.indices,
                                output2.conf1,
                                output2.conf2 ))
write.table(headers.2nd    , file=sobolout2, append=FALSE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
write.table(output.2nd     , file=sobolout2, append=TRUE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)

