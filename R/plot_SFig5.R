library(ggplot2)
library(reshape2)

cases <- c('base', 'low', 'high', 'alt_zc')
case_labels <- c('Base', 'Low Fossil Fuel', 'High Fossil Fuel', 'Delayed Zero-Carbon')

n_samp <- 1e5

sample_param <- function(case, parname, n) {
  mcmc_out <- readRDS(paste0('output/mcmc_', case, '.rds'))
  mcmc_length <- nrow(mcmc_out[[1]]$samples)
  burnin <- 5e5
  post <- do.call(rbind, lapply(mcmc_out[1:4], function(l) l$samples[(burnin+1):mcmc_length,]))
  parnames <- colnames(post)
  # obtain ensemble of posterior samples
  idx <- sample(1:nrow(post), nsamp, replace=TRUE)
  post[idx, parname]
}

# sample from posterior for each case
samps <- lapply(cases, sample_param, parname='tau4', n=n_samp)
# melt samples
samp_melt <- melt(samps)
samp_melt[, 'Scenario'] <- factor(cases[samp_melt$L1], levels=cases, labels=case_labels)

cbbpsqualitative <- c("#000000", "#e79f00", "#9ad0f3", "#CC79A7", "#0072B2", "#009E73", "#F0E442", "#D55E00")

# plot distributions
# set theme
th <- theme_bw(base_size=10) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
theme_set(th)

p <- ggplot(samp_melt) + stat_density(aes(x=value, group=Scenario, fill=Scenario), alpha=0.3, geom='area', position='identity') + scale_x_continuous(expression(tau[4]~(year)), expand=c(0.001, 0.001)) + scale_y_continuous(expand=c(0.001, 0.001)) + scale_fill_manual('Model Scenario', values=cbbpsqualitative) + theme(legend.position='bottom', legend.margin=margin(l=-0.75, b=0.5, unit='cm'), plot.margin=margin(l=0.1, r=0.12, t=0.1, unit='in')) + guides(fill = guide_legend(nrow=2, byrow=FALSE))

pdf('figures/tau4-dist.pdf', height=4, width=4)
p
dev.off()

png('figures/tau4-dist.png', height=4, width=4, units='in', res=600)
p
dev.off()
