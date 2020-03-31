library(reshape2)
library(ggplot2)
library(IAMUQ)
source('R/calib_priors.R')

nsamp <- 1e5 # set desired number of samples

appendices <- c('', '-gwp', '-co2', '-gwp-co2')
app_names <- c('None', 'Christensen et al (2018)', 'Ho et al (2019)', 'Both')

scenarios <- c('base', 'low', 'high', 'alt_zc')

for (i in 1:length(scenarios)) {
  post_samps <- vector('list', length(appendices))
  names(post_samps) <- app_names
  
  for (app in appendices) {
    # read MCMC file
    mcmc_out <- readRDS(paste0('output/mcmc_', scenarios[i],  app, '.rds'))
  
    mcmc_length <- nrow(mcmc_out[[1]]$samples)
    burnin <- 5e5
    post <- do.call(rbind, lapply(mcmc_out[1:4], function(l) l$samples[(burnin+1):mcmc_length,]))
    parnames <- colnames(post)
    # obtain ensemble of posterior samples
    idx <- sample(1:nrow(post), nsamp, replace=TRUE)
    samps <- post[idx, ]
    post_samps[[match(app, appendices)]] <- as.data.frame(samps)
  }

  # combine samples into a list for melting
  all_melt <- melt(post_samps)
  colnames(all_melt)[3] <- 'Inversion'
  colnames(all_melt)[1] <- 'Variable'

# plot
# set theme
th <- theme_bw(base_size=10) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
theme_set(th)

# create labeller for facet labels to convert variable code names to mathematical symbols
var_to_sym <- c('psi1' = expression(psi[1]),
                'psi2' = expression(psi[2]),
                'psi3' = expression(psi[3]),
                'P0' = expression(P[0]) ,
                'lambda' = expression(lambda),
                's' = expression(s),
                'delta' = expression(delta),
                'alpha' = expression(alpha),
                'As' = expression(A[s]),
                'pi' = expression(pi),
                'A0' = expression(A[0]),
                'rho2' = expression(rho[2]),
                'rho3'= expression(rho[3]),
                'tau2' = expression(tau[2]),
                'tau3' = expression(tau[3]),
                'tau4' = expression(tau[4]),
                'kappa' = expression(kappa),
                'sigma_pop' = expression(sigma[1]) ,
                'sigma_prod' = expression(sigma[2]) ,
                'sigma_emis' = expression(sigma[3]),
                'a_11' = expression(a[11]),
                'a_22' = expression(a[22]),
                'a_33' = expression(a[33]),
                'a_21' = expression(a[21]),
                'a_31' = expression(a[31]),
                'a_12' = expression(a[12]),
                'a_13' = expression(a[13]),
                'a_23' = expression(a[23]),
                'a_32' = expression(a[32]),
                'eps_pop' = expression(epsilon[1]),
                'eps_prod' = expression(epsilon[2]),
                'eps_emis' = expression(epsilon[3])
)

  levels(all_melt$Variable) <- var_to_sym[levels(all_melt$Variable)]
  all_melt$Inversion <- factor(all_melt$Inversion, levels=app_names)
 
  p <- ggplot(all_melt) + stat_density(aes(x=value, color=Inversion), geom='line', position='identity') + facet_wrap(vars(Variable), scales='free', labeller=label_parsed, ncol=4) + scale_color_brewer(palette='Dark2') + theme(legend.position='bottom') + scale_y_continuous('Density') + scale_x_continuous('Parameter Value')

  pdf(paste0('figures/expert-dist-', scenarios[i], '.pdf'), height=7, width=7)
  print(p)
  dev.off()

  png(paste0('figures/expert-dist-', scenarios[i], '.png'), height=7, width=7, res=600, units='in')
  print(p)
  dev.off()
}