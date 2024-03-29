library(ggplot2)
library(GGally)

nsamp <- 1e5 # set desired number of samples

mcmc_out <- readRDS('output/mcmc_base-gwp-co2-pop.rds')
mcmc_length <- nrow(mcmc_out[[1]]$samples)
burnin <- 5e5
post <- do.call(rbind, lapply(mcmc_out[1:4], function(l) l$samples[(burnin+1):mcmc_length,]))
parnames <- colnames(post)
# obtain ensemble of posterior samples
idx <- sample(1:nrow(post), nsamp, replace=TRUE)
samps <- post[idx, ]
post_samps <- as.data.frame(samps)

# plot
# set theme
th <- theme_bw(base_size=10) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
theme_set(th)

# create heat map function for the density
heat_dens <- function(data, mapping, colors=heat.colors(12)) {
  ggplot(data = data, mapping = mapping) +
    stat_density_2d(aes(fill = ..level..), geom='polygon') +
    scale_fill_gradientn(colors=colors)
}

# create labeller for facet labels to convert variable code names to mathematical symbols
var_to_sym <- c('psi1' = 'psi[1]',
                'psi2' = 'psi[2]',
                'psi3' = 'psi[3]',
                'P0' = 'P[0]' ,
                'lambda' = 'lambda',
                's' = 's',
                'delta' = 'delta',
                'alpha' = 'alpha',
                'As' = 'A[s]',
                'pi' = 'pi',
                'A0' = 'A[0]',
                'rho2' = 'rho[2]',
                'rho3'= 'rho[3]',
                'tau2' = 'tau[2]',
                'tau3' = 'tau[3]',
                'tau4' = 'tau[4]',
                'kappa' = 'kappa',
                'sigma_pop' = 'sigma[1]' ,
                'sigma_prod' = 'sigma[2]',
                'sigma_emis' = 'sigma[3]',
                'a_11' = 'a[11]',
                'a_22' = 'a[22]',
                'a_33' = 'a[33]',
                'a_21' = 'a[21]',
                'a_31' = 'a[31]',
                'a_12' = 'a[12]',
                'a_13' = 'a[13]',
                'a_23' = 'a[23]',
                'a_32' = 'a[32]',
                'eps_pop' = 'epsilon[1]',
                'eps_prod' = 'epsilon[2]',
                'eps_emis' = 'epsilon[3]'
)

# create pairs plot
p <- ggpairs(post_samps,
      lower = list(continuous=heat_dens),
      upper = list(continuous='cor'),
      diag = list(continuus='barDiag'),
      columnLabels = var_to_sym[colnames(post_samps)],
      labeller = 'label_parsed',
      switch='y') +
     theme(axis.text.x = element_text(angle=90, hjust = 1))

pdf('figures/pairs.pdf', height=28, width=28)
p
dev.off()

png('figures/pairs.png', height=28, width=28, units='in', res=300)
p
dev.off()
