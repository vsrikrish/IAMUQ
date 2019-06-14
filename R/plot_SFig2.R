library(reshape2)
library(ggplot2)
library(IAMUQ)
source('R/calib_priors.R')

nsamp <- 1e5 # set desired number of samples

mcmc_out <- readRDS('output/mcmc_base.rds')
mcmc_length <- nrow(mcmc_out[[1]]$samples)
burnin <- 5e5
post <- do.call(rbind, lapply(mcmc_out[1:4], function(l) l$samples[(burnin+1):mcmc_length,]))
parnames <- colnames(post)
# obtain ensemble of posterior samples
idx <- sample(1:nrow(post), nsamp, replace=TRUE)
samps <- post[idx, ]
post_samps <- as.data.frame(samps)

# set up prior list
prior_df <- set_prior_params(parnames)
priors <- create_prior_list(prior_df)
# sample from priors
pri_samps <- as.data.frame(do.call(cbind, lapply(priors,
    function(p) do.call(match.fun(p[['rand.fun']]),
                  c(list(n=nsamp),
                  p[-which(names(p) %in% c('type', 'dens.fun', 'quant.fun', 'rand.fun'))]))
)))

# combine samples into a list for melting
all_samps <- list('Posterior'=post_samps, 'Prior'=pri_samps)
all_melt <- melt(all_samps)
colnames(all_melt)[3] <- 'Distribution'
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

p <- ggplot(all_melt) + stat_density(aes(x=value, color=Distribution), geom='line', position='identity') + facet_wrap(vars(Variable), scales='free', labeller=label_parsed, ncol=4) + scale_color_brewer(palette='Dark2') + theme(legend.position='bottom') + scale_y_continuous('Density') + scale_x_continuous('Parameter Value')

pdf('figures/FigS2-dist.pdf', height=8, width=8)
p
dev.off()

png('figures/FigS2-dist.png', height=8, width=8, res=600, units='in')
p
dev.off()
