library(ggplot2)
library(reshape2)
library(IAMUQ)
library(grid)
library(gridExtra)
library(gtable)

nsamp <- 1e5

yrs <- 1900:2019
dat <- lapply(iamdata, function(l) l[l$year %in% 1820:2019, ])

sim_out <- readRDS('output/hindcast_base-gwp-co2-pop.rds')

mod_out_names <- c('P', 'Q', 'C')
ylab <- c('Billion Persons', 'Trillions 2011US$', expression('Gt'~CO[2]~'/yr'))
titles <- c('Population', 'Gross World Product', 'Emissions')

p_ci <- list()
cols <- c('ci_base'='red', 'median_base'='red', 'point'='black')
th <- theme_bw(base_size=10) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
theme_set(th)

for (i in 1:length(dat)) {
  pred <- do.call(rbind, lapply(sim_out, function(ll) ll$out[ll$out$year %in% yrs, mod_out_names[i]]))
  q <- apply(pred, 2, quantile, probs=c(0.05, 0.5, 0.95)) # compute 95% CI and medians for each year
  q <- data.frame(cbind(yrs, t(q))) # append year to quantile df
  colnames(q) <- c('year', 'lb', 'median', 'ub')

# plot CIs and data
  p_ci[[i]] <- ggplot() + 
    geom_ribbon(data=q, aes(x=year, ymin=lb, ymax=ub, fill='red'), alpha=0.2) + 
    geom_point(data=dat[[i]], aes(x=year, y=value, color='black'), size=0.75) + 
    geom_line(data=q, aes(x=year, y=median, color='red')) + 
    scale_color_manual('', labels=c('Observation', 'Hindcast Median'), values=c('black', 'red')) +  
    scale_fill_manual('', labels=c('Hindcast 90% Interval'), values='red') + 
    scale_x_continuous('Year', expand=c(0, 0), limits=c(1820, 2019)) + 
    scale_y_continuous(ylab[i], expand=c(0, 0)) + 
    theme(legend.spacing.y = unit(-0.35, "cm"), plot.margin=margin(0.2, 0.2, 0, 0.1, unit='in'), legend.box='vertical') + 
      labs(tag=letters[i]) +
      guides(color=guide_legend(override.aes = list(linetype=c(0, 1), shape=c(16, NA))))
}

# extract legend
g <- ggplot_gtable(ggplot_build(p_ci[[1]]))
legidx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
legend <- g$grobs[[legidx]]

p_ci <- lapply(p_ci, function(p) p + theme(legend.position='none'))

g <- grid.arrange(p_ci[[1]], p_ci[[2]], p_ci[[3]], legend, ncol=2)

pdf('figures/hindcast_calibrated.pdf', width=6, height=6)
grid.draw(g)
dev.off()


png('figures/hindcast_calibrated.png', width=6, height=6, units='in', res=600)
grid.draw(g)
dev.off()
