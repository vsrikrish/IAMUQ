library(ggplot2)
library(reshape2)
library(IAMUQ)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(dplyr)

cases <- c('base', 'short', 'vshort')

nsamp <- 1e5

yrs <- 2020:2100
dat <- lapply(iamdata, function(l) l[l$year %in% 1820:2019, ])

sim_out <- vector('list', length(cases))
names(sim_out) <- cases

for (case in cases) {
  # get posterior samples
  sim_out[[case]] <- readRDS(paste0('output/sim_', case, '-gwp-co2-pop.rds'))
  sim_out[[case]] <- lapply(sim_out[[case]], function(l) l$out)
}


mod_out_names <- c('P', 'Q', 'C')
ylab <- c('Billions', 'Trillions 2011US$', expression(Gt~CO[2]~'/yr'))
titles <- c('Population', 'Gross World Product', 'Emissions')

p_ci <- list()
cols <- brewer.pal(3, 'Dark2')
th <- theme_bw(base_size=10) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
theme_set(th)

for (i in 1:length(dat)) {
  pred <- lapply(sim_out, function(ll) do.call(rbind, lapply(ll, function(l) l[l$year %in% yrs, mod_out_names[i]])))
  q <- lapply(pred, function(l) apply(l, 2, quantile, probs=c(0.05, 0.5, 0.95))) # compute 95% CI and medians for each year
  q <- lapply(q, function(l) as.data.frame(cbind(yrs, t(l)))) # append year to quantile df
  q <- lapply(q, setNames, nm=c('year', 'lb', 'median', 'ub'))
  q_all <- bind_rows(q, .id="case")
  q_all$case <- ordered(q_all$case, levels=c('base', 'short', 'vshort'))

# plot CIs and data
  p_ci[[i]] <- ggplot() + geom_ribbon(data=q_all, aes(x=year, ymin=lb, ymax=ub, fill=case), alpha=0.3) + 
    geom_point(data=dat[[i]], aes(x=year, y=value), color='black', size=0.75) + 
    geom_line(data=q_all, aes(x=year, y=median, color=case)) + 
    scale_color_manual('Calibration Period', labels=c('1820-2019', '1950-2019', '2000-2019'), values=cols) +  
    scale_fill_manual('Calibration Period', labels=c('1820-2019', '1950-2019', '2000-2019'), values=cols) + 
    scale_x_continuous('Year', expand=c(0.001, 0.001), limits=c(2000, 2100)) + 
    scale_y_continuous(ylab[i], expand=c(0.001, 0.001), sec.axis=sec_axis(~., labels=NULL)) + 
    theme(legend.spacing = unit(-0.15, "cm"), plot.margin=margin(0.2, 0.2, 0, 0.1, unit='in')) + 
      labs(tag=letters[i])
}


# extract legend
g <- ggplot_gtable(ggplot_build(p_ci[[1]]))
legidx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
legend <- g$grobs[[legidx]]

p_ci <- lapply(p_ci, function(p) p + theme(legend.position='none'))

g <- grid.arrange(p_ci[[1]], p_ci[[2]], p_ci[[3]], legend, ncol=2)

pdf('figures/hindcast_length.pdf', width=6, height=6)
grid.draw(g)
dev.off()


png('figures/hindcast_length.png', width=6, height=6, units='in', res=600)
grid.draw(g)
dev.off()
