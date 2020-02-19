library(ggplot2)
library(reshape2)
library(IAMUQ)
library(grid)
library(gridExtra)

cases <- c('base', 'short')

nsamp <- 1e5

yrs <- 2015:2100
dat <- lapply(iamdata, function(l) l[l$year %in% 1820:2014, ])

sim_out <- vector('list', length(cases))
names(sim_out) <- cases

for (case in cases) {
  # get posterior samples
  sim_out[[case]] <- readRDS(paste0('output/sim_', case, '-gwp-co2.rds'))
  sim_out[[case]] <- lapply(sim_out[[case]], function(l) l$out)
}


mod_out_names <- c('P', 'Q', 'C')
ylab <- c('Billions', 'Trillions 2011US$', 'Gt C/yr')
titles <- c('Population', 'Gross World Product', 'Emissions')

p_ci <- list()
cols <- c('ci_base'='red', 'ci_short'='blue', 'median_base'='red', 'median_short'='blue', 'point'='black')
th <- theme_bw(base_size=10) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.title=element_text(size=10), legend.text=element_text(size=9))
theme_set(th)

for (i in 1:length(dat)) {
  pred <- lapply(sim_out, function(ll) do.call(rbind, lapply(ll, function(l) l[l$year %in% yrs, mod_out_names[i]])))
  q <- lapply(pred, function(l) apply(l, 2, quantile, probs=c(0.05, 0.5, 0.95))) # compute 95% CI and medians for each year
  q <- lapply(q, function(l) as.data.frame(cbind(yrs, t(l)))) # append year to quantile df
  q <- lapply(q, setNames, nm=c('year', 'lb', 'median', 'ub'))

# plot CIs and data
  p_ci[[i]] <- ggplot() + geom_ribbon(data=q[['base']], aes(x=year, ymin=lb, ymax=ub, fill='ci_base'), alpha=0.2) + geom_ribbon(data=q[['short']], aes(x=year, ymin=lb, ymax=ub, fill='ci_short'), alpha=0.2) + geom_point(data=dat[[i]], aes(x=year, y=value, color='point'), size=0.75) + geom_line(data=q[['base']], aes(x=year, y=median, color='median_base')) + geom_line(data=q[['short']], aes(x=year, y=median, color='median_short')) + scale_color_manual('', labels=c('Median (1820-2014)', 'Median (1950-2014)', 'Observation'), values=cols) +  scale_fill_manual('', labels=c('90% CI (1820-2014)', '90% CI (1950-2014)'), values=cols) + scale_x_continuous('Year', expand=c(0.001, 0.001), limits=c(1950, 2100)) + scale_y_continuous(ylab[i], expand=c(0.001, 0.001)) + guides(color = guide_legend(override.aes = list(linetype=c(1, 1, NA), shape=c(NA, NA, 16)))) + theme(legend.spacing = unit(-0.15, "cm"), plot.margin=margin(0.2, 0.2, 0, 0.1, unit='in'), legend.box='vertical') + ggtitle(paste0(letters[i], ') ', titles[i]))
}

grid_arrange_shared_legend <- function(..., nrow = length(list(...)), ncol = 1, position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)

}


g <- grid_arrange_shared_legend(p_ci[[1]], p_ci[[2]], p_ci[[3]], ncol=3, nrow=1, position='bottom')

pdf('figures/hindcast_length.pdf', width=6, height=4)
grid.draw(g)
dev.off()


png('figures/hindcast_length.png', width=6, height=4, units='in', res=600)
grid.draw(g)
dev.off()

p_ci <- lapply(p_ci, function(p) p + theme_bw(base_size=16))

g <- grid_arrange_shared_legend(p_ci[[1]], p_ci[[2]], p_ci[[3]], ncol=3, nrow=1, position='bottom')

png('figures/hindcast_slide.png', width=8, height=6, units='in', res=600)
grid.draw(g)
dev.off()
