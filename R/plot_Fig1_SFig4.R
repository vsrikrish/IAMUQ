library(ggplot2)
library(reshape2)
library(BAUcalib)
library(grid)
library(gridExtra)

model_pred <- list(base=readRDS('output/reject_sim_base.rds'), short=readRDS('output/reject_sim_short.rds'))

yrs <- 1900:2100
dat <- lapply(baudata, function(l) l[l$year %in% yrs, ])
dat <- lapply(dat, setNames, nm=c('year', 'value'))
mod_out_names <- c('P', 'Q', 'C')
ylab <- c('Billions', 'Trillions 2011US$', 'Gt C/yr')
titles <- c('Population', 'Gross World Product', 'Emissions')

p_ci <- list()
cols <- c('ci_base'='red', 'ci_short'='blue', 'median_base'='red', 'median_short'='blue', 'point'='black')
th <- theme_bw(base_size=10) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
theme_set(th)

for (i in 1:length(dat)) {
  pred <- lapply(model_pred, function(ll) do.call(rbind, lapply(ll, function(l) l[l$year %in% yrs, mod_out_names[i]])))
  q <- lapply(pred, function(l) apply(l, 2, quantile, probs=c(0.05, 0.5, 0.95))) # compute 95% CI and medians for each year
  q <- lapply(q, function(l) as.data.frame(cbind(yrs, t(l)))) # append year to quantile df
  q <- lapply(q, setNames, nm=c('year', 'lb', 'median', 'ub'))

# plot CIs and data
  p_ci[[i]] <- ggplot() + geom_ribbon(data=q[['base']], aes(x=year, ymin=lb, ymax=ub, fill='ci_base'), alpha=0.2) + geom_ribbon(data=q[['short']], aes(x=year, ymin=lb, ymax=ub, fill='ci_short'), alpha=0.2) + geom_point(data=dat[[i]], aes(x=year, y=value, color='point'), size=0.75) + geom_line(data=q[['base']], aes(x=year, y=median, color='median_base')) + geom_line(data=q[['short']], aes(x=year, y=median, color='median_short')) + scale_color_manual('', labels=c('Median (1820-2014)', 'Median (1950-2014)', 'Observation'), values=cols) +  scale_fill_manual('', labels=c('90% CI (1820-2014)', '90% CI (1950-2014)'), values=cols) + scale_x_continuous('Year', expand=c(0.001, 0.001)) + scale_y_continuous(ylab[i], expand=c(0.001, 0.001)) + guides(color = guide_legend(override.aes = list(linetype=c(1, 1, NA), shape=c(NA, NA, 16)))) + theme(legend.spacing = unit(-0.15, "cm"), plot.margin=margin(0.2, 0.2, 0, 0.1, unit='in')) + ggtitle(paste0(letters[i], ') ', titles[i]))
}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

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

pdf('figures/FigS4-hindcast.pdf', width=8, height=3)
grid.draw(g)
dev.off()

png('figures/FigS4-hindcast.png', width=8, height=3, units='in', res=600)
grid.draw(g)
dev.off()

for (i in 1:length(dat)) {
  pred <- do.call(rbind, lapply(model_pred[['base']], function(l) l[l$year %in% yrs, mod_out_names[i]]))
  q <- apply(pred, 2, quantile, probs=c(0.05, 0.5, 0.95)) # compute 95% CI and medians for each year
  q <- as.data.frame(cbind(yrs, t(q))) # append year to quantile df
  colnames(q) <- c('year', 'lb', 'median', 'ub')

# plot CIs and data
  p_ci[[i]] <- ggplot() + geom_ribbon(data=q, aes(x=year, ymin=lb, ymax=ub, fill='ci_base'), alpha=0.2) + geom_point(data=dat[[i]], aes(x=year, y=value, color='point'), size=0.75) + geom_line(data=q, aes(x=year, y=median, color='median_base')) + scale_color_manual('', labels=c('Median', 'Observation'), values=cols) +  scale_fill_manual('', labels=c('90% CI'), values=cols) + scale_x_continuous('Year', expand=c(0.001, 0.001)) + scale_y_continuous(ylab[i], expand=c(0.001, 0.001)) + guides(color = guide_legend(override.aes = list(linetype=c(1, NA), shape=c(NA, 16)))) + theme(legend.spacing = unit(-0.15, "cm"), plot.margin=margin(0.2, 0.2, 0, 0.1, unit='in')) + ggtitle(paste0(letters[i], ') ', titles[i]))


  # compute hit rate and print
  print(sum(vapply(1:nrow(dat[[i]]), function(j) (q[q$year == dat[[i]]$year[j], 'lb'] <= dat[[i]]$value[j]) && (q[q$year == dat[[i]]$year[j], 'ub'] >= dat[[i]]$value[j]), logical(1)), na.rm=T)/nrow(dat[[i]]) * 100)
}

g <- grid_arrange_shared_legend(p_ci[[1]], p_ci[[2]], p_ci[[3]], ncol=3, nrow=1, position='bottom')

pdf('figures/Fig1-hindcast.pdf', width=7, height=3.5)
grid.draw(g)
dev.off()

png('figures/Fig1-hindcast.png', width=7, height=3.5, units='in', res=600)
grid.draw(g)
dev.off()



th <- theme_bw(base_size=15) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.title=element_text(size=15))
theme_set(th)

for (i in 1:length(dat)) {
  pred <- do.call(rbind, lapply(model_pred[['base']], function(l) l[l$year %in% yrs, mod_out_names[i]]))
  q <- apply(pred, 2, quantile, probs=c(0.05, 0.5, 0.95)) # compute 95% CI and medians for each year
  q <- as.data.frame(cbind(yrs, t(q))) # append year to quantile df
  colnames(q) <- c('year', 'lb', 'median', 'ub')

  # plot CIs and data
  p_ci[[i]] <- ggplot() + geom_ribbon(data=q, aes(x=year, ymin=lb, ymax=ub, fill='ci_base'), alpha=0.2) + geom_point(data=dat[[i]], aes(x=year, y=value, color='point'), size=0.75) + geom_line(data=q, aes(x=year, y=median, color='median_base')) + scale_color_manual('', labels=c('Median', 'Observation'), values=cols) +  scale_fill_manual('', labels=c('90% CI'), values=cols) + scale_x_continuous('Year', expand=c(0.001, 0.001)) + scale_y_continuous(ylab[i], expand=c(0.001, 0.001)) + guides(color = guide_legend(override.aes = list(linetype=c(1, NA), shape=c(NA, 16)))) + theme(legend.spacing = unit(-0.15, "cm"), plot.margin=margin(0.1, 0.2, 0, 0.1, unit='in')) + ggtitle(paste0(letters[i], ') ', titles[i]))

}


g <- grid_arrange_shared_legend(p_ci[[1]], p_ci[[2]], p_ci[[3]], ncol=3, nrow=1, position='bottom')

png('figures/hindcast-slide.png', width=9, height=5, units='in', res=600)
grid.draw(g)
dev.off()
