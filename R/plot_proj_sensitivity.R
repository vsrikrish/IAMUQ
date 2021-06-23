library(IAMUQ)
library(reshape2)
library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)
library(dplyr)

yrs <- 2015:2100

dat <- lapply(iamdata, function(l) {l[l$year %in% 1820:2014,]})
obs <- dat[['emissions']][dat[['emissions']]$year %in% 2000:2014, ]

scenarios <- c('base', 'alt_multi')
scen_labels <- c('Default Priors', 'Alternate Priors')

cbbpsqualitative <- c("#000000", "#e79f00", "#9ad0f3", "#CC79A7", "#0072B2", "#009E73", "#F0E442", "#D55E00")

get_emissions <- function(sim_out, yrs) {
  do.call(cbind, lapply(sim_out, function(l) l$out$C[l$out$year %in% yrs]))
}

emis_q <- function(emis, yrs) {
  co2_q <- data.frame(Year=yrs, t(apply(emis, 1, quantile, probs=c(0.05, 0.95))))
  colnames(co2_q)[2:3] <- c('lower', 'upper')

  co2_q
}

emis_dist <- function(emis, yr) {
    data.frame(value=emis[which(yrs == yr), ])
}

sim_out <- vector('list', length(scenarios))
names(sim_out) <- scenarios

for (scen in scenarios) {
  # read simulation output file
  sim_out[[scen]] <- readRDS(paste0('output/sim_', scen, '-gwp-co2.rds'))
}

emis <- lapply(sim_out, get_emissions, yrs=yrs)

co2_q <- lapply(emis, emis_q, yrs=yrs)
co2_q <- bind_rows(co2_q, .id='scenario')
co2_q$case <- factor(co2_q$scenario, levels=scenarios, labels=scen_labels, ordered=TRUE)

co2_2100 <- lapply(emis, emis_dist, yr=2100)
co2_2100 <- bind_rows(co2_2100, .id='scenario')
co2_2100$case <- factor(co2_2100$scenario, levels=scenarios, labels=scen_labels, ordered=TRUE)

p_series <- ggplot() + geom_ribbon(data=co2_q, 
                                   aes(x=Year, ymin=lower, ymax=upper, fill=case), color=NA, alpha=0.3) + 
  geom_point(data=obs, aes(x=year, y=value), color='black', size=1) + 
  scale_y_continuous(expression(CO[2]~Emissions~(Gt~C/yr)), limits=c(0, 40), expand=c(0, 0)) + 
  scale_color_manual('Emissions Scenario', values=tol9qualitative) + 
  scale_linetype_discrete('Marker') + 
  scale_fill_manual('Model Scenario', values=cbbpsqualitative) + 
  scale_x_continuous('Year', limits=c(2000, 2100), breaks=seq(2000, 2100, by=20), expand=c(0, 0)) + 
  theme_classic(base_size=10) +
  theme(plot.margin=unit(c(7, -0.35, 5, 2), 'mm'), legend.position='bottom', 
        legend.box='vertical', legend.box.just = 'left', legend.spacing = unit(0, 'cm'),
        legend.key.size=unit(0.5, 'cm'), legend.margin=margin(c(0.5, 0, 0, 0))) + 
  guides(color = guide_legend(order=1, ncol=2, byrow=FALSE, override_aes=list(size=5)), 
         linetype = guide_legend(order=2), fill = guide_legend(order=3, nrow=2, byrow=FALSE))

p_marg <- ggplot() + stat_density(data=co2_2100, 
                                  aes(x=value, fill=case, color=case), geom='line', 
                                  position='identity') + 
  scale_x_continuous(limits=c(0, 40), expand=c(0, 0)) + 
  scale_y_continuous(expand=c(0, 0)) + 
  coord_flip() + 
  scale_fill_manual('', values=cbbpsqualitative) + 
  scale_color_manual('', values=cbbpsqualitative) + 
  theme_classic(base_size=5) +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), 
        panel.border=element_blank(), plot.margin=unit(c(7, 1, 13.2, -0.05), 'mm'), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), 
        axis.line.x=element_blank(), legend.position='none')

p <- gtable_row('scenarios', 
                grobs=list(ggplotGrob(p_series + theme(legend.position='none')), 
                           ggplotGrob(p_marg)), height=unit(3, 'in'), widths=unit(c(2.7, 0.8), 'in'), 
                z=c(2,1))

# extract legend from time series plot
g <- ggplotGrob(p_series)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$width)

fig <- arrangeGrob(p, legend, ncol=1,
          heights=unit.c(unit(1, 'npc') - lheight, lheight)
         )

pdf('figures/emissions_proj_sensitivity.pdf', height=3.5, width=3.5)
grid.draw(fig)
dev.off()

png('figures/emissions_proj_sensitivity.png', height=3.5, width=3.5, units='in', res=600)
grid.draw(fig)
dev.off()
