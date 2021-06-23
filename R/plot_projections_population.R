library(IAMUQ)
library(readxl)
library(reshape2)
library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)
library(dplyr)

ssp_dat <- read_excel('data/ssp_pop.xlsx')
ssp_dat <- as.data.frame(ssp_dat)[1:7,]
ssp_dat[, 'Scenario'] <- gsub('\\s.*', '', ssp_dat[, 'Scenario'])
ssp_dat[, 'Scenario'] <- gsub('([a-z]*)(\\d)$', '\\1\\.\\2', ssp_dat[,'Scenario'])
ssp_dat <- ssp_dat[, c(2, 6:(ncol(ssp_dat)-2))]

ssp_melt <- melt(ssp_dat, id.vars=c('Scenario'))
colnames(ssp_melt)[2] <- 'Year'
ssp_melt[, 'Year'] <- as.numeric(levels(ssp_melt[, 'Year']))[ssp_melt[, 'Year']]

yrs <- 2015:2100

dat <- lapply(iamdata, function(l) {l[l$year %in% 1820:2014,]})
obs <- dat[['pop']][dat[['pop']]$year %in% 2000:2014, ]

tol9qualitative=c("#88CCEE", "#44AA99", "#117733", "#332288", "#AA4499",  "#CC6677", "#882255", "#6699CC", "#999933")

cbbpsqualitative <- c("#000000", "#e79f00", "#9ad0f3", "#CC79A7", "#0072B2", "#009E73", "#F0E442", "#D55E00")

scenarios <- c('base', 'low', 'high', 'del_zc')
scen_labels <- c('Standard', 'Low Fossil Fuel', 'High Fossil Fuel', 'Delayed Zero-Carbon')

get_pop <- function(sim_out, yrs) {
  do.call(cbind, lapply(sim_out, function(l) l$out$P[l$out$year %in% yrs]))
}

out_q <- function(out, yrs) {
  out_q <- data.frame(Year=yrs, t(apply(out, 1, quantile, probs=c(0.05, 0.95))))
  colnames(out_q)[2:3] <- c('lower', 'upper')
  
  out_q
}

out_dist <- function(out, yr) {
  data.frame(value=out[which(yrs == yr), ])
}

sim_out <- vector('list', length(scenarios))
names(sim_out) <- scenarios

for (scen in scenarios) {
  # read simulation output file
  sim_out[[scen]] <- readRDS(paste0('output/sim_', scen, '-gwp-co2.rds'))
}

pop <- lapply(sim_out, get_pop, yrs=yrs)

pop_q <- lapply(pop, out_q, yrs=yrs)
pop_q <- bind_rows(pop_q, .id='scenario')
pop_q$case <- factor(pop_q$scenario, levels=scenarios, labels=scen_labels)

pop_2100 <- lapply(pop, out_dist, yr=2100)
pop_2100 <- bind_rows(pop_2100, .id='scenario')
pop_2100$case <- factor(pop_2100$scenario, levels=scenarios, labels=scen_labels)

p_series <- ggplot() + geom_ribbon(data=pop_q, 
                                   aes(x=Year, ymin=lower, ymax=upper, fill=case), color=NA, alpha=0.3) + 
  geom_line(data=ssp_melt, aes(x=Year, y=value / 1000, color=Scenario)) +  
  geom_point(data=obs, aes(x=year, y=value), color='black', size=1) + 
  scale_y_continuous('Global Population (billion persons)', limits=c(6, 11.75), expand=c(0, 0)) + 
  scale_color_manual('SSP Scenario', values=tol9qualitative) + 
  scale_linetype_discrete('Marker') + 
  scale_fill_manual('Model Scenario', values=cbbpsqualitative) + 
  scale_x_continuous('Year', limits=c(2000, 2100), breaks=seq(2000, 2100, by=20), expand=c(0, 0)) + 
  theme_classic(base_size=10) +
  theme(plot.margin=unit(c(7, -0.35, 5, 2), 'mm'), legend.position='bottom', 
        legend.box='vertical', legend.box.just = 'left', legend.spacing = unit(0, 'cm'),
        legend.key.size=unit(0.5, 'cm'), legend.margin=margin(c(0.5, 0, 0, 0))) + 
  guides(color = guide_legend(order=1, ncol=2, byrow=FALSE, override_aes=list(size=5)), 
         linetype = guide_legend(order=2), fill = guide_legend(order=3, nrow=2, byrow=FALSE))

p_marg <- ggplot() + stat_density(data=pop_2100, 
                                  aes(x=value, fill=case, color=case), geom='line', 
                                  position='identity') + 
  scale_x_continuous(limits=c(6, 11.5), expand=c(0, 0)) + 
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

pdf('figures/population-SSP.pdf', height=4, width=3.5)
grid.draw(fig)
dev.off()

png('figures/population-SSP.png', height=4, width=3.5, units='in', res=600)
grid.draw(fig)
dev.off()
