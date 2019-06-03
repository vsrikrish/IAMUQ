library(BAUcalib)
library(readr)
library(reshape2)
library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)
library(dplyr)

ssp_dat <- read_csv('data/ssp_db.csv', col_types='cccccnnnnnnnnnnnl')
ssp_dat <- as.data.frame(ssp_dat)
ssp_dat[, 'Scenario'] <- gsub('-Baseline', '', ssp_dat[, 'Scenario'])
ssp_dat <- ssp_dat[ssp_dat[, ncol(ssp_dat)], c(2, 6:(ncol(ssp_dat)-1))]
ssp_dat[, '2000'] <- NA
ssp_dat <- ssp_dat[, c(1, ncol(ssp_dat), 2:(ncol(ssp_dat)-1))]

rcp_dat <- read_csv('data/rcp_db.csv', col_types='ccccnnnnnnnnnnnnc')
rcp_dat <- as.data.frame(rcp_dat)
rcp_dat <- rcp_dat[, c(2, 6:(ncol(rcp_dat))-1)]

scen_dat <- rbind(ssp_dat, rcp_dat)
scen_dat[, 'Run'] <- 1:nrow(scen_dat)

scen_melt <- melt(scen_dat, id.vars=c('Scenario', 'Run'))
colnames(scen_melt)[3] <- 'Year'
scen_melt[, 'Year'] <- as.numeric(levels(scen_melt[, 'Year']))[scen_melt[, 'Year']]
scen_melt[grepl('SSP', scen_melt[, 'Scenario']), 'value'] <- scen_melt[grepl('SSP', scen_melt[, 'Scenario']), 'value'] / (3.67 * 1000)


yrs <- 1700:2100
obs <- baudata[['emissions']][baudata[['emissions']]$year %in% 2000:2014, ]

tol9qualitative=c("#88CCEE", "#332288", "#117733", "#44AA99", "#AA4499",  "#CC6677", "#882255", "#6699CC", "#999933")

cbbpsqualitative <- c("#000000", "#e79f00", "#9ad0f3", "#CC79A7", "#0072B2", "#009E73", "#F0E442", "#D55E00")

cases <- c('base', 'low', 'high', 'alt_zc')
case_labels <- c('Base', 'Low Fossil Fuel', 'High Fossil Fuel', 'Delayed Zero-Carbon')

get_emissions <- function(case) {
  model_pred <- readRDS(paste0('output/reject_sim_', case, '.rds'))

  do.call(cbind, lapply(model_pred, function(l) l$C))
}

emis_q <- function(emis, yrs) {
  co2_q <- data.frame(Year=yrs, t(apply(emis, 1, quantile, probs=c(0.05, 0.95))))
  colnames(co2_q)[2:3] <- c('lower', 'upper')

  co2_q
}

emis_dist <- function(emis, yr) {
    data.frame(value=emis[which(yrs == yr), ])
}

emis <- lapply(cases, get_emissions)
names(emis) <- cases

co2_q <- lapply(emis, emis_q, yrs=yrs)
co2_q <- bind_rows(co2_q, .id='case')
co2_q$case <- factor(co2_q$case, levels=cases, labels=case_labels)

co2_2100 <- lapply(emis, emis_dist, yr=2100)
co2_2100 <- bind_rows(co2_2100, .id='case')
co2_2100$case <- factor(co2_2100$case, levels=cases, labels=case_labels)

th <- theme_bw(base_size=10) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_text(size=9))
theme_set(th)

p_series <- ggplot() + geom_ribbon(data=co2_q, aes(x=Year, ymin=lower, ymax=upper, fill=case), color=NA, alpha=0.3) + geom_line(data=scen_melt, aes(x=Year, y=value, group=factor(Run), color=Scenario)) +  geom_point(data=obs, aes(x=year, y=emissions), color='black', size=1) + scale_y_continuous(expression(CO[2]~Emissions~(Gt~C/yr)), limits=c(0, 40), expand=c(0.001, 0.001)) + scale_color_manual('Emissions Scenario', values=tol9qualitative) + scale_linetype_discrete('Marker') + scale_fill_manual('Model Scenario', values=cbbpsqualitative) + scale_x_continuous('Year', limits=c(2000, 2100), breaks=seq(2000, 2100, by=20), expand=c(0.001, 0.001)) + theme(plot.margin=unit(c(0.6, 0, 0.5, 0.08), 'in'), legend.position='bottom', legend.box='vertical', legend.box.just = 'left', legend.spacing = unit(-0.2, 'cm')) + guides(color = guide_legend(order=1, nrow=3, byrow=FALSE, override_aes=list(size=5)), linetype = guide_legend(order=2), fill = guide_legend(order=3, nrow=2, byrow=FALSE))

p_marg <- ggplot() + stat_density(data=co2_2100, aes(x=value, fill=case, color=case), alpha=0.3, geom='area', position='identity') + scale_x_continuous(limits=c(0, 40), expand=c(0.001, 0.001)) + scale_y_continuous(expand=c(0.001, 0.001)) + coord_flip() + scale_fill_manual('', values=cbbpsqualitative) + scale_color_manual('', values=cbbpsqualitative) + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), panel.border=element_blank(), plot.margin=unit(c(0.605, 0.2, 0.815, -0.02), 'in'), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.line.y=element_line(color='black'), legend.position='none')

p <- gtable_row('scenarios', grobs=list(ggplotGrob(p_series + theme(legend.position='none')), ggplotGrob(p_marg)), height=unit(3, 'in'), widths=unit(c(2.7, 0.8), 'in'), z=c(2,1))


# extract legend from time series plot
g <- ggplotGrob(p_series)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$width)

fig <- arrangeGrob(p, legend, ncol=1,
          heights=unit.c(unit(1, 'npc') - lheight, lheight)
         )

pdf('figures/Fig2-emissions.pdf', height=3.5, width=3.5)
grid.draw(fig)
dev.off()

png('figures/Fig2-emissions.png', height=3.5, width=3.5, units='in', res=600)
grid.draw(fig)
dev.off()

th <- theme_bw(base_size=18) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
theme_set(th)

p_series <- ggplot() + geom_ribbon(data=co2_q, aes(x=Year, ymin=lower, ymax=upper, fill=case), color=NA, alpha=0.3) + geom_line(data=scen_melt, aes(x=Year, y=value, group=factor(Run), color=Scenario)) +  geom_point(data=obs, aes(x=year, y=emissions), color='black', size=1.5) + scale_y_continuous(expression(CO[2]~Emissions~(Gt~C/yr)), limits=c(0, 40), expand=c(0.001, 0.001)) + scale_color_manual('Emissions Scenario', values=tol9qualitative) + scale_linetype_discrete('Marker') + scale_fill_manual('Model Scenario', values=cbbpsqualitative) + scale_x_continuous(limits=c(2000, 2100), breaks=seq(2000, 2100, by=20), expand=c(0.001, 0.001)) + theme(plot.margin=unit(c(0.5, 0, 0.2, 0.1), 'in'), legend.position='bottom', legend.box='vertical', legend.box.just = 'left', legend.spacing = unit(-0.2, 'cm')) + guides(color = guide_legend(order=1, nrow=3, byrow=FALSE, override_aes=list(size=5)), linetype = guide_legend(order=2), fill = guide_legend(order=3, nrow=2, byrow=FALSE))

p_marg <- ggplot() + stat_density(data=co2_2100, aes(x=value, fill=case, color=case), alpha=0.3, geom='area', position='identity') + scale_x_continuous(limits=c(0, 40), expand=c(0.001, 0.001)) + scale_y_continuous(expand=c(0.001, 0.001)) + coord_flip() + scale_fill_manual('', values=cbbpsqualitative) + scale_color_manual('', values=cbbpsqualitative) + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), panel.border=element_blank(), plot.margin=unit(c(0.5, 0.2, 0.7, -0.05), 'in'), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.line.y=element_line(color='black'), legend.position='none')

p <- gtable_row('scenarios', grobs=list(ggplotGrob(p_series + theme(legend.position='none')), ggplotGrob(p_marg)), height=unit(5.5, 'in'), widths=unit(c(5, 1.5), 'in'), z=c(2,1))


# extract legend from time series plot
g <- ggplotGrob(p_series)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$width)

fig <- arrangeGrob(p, legend, ncol=1,
          heights=unit.c(unit(1, 'npc') - lheight, lheight)
         )

png('figures/ssp-slide.png', height=7, width=6.5, units='in', res=600)
grid.draw(fig)
dev.off()
