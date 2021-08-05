library(IAMUQ)
library(readxl)
library(reshape2)
library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)
library(dplyr)

ssp_dat <- read_excel('data/cmip6_co2.xlsx')
ssp_dat <- as.data.frame(ssp_dat)[1:7,]
ssp_dat[, 'Scenario'] <- gsub('\\s.*', '', ssp_dat[, 'Scenario'])
ssp_dat[, 'Scenario'] <- gsub('([a-z]*)(\\d)$', '\\1\\.\\2', ssp_dat[,'Scenario'])
ssp_dat <- ssp_dat[, c(2, 6:(ncol(ssp_dat)-2))]

ssp_melt <- melt(ssp_dat, id.vars=c('Scenario'))
colnames(ssp_melt)[2] <- 'Year'
ssp_melt[, 'Year'] <- as.numeric(levels(ssp_melt[, 'Year']))[ssp_melt[, 'Year']]
ssp_melt[, 'value'] <- ssp_melt[, 'value'] /1000

yrs <- 2020:2100

dat <- lapply(iamdata, function(l) {l[l$year %in% 1820:2019,]})
obs <- dat[['emissions']][dat[['emissions']]$year %in% 2000:2019, ]

tol9qualitative=c("#88CCEE", "#44AA99", "#117733", "#332288", "#AA4499",  "#CC6677", "#882255", "#6699CC", "#999933")

cbbpsqualitative <- c("#000000", "#e79f00", "#9ad0f3", "#CC79A7", "#0072B2", "#009E73", "#F0E442", "#D55E00")

scenarios <- c('base', 'low', 'high', 'del_zc')
scen_labels <- c('Standard', 'Low Fossil Fuel', 'High Fossil Fuel', 'Delayed Zero-Carbon')

appendices <- c('', '-gwp', '-co2', '-pop', '-gwp-co2-pop')
app_names <- c('None', 'Economic~Growth', 'CO[2]~Emissions', 'Population', 'All')
app_titles <- parse(text=app_names)

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

p <- vector('list', length(appendices))
tmarg <- c(16.5, 16.5, 17.5, 17.5, 16.5) # set top margins for the marginal plot; this changes based on the types of characters in the plot titles

for (i in 1:length(appendices)) {
  sim_out <- vector('list', length(scenarios))
  names(sim_out) <- scenarios
  
  for (scen in scenarios) {
    # read simulation output file
    sim_out[[scen]] <- readRDS(paste0('output/sim_', scen, appendices[i], '.rds'))
  }
  
  emis <- lapply(sim_out, get_emissions, yrs=yrs)
  
  co2_q <- lapply(emis, emis_q, yrs=yrs)
  co2_q <- bind_rows(co2_q, .id='scenario')
  co2_q$case <- factor(co2_q$scenario, levels=scenarios, labels=scen_labels)
  
  co2_2100 <- lapply(emis, emis_dist, yr=2100)
  co2_2100 <- bind_rows(co2_2100, .id='scenario')
  co2_2100$case <- factor(co2_2100$scenario, levels=scenarios, labels=scen_labels)
  
  p_series <- ggplot() + geom_ribbon(data=co2_q, 
                                     aes(x=Year, ymin=lower, ymax=upper, fill=case), color=NA, alpha=0.3) + 
    geom_line(data=ssp_melt, aes(x=Year, y=value, color=Scenario)) +  
    geom_point(data=obs, aes(x=year, y=value), color='black', size=1) + 
    scale_y_continuous(expression(CO[2]~Emissions~(Gt~CO[2]/yr)), limits=c(0, 120), expand=c(0, 0), sec.axis=sec_axis(~., labels=NULL)) + 
    scale_color_manual('Emissions Scenario', values=tol9qualitative) + 
    scale_linetype_discrete('Marker') + 
    scale_fill_manual('Model Scenario', values=cbbpsqualitative) + 
    scale_x_continuous('Year', limits=c(2000, 2100), breaks=seq(2000, 2100, by=20), expand=c(0, 0)) + 
    theme_classic(base_size=10) +
    theme(plot.margin=unit(c(7, -0.35, 5, 2), 'mm'), legend.position='bottom', 
          legend.box='vertical', legend.box.just = 'left', legend.spacing = unit(0, 'cm'),
          legend.key.size=unit(0.5, 'cm'), legend.margin=margin(c(0.5, 0, 0, 0))) + 
    guides(color = guide_legend(order=1, ncol=2, byrow=FALSE, override_aes=list(size=5)), 
           linetype = guide_legend(order=2), fill = guide_legend(order=3, nrow=2, byrow=FALSE)) +
    labs(tag=letters[i]) + ggtitle(app_titles[i])
  
  p_marg <- ggplot() + stat_density(data=co2_2100, 
                                    aes(x=value, fill=case, color=case), geom='line', 
                                    position='identity') + 
    scale_x_continuous(limits=c(0, 120), expand=c(0, 0)) + 
    scale_y_continuous(expand=c(0, 0)) + 
    coord_flip() + 
    scale_fill_manual('', values=cbbpsqualitative) + 
    scale_color_manual('', values=cbbpsqualitative) + 
    theme_classic(base_size=10) +
    theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank(),  
          panel.border=element_blank(), plot.margin=unit(c(tmarg[i], 1, 12.5, -0.3), 'mm'), 
          axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), 
          axis.line.x=element_blank(), legend.position='none')
  
  p[[i]] <- gtable_row('scenarios', 
                  grobs=list(ggplotGrob(p_series + theme(legend.position='none')), 
                             ggplotGrob(p_marg)), height=unit(3, 'in'), widths=unit(c(2.7, 0.8), 'in'), 
                  z=c(2,1))

}
  
# extract legend from time series plot
g <- ggplot_gtable(ggplot_build(p_series))
legidx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
legend <- g$grobs[[legidx]]

fig <- arrangeGrob(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], legend, nrow=3, ncol=2)

pdf('figures/emissions-all-expert.pdf', height=8, width=7)
grid.draw(fig)
dev.off()

png('figures/emissions-all-expert.png', height=8, width=7, units='in', res=300)
grid.draw(fig)
dev.off()
