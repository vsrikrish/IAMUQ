library(IAMUQ)
library(readxl)
library(reshape2)
library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)
library(dplyr)

yrs <- 2015:2100

ssp_dat <- read_excel('data/cmip6_co2.xlsx')
ssp_dat <- as.data.frame(ssp_dat)[1:7,]
ssp_dat[, 'Scenario'] <- gsub('\\s.*', '', ssp_dat[, 'Scenario'])
ssp_dat[, 'Scenario'] <- gsub('([a-z]*)(\\d)$', '\\1\\.\\2', ssp_dat[,'Scenario'])
ssp_dat <- ssp_dat[, c(2, 6:(ncol(ssp_dat)-2))]

ssp_cum <- apply(ssp_dat[,-1], 1, function(s) approx(x=as.numeric(colnames(ssp_dat)[-1]), y=s, xout=2018:2100))
ssp_cum <- lapply(ssp_cum, function(l) sum(l$y))
names(ssp_cum) <- ssp_dat[,'Scenario']
ssp_cum_melt <- melt(ssp_cum)
ssp_cum_melt$value <- ssp_cum_melt$value / (3.67 * 1000)
colnames(ssp_cum_melt)[2] <- 'Scenario'

ssp_melt <- melt(ssp_dat, id.vars=c('Scenario'))
colnames(ssp_melt)[2] <- 'Year'
ssp_melt[, 'Year'] <- as.numeric(levels(ssp_melt[, 'Year']))[ssp_melt[, 'Year']]
ssp_melt[, 'value'] <- ssp_melt[, 'value'] / (3.67 * 1000)
ssp_melt <- ssp_melt[ssp_melt$Year == 2100, ]


scenarios <- c('base', 'low', 'high', 'del_zc')
scen_labels <- c('Standard', 'Low Fossil Fuel', 'High Fossil Fuel', 'Delayed Zero-Carbon')

get_emissions <- function(sim_out) {
  do.call(cbind, lapply(sim_out, function(l) l$out$C))
}

emis_q <- function(emis, yrs) {
  co2_q <- data.frame(Year=yrs, t(apply(emis, 1, quantile, probs=c(0.05, 0.95))))
  colnames(co2_q)[2:3] <- c('lower', 'upper')

  co2_q
}

emis_dist <- function(emis, yrs, yr) {
    data.frame(value=emis[which(yrs == yr), ])
}

sim_out <- vector('list', length(scenarios))
names(sim_out) <- scenarios

for (scen in scenarios) {
  # read simulation output file
  sim_out[[scen]] <- readRDS(paste0('output/sim_', scen, '-gwp-co2.rds'))
}

emis <- lapply(sim_out, get_emissions)

co2_2100 <- lapply(emis, emis_dist, yrs=2015:2100, yr=2100)
co2_2100 <- bind_rows(co2_2100, .id='scenario')
co2_2100$case <- factor(co2_2100$scenario, levels=scenarios, labels=scen_labels)

co2_cum <- lapply(emis, cum_co2, yrs=2015:2100, start=2018)
co2_cum <- lapply(co2_cum, function(l) data.frame(value=l))
co2_cum <- bind_rows(co2_cum, .id='scenario')
co2_cum$case <- factor(co2_cum$scenario, levels=scenarios, labels=scen_labels)

cbbpsqualitative <- c("#000000", "#e79f00", "#9ad0f3", "#CC79A7", "#0072B2", "#009E73", "#F0E442", "#D55E00")

tol9qualitative=c("#88CCEE", "#332288", "#117733", "#44AA99", "#AA4499",  "#CC6677", "#882255", "#6699CC", "#999933")

#p <- ggplot() + geom_vline(data=ssp_melt, aes(xintercept=value), linetype='dotted') + stat_ecdf(data=co2_2100, aes(x=value, color=case))  + geom_text(data=ssp_melt, aes(x=value, y=1, label=Scenario), angle=90, size=2.5, fontface='bold', nudge_x=-1, nudge_y=ifelse(ssp_melt$value <= 20, -0.05, -0.9)) + scale_x_continuous(expression(CO[2]~Emissions~'in'~2100~(Gt~C/yr)), limits=c(0, 50)) + scale_y_continuous('Cumulative Probability') + scale_color_manual('Model Scenario', values=cbbpsqualitative) + theme(legend.position='bottom', legend.box='vertical', legend.box.just = 'left', legend.spacing = unit(-0.2, 'cm')) + guides(color = guide_legend(order=1, nrow=2, byrow=FALSE, override_aes=list(size=5)))


#pdf('figures/co2_2100-cdf.pdf', width=4.5, height=5)
#p
#dev.off()

#geom_text(data=ssp_cum_melt, aes(x=value, y=1, label=Scenario), angle=90, size=2.5, fontface='bold', nudge_x=c(-40, 40, 40, -40, -40, -40, -40), 
#          nudge_y=ifelse(ssp_cum_melt$value <= 1000, -0.05, -0.9)) + 

budget <- 1500 / 3.67

ssp_vjust <- c(-1, 1, 1, 1, -1, -1.25, 0)  
ssp_hjust <- c(0, -0.45, -0.5, 0.25, 0, -0.4, 0)

p <- ggplot() + 
  geom_rect(data=data.frame(xmin=0, xmax=budget, ymin=0, ymax=1), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='darkgreen', alpha=0.2) +
  geom_vline(data=ssp_cum_melt, aes(xintercept=value), linetype='solid', color='grey') + 
  stat_ecdf(data=co2_cum, aes(x=value, color=case))  + 
  scale_x_continuous(expression(CO[2]~Emissions~'from'~2018~'-'~2100~(Gt~C)), limits=c(0, 2500), expand=c(0, 0)) + 
  scale_y_continuous('Cumulative Probability', expand=c(0, 0),
                     labels=function(x) ifelse(x==as.integer(x), as.character(round(x)), as.character(round(x, 2)))) + 
  scale_color_manual('Model\nScenario', values=cbbpsqualitative) + 
  theme_classic(base_size=10) +
  theme(legend.position='bottom', legend.box='vertical', legend.box.just = 'left', 
        legend.spacing = unit(-0.2, 'cm'), plot.margin=unit(c(4, 1, 0, 1), 'lines'),
        legend.margin = margin(c(0, 0, 0, -15), 'lines')) +
#        legend.text = element_text(size=5), axis.text=element_text(size=5)) + 
  guides(color = guide_legend(order=1, nrow=2, byrow=FALSE))#, override_aes=list(size=5)))

for (i in 1:nrow(ssp_cum_melt)) {
  p <- p + annotation_custom(grob=grid::textGrob(label=ssp_cum_melt[i, 'Scenario'], hjust=0.5 + ssp_hjust[i], vjust=ssp_vjust[i], gp=gpar(fontsize=8)), 
                             xmin=ssp_cum_melt[i, 'value'], xmax=ssp_cum_melt[i, 'value'], ymin=1.1, ymax=1.1)
  if (i == 4) {
    p <- p + annotation_custom(grob=grid::linesGrob(arrow=grid::arrow(type='closed', ends='first', length=unit(1.5, 'mm')), gp=gpar(fill='black')),
                               xmax=ssp_cum_melt[i, 'value'], xmin=ssp_cum_melt[i, 'value'], 
                               ymin=1, ymax=1.09 - 0.03*ssp_vjust[i])
  } else if (ssp_vjust[i] == -1) {
    p <- p + annotation_custom(grob=grid::linesGrob(arrow=grid::arrow(type='closed', ends='first', length=unit(1.5, 'mm')), gp=gpar(fill='black')),
                               xmax=ssp_cum_melt[i, 'value'], xmin=ssp_cum_melt[i, 'value'] - 200*ssp_hjust[i], 
                               ymin=1, ymax=1.07 - 0.04*ssp_vjust[i])
  } else {
    p <- p + annotation_custom(grob=grid::linesGrob(arrow=grid::arrow(type='closed', ends='first', length=unit(1.5, 'mm')), gp=gpar(fill='black')),
                             xmax=ssp_cum_melt[i, 'value'], xmin=ssp_cum_melt[i, 'value'] - 200*ssp_hjust[i], 
                             ymin=1, ymax=1.09 - 0.04*ssp_vjust[i])
  }
}

p <- p + annotation_custom(grob=grid::linesGrob(arrow=grid::arrow(type='open', ends='first', length=unit(3, 'mm')), gp=gpar(col='blue')),
                           xmin=500, xmax=2000, ymin=1.22, ymax=1.22)
p <- p + annotation_custom(grob=grid::textGrob(label=expression('Increasing Probability of Achieving 2'*degree*C~'Target'), gp=gpar(col='blue', fontsize=8)),
                           xmin=1000, xmax=1500, ymin=1.27, ymax=1.27)


p <- p + coord_cartesian(clip = "off")
  
pdf('figures/co2_cum-cdf.pdf', width=3.5, height=4.1)
print(p)
dev.off()

png('figures/co2_cum-cdf.png', width=89, height=100, units='mm', res=600)
print(p)
dev.off()
