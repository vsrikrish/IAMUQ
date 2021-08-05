library(IAMUQ)
library(readxl)
library(reshape2)
library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)
library(dplyr)

## read in CO2 emissions from the SSP database
ssp_dat <- read_excel('data/cmip6_co2.xlsx')
ssp_dat <- as.data.frame(ssp_dat)[1:7,]
ssp_dat[, 'Scenario'] <- gsub('\\s.*', '', ssp_dat[, 'Scenario'])
ssp_dat[, 'Scenario'] <- gsub('([a-z]*)(\\d)$', '\\1\\.\\2', ssp_dat[,'Scenario'])
ssp_dat <- ssp_dat[, c(2, 6:(ncol(ssp_dat)-2))]

co2_melt <- melt(ssp_dat, id.vars=c('Scenario'))
colnames(co2_melt)[2] <- 'Year'
co2_melt[, 'Year'] <- as.numeric(levels(co2_melt[, 'Year']))[co2_melt[, 'Year']]
co2_melt[, 'value'] <- co2_melt[, 'value'] / 1000
co2_melt$variable <- 'C'

## read in population from the SSP database
ssp_dat <- read_excel('data/ssp_pop.xlsx')
ssp_dat <- as.data.frame(ssp_dat)[1:7,]
ssp_dat[, 'Scenario'] <- gsub('\\s.*', '', ssp_dat[, 'Scenario'])
ssp_dat[, 'Scenario'] <- gsub('([a-z]*)(\\d)$', '\\1\\.\\2', ssp_dat[,'Scenario'])
ssp_dat <- ssp_dat[, c(2, 6:(ncol(ssp_dat)-1))]

pop_melt <- melt(ssp_dat, id.vars=c('Scenario'))
colnames(pop_melt)[2] <- 'Year'
pop_melt[, 'Year'] <- as.numeric(levels(pop_melt[, 'Year']))[pop_melt[, 'Year']]
pop_melt$variable <- 'P'
pop_melt$value <- pop_melt$value / 1000 # convert to billion persons

## read in GDP from the SSP database
ssp_dat <- read_excel('data/ssp_gdp.xlsx')
ssp_dat <- as.data.frame(ssp_dat)[1:7,]
ssp_dat[, 'Scenario'] <- gsub('\\s.*', '', ssp_dat[, 'Scenario'])
ssp_dat[, 'Scenario'] <- gsub('([a-z]*)(\\d)$', '\\1\\.\\2', ssp_dat[,'Scenario'])
ssp_dat <- ssp_dat[, c(2, 6:(ncol(ssp_dat)-1))]

gdp_melt <- melt(ssp_dat, id.vars=c('Scenario'))
colnames(gdp_melt)[2] <- 'Year'
gdp_melt[, 'Year'] <- as.numeric(levels(gdp_melt[, 'Year']))[gdp_melt[, 'Year']]
gdp_melt[, 'value'] <- gdp_melt$value * 1.15 # adjust for inflation between 2005 and 2011
gdp_melt$variable <- 'Q'
gdp_melt$value <- gdp_melt$value / 1000 # convert to trillions USD

## combine the SSP data into one melted dataframe
ssp_melt <- do.call(rbind, list(co2_melt, pop_melt, gdp_melt))

yrs <- 2020:2100

dat <- lapply(iamdata, function(l) {l[l$year %in% 1820:2019,]})
obs <- lapply(dat, function(l) l[l$year %in% 2000:2019,])

tol9qualitative=c("#88CCEE", "#44AA99", "#117733", "#332288", "#AA4499",  "#CC6677", "#882255", "#6699CC", "#999933")

cbbpsqualitative <- c("#000000", "#e79f00", "#9ad0f3", "#CC79A7", "#0072B2", "#009E73", "#F0E442", "#D55E00")

scenarios <- c('base', 'low', 'high', 'del_zc')
scen_labels <- c('Standard', 'Low Fossil Fuel', 'High Fossil Fuel', 'Delayed Zero-Carbon')

get_output <- function(sim_out, var, yrs) {
  do.call(cbind, lapply(sim_out, function(l) l$out[l$out$year %in% yrs, var]))
}

out_q <- function(out, yrs) {
  q <- data.frame(Year=yrs, t(apply(out, 1, quantile, probs=c(0.05, 0.95))))
  colnames(q)[2:3] <- c('lower', 'upper')
  
  q
}

out_dist <- function(out, yr) {
  data.frame(value=out[which(yrs == yr), ])
}

sim_out <- vector('list', length(scenarios))
names(sim_out) <- scenarios

for (scen in scenarios) {
  # read simulation output file
  sim_out[[scen]] <- readRDS(paste0('output/sim_', scen, '-gwp-co2-pop.rds'))
}

output_names <- c('P', 'Q', 'C')
output_labels <- c('Billion Persons', 'Trillions 2011US$', expression(CO[2]~Emissions~(Gt~CO[2]/yr)))
output_ylims <- data.frame(P=c(6, 13), Q=c(0, 750), C=c(0, 120))

p <- vector('list', length(output_names))

for (i in 1:length(output_names)) {
    out <- lapply(sim_out, get_output, var=output_names[i], yrs=yrs)

    q <- lapply(out, out_q, yrs=yrs)
    q <- bind_rows(q, .id='scenario')
    q$case <- factor(q$scenario, levels=scenarios, labels=scen_labels)

    marg_2100 <- lapply(out, out_dist, yr=2100)
    marg_2100 <- bind_rows(marg_2100, .id='scenario')
    marg_2100$case <- factor(marg_2100$scenario, levels=scenarios, labels=scen_labels)

    p_series <- ggplot() + geom_ribbon(data=q, 
                                   aes(x=Year, ymin=lower, ymax=upper, fill=case), color=NA, alpha=0.3) + 
        geom_line(data=ssp_melt[ssp_melt$variable == output_names[i], ], aes(x=Year, y=value, color=Scenario)) +  
        geom_point(data=obs[[i]], aes(x=year, y=value), color='black', size=1) + 
        scale_y_continuous(output_labels[i], limits=as.numeric(output_ylims[, i]), expand=c(0, 0), sec.axis=sec_axis(~., labels=NULL)) + 
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
        labs(tag=letters[i])

    p_marg <- ggplot() + stat_density(data=marg_2100, 
                                      aes(x=value, fill=case, color=case), geom='line', 
                                      position='identity') + 
        scale_x_continuous(limits=as.numeric(output_ylims[, i]), expand=c(0, 0)) + 
        scale_y_continuous(expand=c(0, 0)) + 
        coord_flip() + 
        scale_fill_manual('', values=cbbpsqualitative) + 
        scale_color_manual('', values=cbbpsqualitative) + 
        theme_classic(base_size=5) +
        theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line = element_blank(),
              panel.border=element_blank(), plot.margin=unit(c(11, 1, 13.2, -0.3), 'mm'), 
              axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), 
              axis.line.x=element_blank(), legend.position='none')

    p[[i]] <- gtable_row('scenarios', 
                    grobs=list(ggplotGrob(p_series + theme(legend.position='none')), 
                               ggplotGrob(p_marg)), height=unit(3, 'in'), widths=unit(c(2.7, 0.8), 'in'), 
                    z=c(2,1))
}

# extract legend from last series plot
g <- ggplot_gtable(ggplot_build(p_series))
legidx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
legend <- g$grobs[[legidx]]

fig <- arrangeGrob(p[[1]], p[[2]], p[[3]], legend, nrow=2, ncol=2)

pdf('figures/model_projections.pdf', height=6, width=7)
grid.draw(fig)
dev.off()

png('figures/model_projections.png', height=6, width=7, units='in', res=600)
grid.draw(fig)
dev.off()
