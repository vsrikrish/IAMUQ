library(IAMUQ)
library(MASS)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(readxl)
library(reshape2)

yrs <- 2015:2100

get_variable <- function(sim_out, var) {
  do.call(cbind, lapply(sim_out, function(l) l$out[,var]))
}

# read simulation output file
sim_out <- readRDS('output/sim_base-gwp-co2.rds')

pop <- get_variable(sim_out, var='P')
gwp <- get_variable(sim_out, var='Q')
emis <- get_variable(sim_out, var='C')

avg_pc_gwp <- avg_gwp_rate(pop, gwp, yrs=2015:2100, start=2018)*100
avg_c_intens <- avg_gwp_rate(gwp, emis, yrs=2015:2100, start=2018)*100
co2_cum <- cum_co2(emis, yrs=2015:2100, start=2018)

dat <- data.frame(GWP=avg_pc_gwp, CI=avg_c_intens, CO2=co2_cum)
dat$CI <- -1*dat$CI
co2col <- brewer.pal(9, 'Reds')

ssp_ci <- as.data.frame(read_excel('data/ssp_ci.xlsx')[c(1,3:5), c(2,8,16)])
ssp_ci$rate <- (exp((log(ssp_ci[,3]) - log(ssp_ci[,2]))/(2100-2020)) - 1) * 100
ssp_gwp <- as.data.frame(read_excel('data/ssp_pc_gwp.xlsx')[c(1,3:5), c(2,8,16)])
ssp_gwp$rate <- (exp((log(ssp_gwp[,3]) - log(ssp_gwp[,2]))/(2100-2020)) - 1) * 100
ssp_names <- c('SSP3-7.0', 'SSP4-6.0', 'SSP2-4.5', 'SSP5-8.5')
ssp_dat <- data.frame(SSP=ssp_names, CI=-1*ssp_ci$rate, GWP=ssp_gwp$rate)
ssp_dat$hoff <- c(0.2, 0, 0, 0)
ssp_dat$hpos <- ssp_dat$CI + ssp_dat$hoff
ssp_dat$voff <- c(-0.5, -0.5, -0.5, 0.5)
ssp_dat$vpos <- ssp_dat$GWP + ssp_dat$voff
ssp_dat$arrstart <- ssp_dat$vpos - ssp_dat$voff / 5
ssp_dat$arrend <- ssp_dat$GWP + ssp_dat$voff / 5
tol9qualitative=c("#117733", "#332288",  "#CC6677", "#44AA99")

# compute 95% contour
dat_kde <- kde2d(dat[,'CI'], dat[,'GWP'], n = 500)
dx <- diff(dat_kde$x[1:2]) 
dy <- diff(dat_kde$y[1:2])
sz <- sort(dat_kde$z)
c1 <- cumsum(sz) * dx * dy
prob <- 0.95
dimnames(dat_kde$z) <- list(dat_kde$x,dat_kde$y)
dc <- melt(dat_kde$z)
dc$prob <- approx(sz,1-c1,dc$value)$y
md <- dc[which.max(dc$value),]

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())

marg_gwp <- ggplot(dat) + stat_density(aes(x=GWP), geom='line', color='blue') +
  theme_classic(base_size= 6) + 
  scale_x_continuous(limits=c(-1.5, 5.5), expand=c(0, 0)) +
  coord_flip() +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), 
        axis.line=element_blank(), legend.position='none',
        plot.margin=unit(c(0, 0, 23.5, 0), 'mm'))

marg_ci <- ggplot(dat) + stat_density(aes(x=CI), geom='line', color='blue') +
  theme_classic(base_size=6) + 
  scale_x_continuous(limits=c(-1.5, 5.5), expand=c(0, 0)) +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), 
        axis.line=element_blank(), legend.position='none',
        plot.margin=unit(c(0, 0, 0, 8.5), 'mm'))

scatplot <- ggplot(dat[dat$CO2 < 3000,], aes(x=CI, y=GWP)) + 
  stat_summary_2d(aes(z=CO2), fun=mean, binwidth=0.075) +
  geom_contour(data=dc, aes(x=Var1, y=Var2, z=prob), breaks=prob, color='grey', alpha=0.8) +
  geom_point(data=md, aes(x=Var1, y=Var2), color='grey', size=2) +
  geom_text(data=md, aes(x=Var1, y=Var2+0.5), color='black', label='Posterior Mode', size=3) +
  geom_segment(data=md, aes(x=Var1, xend=Var1, y=Var2+0.4, yend=Var2+0.1),
               size=1, arrow=arrow(length=unit(0.5, 'mm'), type='closed'), 
               color='black', arrow.fill='black') +
  geom_abline(slope=1, intercept=0, linetype='dotted') +
  geom_point(data=ssp_dat, aes(x=CI, y=GWP), color='black', size=2) +
  geom_text(data=ssp_dat, aes(x=hpos, y=vpos, label=SSP), size=3) +
  geom_segment(data=ssp_dat, aes(x=hpos, xend=CI+c(0.05, 0, 0, 0), y=arrstart, yend=arrend),
               size=1, arrow=arrow(length=unit(0.5, 'mm'), type='closed')) +
  theme_classic(base_size=10) +
  scale_x_continuous('Rate of Carbon Intensity Decrease 2018-2100 (%)', limits=c(-1.5, 5.5), expand=c(0, 0)) +
  scale_y_continuous('Rate of Per-Capita GWP Growth 2018-2100 (%)', limits=c(-1.5, 5.5), expand=c(0, 0)) +
  scale_fill_gradientn(expression('Cumulative'~CO[2]~'Emissions 2018-2100 (GtC)'), 
                       colors=co2col) +
  theme(legend.position='bottom', 
        legend.key.height=unit(2, 'mm'),
        legend.box='vertical', legend.box.just='left') +
  annotate('text', x=0, y=4.2, color='orangered4', size=3, label='Runaway Emissions') +
  annotate('text', x=4.75, y=3.25, color='green4', size=3, label='Approximate\nStabilization') 


fig <- arrangeGrob(marg_ci, empty, scatplot, marg_gwp, nrow=2, ncol=2,
                   widths=c(8, 1), heights=c(1, 8))

png('figures/growth-scatter.png', height=120, width=120, units='mm', res=300)
grid.draw(fig)
dev.off()

pdf('figures/growth-scatter.pdf', height=4.7, width=4.7)
grid.draw(fig)
dev.off()
