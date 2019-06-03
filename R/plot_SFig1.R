library(ggplot2)
library(grid)
library(gridExtra)
library(BAUcalib)

pars <- readRDS('output/map_iid.rds')$optim$bestmem
parnames <- names(pars)

model_out <- run_model(pars, parnames)

dat <- lapply(baudata, function(l) {l[l$year %in% 1820:2014,]})

res <- residuals(pars, parnames, model_out, dat)
res <- lapply(setNames(names(res), names(res)), function(n) {data.frame(year=dat[[n]]$year, residuals=res[[n]])})

r <- do.call(cbind, lapply(res, function(l) l$residuals))
colnames(r) <- c('Population', 'GWP', 'Emissions')
r_acf <- acf(r, plot=F)
print(r_acf)

pop <- merge(dat$pop, model_out[,c('year', 'P')], by='year')
prod <- merge(dat$prod, model_out[,c('year', 'Q')], by='year')
emis <- merge(dat$emissions, model_out[,c('year', 'C')], by='year')

th <- theme_bw(base_size=10) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

theme_set(th)

p <- list()
p[[1]] <- ggplot(pop) + geom_point(aes(x=year, y=pop, color='Observations')) + geom_line(aes(x=year, y=P, color='Best Fit')) + ggtitle('Population') + scale_x_continuous('Year') + scale_y_continuous('Billions') + scale_color_manual(name='', values=c('Observations'='black', 'Best Fit'='red'), breaks=c('Observations', 'Best Fit'), guide=guide_legend(override.aes = list(linetype=c(NA, 1), shape=c(16, NA)))) + theme(legend.position=c(0.4, 0.8))
p[[2]] <- ggplot(prod) + geom_point(aes(x=year, y=prod, color='Observations')) + geom_line(aes(x=year, y=Q, color='Best Fit'))+ ggtitle('Gross World Product') + scale_x_continuous('Year') + scale_y_continuous('Trillions 2011US$')  + scale_color_manual(name='', values=c('Observations'='black', 'Best Fit'='red'), breaks=c('Observations', 'Best Fit'), guide=guide_legend(override.aes = list(linetype=c(NA, 1), shape=c(16, NA)))) + theme(legend.position=c(0.4, 0.8))
p[[3]] <- ggplot(emis) + geom_point(aes(x=year, y=emissions, color='Observations')) + geom_line(aes(x=year, y=C, color='Best Fit')) + ggtitle('Carbon Emissions') + scale_x_continuous('Year') + scale_y_continuous('GtC/yr') + scale_color_manual(name='', values=c('Observations'='black', 'Best Fit'='red'), breaks=c('Observations', 'Best Fit'), guide=guide_legend(override.aes = list(linetype=c(NA, 1), shape=c(16, NA)))) + theme(legend.position=c(0.4, 0.8))
p[[4]] <- ggplot(res[['pop']]) + geom_point(aes(x=year, y=residuals)) + ggtitle('Residuals') + scale_x_continuous('Year') + scale_y_continuous('log Billions')
p[[5]] <- ggplot(res[['prod']]) + geom_point(aes(x=year, y=residuals)) + ggtitle('Residuals') + scale_x_continuous('Year') + scale_y_continuous('log Trillions 2011US$')
p[[6]] <- ggplot(res[['emissions']]) + geom_point(aes(x=year, y=residuals)) + ggtitle('Residuals')  + scale_x_continuous('Year') + scale_y_continuous('log GtC/yr')

pdf('figures/iid-resid.pdf', height=6.5, width=6.5)
do.call("grid.arrange", c(p, ncol=3))
dev.off()

png('figures/iid-resid.png', height=6.5, width=6.5, units='in',res=600)
do.call("grid.arrange", c(p, ncol=3))
dev.off()

