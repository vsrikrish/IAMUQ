library(ggplot2)
library(grid)
library(gridExtra)

model_out <- readRDS('output/mle-iid.rds')

# compute residuals for each output
r <- list()
pop <- merge(model_out[['dat']][['pop']], model_out[['bestfit']][, c('year', 'P')], by='year')
r[['pop']] <- data.frame(year=pop$year, residuals=pop$pop - pop$P)
prod <- merge(model_out[['dat']][['prod']], model_out[['bestfit']][, c('year', 'Q')], by='year')
r[['prod']] <- data.frame(year=prod$year, residuals=prod$prod - prod$Q)
emis <- merge(model_out[['dat']][['emissions']], model_out[['bestfit']][, c('year', 'C')], by='year')
r[['emissions']] <- data.frame(year=emis$year, residuals=emis$emissions - emis$C)

th <- theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

theme_set(th)

p <- list()
p[[1]] <- ggplot(pop) + geom_point(aes(x=year, y=pop)) + geom_line(aes(x=year, y=P), color='red') + ggtitle('Population') + scale_x_continuous('Year') + scale_y_continuous('Billions')
p[[2]] <- ggplot(prod) + geom_point(aes(x=year, y=prod)) + geom_line(aes(x=year, y=Q), color='red')+ ggtitle('GWP') + scale_x_continuous('Year') + scale_y_continuous('Trillions 2011US$')
p[[3]] <- ggplot(emis) + geom_point(aes(x=year, y=emissions)) + geom_line(aes(x=year, y=C), color='red') + ggtitle('Carbon Emissions') + scale_x_continuous('Year') + scale_y_continuous('Gt C')
p[[4]] <- ggplot(r[['pop']][-1,]) + geom_point(aes(x=year, y=residuals)) + ggtitle('Residuals') + scale_x_continuous('Year') + scale_y_continuous('Billions')
p[[5]] <- ggplot(r[['prod']][-1,]) + geom_point(aes(x=year, y=residuals)) + ggtitle('Residuals') + scale_x_continuous('Year') + scale_y_continuous('Trillions 2011US$')
p[[6]] <- ggplot(r[['emissions']]) + geom_point(aes(x=year, y=residuals)) + ggtitle('Residuals')  + scale_x_continuous('Year') + scale_y_continuous('Gt C')

pdf('figures/iid-resid.pdf')
do.call("grid.arrange", c(p, ncol=3))
dev.off()