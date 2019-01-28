library(ggplot2)
library(grid)
library(gridExtra)

# compute residuals for each output
model_out <- readRDS('temp-ar.rds')
a <- model_out[['mle']]$optim$bestmem[match(c('a_pop', 'a_prod', 'a_emis'), model_out[['parnames']])] # AR coefficient
names(a) <- c('pop', 'prod', 'emissions')

r <- list()
pop <- merge(model_out[['dat']][['pop']], model_out[['bestfit']][, c('year', 'P')], by='year')
res <- data.frame('model'= pop$P[-1] - a['pop']*pop$P[-length(pop$P)],
                  'obs' = pop$pop[-1] - a['pop']*pop$pop[-length(pop$pop)]
)
r[['pop']] <- data.frame(year=pop$year[-1], residuals=res$model - res$obs)
prod <- merge(model_out[['dat']][['prod']], model_out[['bestfit']][, c('year', 'Q')], by='year')
res <- data.frame('model'= prod$Q[-1] - a['prod']*prod$Q[-length(prod$Q)],
                  'obs' = prod$prod[-1] - a['prod']*prod$prod[-length(prod$prod)]
)
r[['prod']] <- data.frame(year=prod$year[-1], residuals=res$model - res$obs)
emis <- merge(model_out[['dat']][['emissions']], model_out[['bestfit']][, c('year', 'C')], by='year')
res <- data.frame('model'= emis$C[-1] - a['emissions']*emis$C[-length(emis$C)],
                  'obs' = emis$emissions[-1] - a['emissions']*emis$emissions[-length(emis$emissions)]
)
r[['emissions']] <- data.frame(year=emis$year[-1], residuals=res$model - res$obs)

th <- theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

theme_set(th)

p <- list()
p[[1]] <- ggplot(pop) + geom_point(aes(x=year, y=pop)) + geom_line(aes(x=year, y=P), color='red') + ggtitle('Population') + scale_x_continuous('Year') + scale_y_continuous('Billions')
p[[2]] <- ggplot(prod) + geom_point(aes(x=year, y=prod)) + geom_line(aes(x=year, y=Q), color='red')+ ggtitle('GWP') + scale_x_continuous('Year') + scale_y_continuous('Trillions 2011US$')
p[[3]] <- ggplot(emis) + geom_point(aes(x=year, y=emissions)) + geom_line(aes(x=year, y=C), color='red') + ggtitle('Carbon Emissions') + scale_x_continuous('Year') + scale_y_continuous('Gt C')
p[[4]] <- ggplot(r[['pop']]) + geom_point(aes(x=year, y=residuals)) + ggtitle('Residuals') + scale_x_continuous('Year') + scale_y_continuous('Billions')
p[[5]] <- ggplot(r[['prod']]) + geom_point(aes(x=year, y=residuals)) + ggtitle('Residuals') + scale_x_continuous('Year') + scale_y_continuous('Trillions 2011US$')
p[[6]] <- ggplot(r[['emissions']]) + geom_point(aes(x=year, y=residuals)) + ggtitle('Residuals')  + scale_x_continuous('Year') + scale_y_continuous('Gt C')

pdf('figures/ar-resid.pdf')
do.call("grid.arrange", c(p, ncol=3))
dev.off()
