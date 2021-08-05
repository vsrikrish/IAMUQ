library(ggplot2)
library(reshape2)

nsamp <- 1e5

yrs <- seq(2020, 2100, 10)

scenarios <- c('base', 'nopen')

for (s in scenarios) {
    sim_out <- readRDS(paste0('output/sim_', s, '-gwp-co2-pop.rds'))
    tech_share_cols <- c('Frac_FossilHi', 'Frac_FossilLo', 'Frac_NonFossil')
    tech_share <- lapply(sim_out, function(l) l$out[l$out$year %in% yrs, c('year', tech_share_cols)])
    tech_share_melt <- melt(tech_share, id.vars='year')
    tech_share_melt$year <- factor(tech_share_melt$year, ordered=TRUE)

    labels <- c(Frac_FossilHi = 'High-Emissions Fossil Tech',
                Frac_FossilLo = 'Low-Emissions Fossil Tech',
                Frac_NonFossil = 'Non-Emitting Tech')

    p <- ggplot(tech_share_melt) + 
        geom_violin(aes(x=year, y=value)) +
        scale_x_discrete('Year') +
        scale_y_continuous('Global Energy Technology Share', labels=scales::percent) +
        facet_wrap(vars(variable), labeller = labeller(variable = labels), ncol=1, strip.position='right') +
        theme_classic(base_size = 12)

    png(paste0('figures/tech_share_violin-', s, '.png'), height=8, width=6, res=300, units='in')
    print(p)
    dev.off()

    pdf(paste0('figures/tech_share_violin-', s, '.pdf'), height=8, width=6)
    print(p)
    dev.off()
}
