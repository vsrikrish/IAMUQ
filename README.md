# IAMUQ

## Overview

IAMUQ (Integrated Assessment Model Uncertainty Quantification) defines a simple integrated assessment model (IAM) and provides scripts for calibration and simulation. This code is intended to accompany

> Srikrishnan, V,  Guan, Y, Tol, RSJ, and Keller, K, 2021. "Probabilistic projections of baseline 21st century CO$_2$ emissions using a simple calibrated integrated assessmentmModel."

This study uses a simple, mechanistically-motivated model, calibrated on century-scale historical data, to generate probabilistic projections of anthropogenic CO<sub>2</sub> emissions through 2100 under several different scenarios. These scenarios are related to deep uncertainties concerning available fossil fuel resources and the rate of global decarbonization. We find that the fossil fuel supply contraints can affect the upper tail of the emissions distribution. These projections make no assumption about future decarbonization or negative emissions trends beyond a continuation of dynamics consistent with historical and current trends, and so considerations of changes to climate policy or technological backsliding will have an effect on the projections.

We also conduct a global sensitivity analysis and find that the interactions between terms, particularly involving economic terms affecting total factor productivity growth and capital stock dynamics, are particularly important for explaining the variability in cumulative emissions for the rest of the century. Another very important factor is the emissions intensity of our stylized fossil technology which roughly corresponds to natural gas and oil.

This code (in `Rcpp`) is intended to facilitate replication or expansion of these results. Code is included for an `IAMUQ` R package, which includes the core functionality for the model (written in `Rcpp`), calibration (using Markov Chain Monte Carlo), and simulation (both written in R). Other scripts are included for calibrating specific scenarios and generating figures.

## Requirements

### `IAMUQ` package installation

The `Rcpp` code (which utilizes RcppArmadillo to access the armadillo linear algebra routines) requires an updated version of GCC. GCC versions before 5 may not work, while GCC 7 definitely does. 

The following R packages are required (and should be installed when the package is installed):

1. `Rcpp` (to compile the model code and link to R)
2. `RcppArmadillo` (for armadillo linear algebra routines and views in Rcpp)
3. `parallel` (for calibration using multiple cores)
4. `DEoptim` (we find a maximum a posteriori (MAP) estimate using differential evolution)
5. `adaptMCMC` (the model is calibrated using adaptive MCMC)
6. `mvtnorm` (a model discrepancy structure results in a multivariate normal likelihood function)
7. `magic` (for constructing block-diagonal matrices).

### Analysis and plotting

The other R codes require the following packages:

1. `parallel` (for simulation across multiple cores)
2. `readr` (to import data files)
3. `dplyr` (to bind list elements into a data frame)
4. `reshape2` (to melt data into long format)
5. `ggplot2` (for plotting)
6. `gridExtra` (additional plot layout functions)
7. `gtable` (additional plot layout functions)
8. `RColorBrewer` (color-blind color palettes)
9. `plotrix` (plotting, used for the Sobol` sensitivity plot)
10. `GGally` (for the ggpairs function)

## Installation

To install the `IAMUQ` package, first clone the repository. Then, use the following commands from the local repository root directory:

1. `R CMD build IAMUQ`
2. `R CMD INSTALL IAMUQ_1.0.tar.gz`

## Data

The original population, GDP, and CO<sub>2</sub> emissions data can be found in `IAMUQ/inst/extdata`. Code to process this data is provided in `IAMUQ/data-raw/read_data.R`. This code produces `IAMUQ/data/iamdata.rda`, which is loaded by the `IAMUQ` package. The raw population and GDP per capita data are taken from:

> Bolt, Jutta and Jan Luiten van Zanden (2020), Maddison Project Database, version 2020. [“Maddison style estimates of the evolution of the world economy. A new 2020 update ”](https://www.rug.nl/ggdc/historicaldevelopment/maddison/publications/wp15.pdf). URL: `<https://www.rug.nl/ggdc/historicaldevelopment/maddison/releases/maddison-project-database-2020?lang=en> (accessed 07-01-2020)`.

To extend this data to 2019, we adjust 2019 GDP PPP from the [World Bank](https://data.worldbank.org/indicator/NY.GDP.MKTP.PP.KD) to 2011 international $ and use 2019 population data from the [United Nations](https://population.un.org/wpp/).

CO<sub>2</sub> emissions data from 1700--1949 are taken from:

> Boden, T. A., Andres, R. J., and Marland, G. 2017. "Global, Regional, and National Fossil-Fuel CO2 Emissions (1751 - 2014) (V. 2017)". doi: `<https://doi.org/10.3334/CDIAC/00001_V2017>`.

Emissions data from 1950--2019 are from:

> Friedlingstein, P., O’Sullivan, M., Jones, M.W., Andrew, R.M., Hauck, J., Olsen, A., Peters, G.P., Peters, W., Pongratz, J., Sitch, S., Le Quéré, C., Canadell, J.G., Ciais, P., Jackson, R.B., Alin, S., Aragão, L.E.O.C., Arneth, A., Arora, V., Bates, N.R., Becker, M., Benoit-Cattin, A., Bittig, H.C., Bopp, L., Bultan, S., Chandra, N., Chevallier, F., Chini, L.P., Evans, W., Florentie, L., Forster, P.M., Gasser, T., Gehlen, M., Gilfillan, D., Gkritzalis, T., Gregor, L., Gruber, N., Harris, I., Hartung, K., Haverd, V., Houghton, R.A., Ilyina, T., Jain, A.K., Joetzjer, E., Kadono, K., Kato, E., Kitidis, V., Korsbakken, J.I., Landschützer, P., Lefèvre, N., Lenton, A., Lienert, S., Liu, Z., Lombardozzi, D., Marland, G., Metzl, N., Munro, D.R., Nabel, J.E.M.S., Nakaoka, S.-I., Niwa, Y., O’Brien, K., Ono, T., Palmer, P.I., Pierrot, D., Poulter, B., Resplandy, L., Robertson, E., Rödenbeck, C., Schwinger, J., Séférian, R., Skjelvan, I., Smith, A.J.P., Sutton, A.J., Tanhua, T., Tans, P.P., Tian, H., Tilbrook, B., van der Werf, G., Vuichard, N., Walker, A.P., Wanninkhof, R., Watson, A.J., Willis, D., Wiltshire, A.J., Yuan, W., Yue, X., Zaehle, S., 2020. "Global Carbon Budget 2020". *Earth System Science Data* 12, 3269–3340. doi: `<https://doi.org/10.5194/essd-12-3269-2020>`.

CMIP6 Representative Concentration Pathway - Shared Socioeconomic Pathway scenario data was downloaded from the [SSP database](https://tntcat.iiasa.ac.at/SspDb/dsd). The resulting data is available in `data/cmip6_co2.xlsx`.

We also use distributions provided or derived by three expert assessments and probabilistic projections. Those are:

> Christensen, P., Gillingham, K., Nordhaus, W., 2018. "Uncertainty in forecasts of long-run economic growth". *Proc Natl Acad Sci USA* 115, 5409–5414. doi: `<https://doi.org/10.1073/pnas.1713628115>`.
> Ho, E., Budescu, D.V., Bosetti, V., van Vuuren, D.P., Keller, K., 2019. "Not all carbon dioxide emission scenarios are equally likely: a subjective expert assessment." *Climatic Change* 155, 545–561. doi: `<https://doi.org/10.1007/s10584-019-02500-y>`.
> United Nations, Department of Economic and Social Affairs, Population Division, 2019. Probabilistic Population Projections, Rev. 1, based on the World Population Prospects 2019, Rev. 1. URL: `<https://population.un.org/wpp/Download/Standard/Population/>`

## Workflow

The calibration and analysis codes are designed to work with a PBS job submission queue on a high-performance computing system and should be easily modifiable to work with a different submission system. To run individual cases without calibration, the maximum a posteriori (MAP) and Markov chain Monte Carlo (MCMC) can be used by passing options through the command line (see `R/map_scenarios.R` and `R/mcmc_scenarios.R`).

### Model calibration

The calibrations in the manuscript can be reproduced using scenarios and priors set up in the repository code (in `R/map_scenarios.R`, `R/mcmc_scenarios.R`, and `R/calib_priors.R`). The calibrations under alternate priors are set up in `R/map_sensitivity.R` and `R/mcmc_sensitivity.R`. Other calibration scenarios can be set up using similar method)s for specifying priors, data, and expert assessment inclusion.

By default, the MAP and MCMC analyses are run using multiple processors. The IAMUQ package contains the `IAMUQ::find_map()` and `IAMUQ::run_mcmc()` functions, which can be run in serial instead (see documentation).

### Projections

To generate projections based on a calibration, use the `IAMUQ::cond_sim_model()` function, as in `R/conditional_simulation.R`.

### Sobol' sensitivity analysis

The global sensitivity analysis is conducted using `R/sobol.R`. This code writes out text files with the parameter sensitivities and confidence intervals.

### Plotting

The code to generate the plots in Srikrishnan, Guan, Tol, and Keller (2021) is available in the repository:

1. Figure 1: `R/plot_hindcast.R`
2. Figure 2: `R/plot_projections_CMIP6.R`
3. Figure 3: `R/plot_cdf.R`
4. Figure 4: `R/plot_hindcast_length.R`
5. Figure 5: R/sobol_plot.R
6. Figures S1 and S6: `R/plot_tech_shares.R`
5. Figure S2: `R/plot_post_scenario.R`
6. Figure S3: `R/plot_projections_CMIP6_all.R`
7. Figure S4: `R/plot_estimate_inversion.R` (this figure uses only the base model scenario, though the script will generate figures for the others as well).
8. Figure S5: `R/plot_projections_CMIP6_techconstraint.R` 
9. Figure S7: `R/plot_growth_rates.R`
10. Figure S8: `R/plot_post_prior.R`
11. Figure S9: `R/plot_post_calibration.R`
12. Figure S10: `R/plot_pairs.R`
13. Figure S11: `R/plot_proj_sensitivity.R`

## Contact (corresponding author)
Vivek Srikrishnan
Email: <viveks@cornell.edu>
