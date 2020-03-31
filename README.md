# IAMUQ

## Overview

IAMUQ (Integrated Assessment Model Uncertainty Quantification) defines a simple integrated assessment model (IAM) and provides scripts for calibration and simulation. This code is intended to accompany 
>Srikrishnan, V,  Guan, Y, Tol, RSJ, and Keller, K (2019): Available fossil fuel resources, decarbonization rates, and economic growth  are key drivers of feasibility to achieve Paris climate targets.

This study analyzes projections of anthropogenic CO<sub>2</sub> emissions through 2100 under several different scenarios. These scenarios are related to deep uncertainties concerning available fossil fuel resources and the rate of global decarbonization. We find that these deep uncertainties are, along economic growth dynamics, the key socioeconomic drivers in determining the achievability of the Paris climate targets.

This code (in Rcpp) is intended to facilitate replication or expansion of these results. Code is included for an `IAMUQ` R package, which includes the core functionality for the model (written in Rcpp), calibration (using Markov Chain Monte Carlo), and simulation (both written in R). Other scripts are included for calibrating specific scenarios and generating figures.

## Requirements

### IAMUQ package installation

The Rcpp code (which utilizes RcppArmadillo to access the armadillo linear algebra routines) requires an updated version of GCC. GCC versions before 5 may not work, while GCC 7 definitely does. 

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

To install the `IAMUQ` package, use the following commands from the local repository root directory:

1. `R CMD build IAMUQ`
2. `R CMD INSTALL IAMUQ_1.0.tar.gz`

## Data

The original Maddison Project and Boden et al (2017) data is available in IAMUQ/inst/extdata. Code to process this data is provided in IAMUQ/data-raw/read_data.R. This code produces IAMUQ/data/iamdata.rda, which is loaded by the `IAMUQ` package.

CMIP6 Representative Concentration Pathway - Shared Socioeconomic Pathway scenario data was downloaded from the [SSP database](https://tntcat.iiasa.ac.at/SspDb/dsd). The resulting data is available in data/cmip6_co2.xlsx.

## Workflow

The calibration and analysis codes are designed to work with a PBS job submission queue on a high-performance computing system and should be easily modifiable to work with a different submission system. To run individual cases without calibration, the maximum a posteriori (MAP) and Markov chain Monte Carlo (MCMC) can be used by passing options through the command line (see R/map_scenarios.R and R/mcmc_scenarios.R).

### Model calibration

The calibrations in the manuscript can be reproduced using scenarios and priors set up in the repository code (in R/map_scenarios.R, R/mcmc_scenarios.R, and R/calib_priors.R). The calibrations under alternate priors are set up in R/map_sensitivity.R and R/mcmc_sensitivity.R. Other calibration scenarios can be set up using similar methods for specifying priors, data, and expert assessment inclusion.

By default, the MAP and MCMC analyses are run using multiple processors. The IAMUQ package contains the `IAMUQ::find_map()` and `IAMUQ::run_mcmc()` functions, which can be run in serial instead (see documentation).

### Projections

To generate projections based on a calibration, use the `IAMUQ::cond_sim_model()` function, as in R/conditional_simulation.R.

### Sobol' sensitivity analysis

The global sensitivity analysis is conducted using R/sobol.R. This code writes out text files with the parameter sensitivities and confidence intervals.

### Plotting

The code to generate the plots in Srikrishnan, Guan, Tol, and Keller (2020) is available in the repository:

1. Figure 1: R/plot_projections_CMIP6.R
2. Figure 2: R/plot_cdf.R
3. Figure 3: R/sobol_plot.R
4. Extended Data Figure 1: R/plot_proj_sensitivity.R
5. Extended Data Figure 2: R/plot_hindcast_length.R
6. Extended Data Figure 3: R/plot_posterior_prior.R
7. Extended Data Figure 4: R/plot_post_scenario.R
8. Extended Data Figure 5: R/plot_estimate_inversion.R (this figure uses only the base model scenario, though the script will generate figures for the others as well).
9. Extended Data Figure 6: R/plot_projections_CMIP6_all.R
10. Extended Data Figure 7: R/plot_growth_rates.R

## Contact (corresponding author)
Vivek Srikrishnan
Email: <vivek@psu.edu>
