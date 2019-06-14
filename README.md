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

## Contact (corresponding author)
Vivek Srikrishnan
Email: <vivek@psu.edu>
