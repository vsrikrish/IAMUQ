##==============================================================================
## sobol_functions.R
##
## Original codes by Calvin Whealton, Cornell University
## https://github.com/calvinwhealton/SensitivityAnalysisPlots
## and
## Perry Oddo, Penn State
##
## Modified (condensed - codes are almost entirely unchanged) for brevity by
## Tony Wong (twong@psu.edu). The only code that was changed is the plotting
## routine 'plotRadCon'; Tony added
## (1) inputs for the first-, total- and second-order % values used for legend
## (2) generate labels for first-, total- and second-order in legend
## (3) write legend in black ink instead of gray
## (4) include '%' sign and cut legend labels off at the decimal (whole numbers
##     only)
##
## Tony also modified the 'sig' test for significance to test for confidence
## interval bounds of the same sign (otherwise, 0 is in CI for sensitivity
## index) and greater than 1%.
##
## Some functions were added by Vivek Srikrishnan (vivek@psu.edu),
## such as sobol_func_eval, sobol_func_wrap, and co2_yr. map_range and
## sample_value are based on Tony Wong's workflow.
##=============================================================================



library(parallel)

# map between [0,1] and a bounded parameter range
map_range <- function(x, bdin, bdout) {
  bdout[1] + (bdout[2] - bdout[1]) * ((x - bdin[1]) / (bdin[2] - bdin[1]))
}

# sample from a Gaussian distribution near parameter values and map to [0, 1]
sample_value <- function(n, parvals, bw, bds) {
  samps <- rnorm(n, mean=parvals, sd=bw) # sample from distribution
  # if any samples are outside of the parameter bounds, resample
  resample <- (samps < bds[1]) | (samps > bds[2])
  while (any(resample)) {
      samps[resample] <- rnorm(sum(resample), mean=parvals[resample], sd=bw)
      resample <- (samps < bds[1]) | (samps > bds[2])
  }
  map_range(samps, bdin=bds, bdout=c(0, 1)) # return after mapping to [0, 1]
}

# evaluate function in parallel
temp_eval_par <- function(vals, parnames, par_bds) {
  export_names <- c('parnames', 'temp_2100')
  # map samples from [0, 1] to parameter values
  samps <- vapply(parnames, function(p) map_range(vals[, match(p, parnames)], bdin=c(0, 1), bdout=par_bds[[p]]), numeric(nrow(vals)))
  
  cl <- makeCluster(detectCores()) # start cluster
  clusterEvalQ(cl, library(BAUcalib))
  clusterExport(cl, export_names) # export function to evaluate
    
  ## apply function to samples
  out <- parApply(cl, samps, 1, temp_wrap, parnames=parnames)
  
  stopCluster(cl) # stop cluster
  
  out - mean(out) # return function output centered at zero, which is required
}

# wrapper for function evaluation
temp_wrap <- function(pars, parnames) {
  # simulate model
  model_out <- run_model(pars, parnames, start=1700, end=2100)
  # call summary function and return
  do.call(temp_2100, list(emis=model_out$C, yrs=model_out$year, tcre=pars[match('TCRE', parnames)]))
}

# compute temperature based on cumulative emissions
temp_2100 <- function(emis, yrs, tcre, start=2014, end=2100, anomaly=0.92) {
  cum_emis <- cum_co2(emis, yrs, start=start, end=end)
  cum_emis * tcre / 1000 + anomaly
}