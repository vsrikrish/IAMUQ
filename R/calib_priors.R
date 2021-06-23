#########################################################
# set_prior_params(parnames):                           #
#   This function sets up a data frame of prior         #
#   parameters with more "intuitive" parameters such as #
#   confidence interval bounds. Those are then          #
#   converted to the parameters required by the density #
#   functions when necessary.                           #
#########################################################

set_prior_params <- function(parnames) {
  
  # initialize prior parameter data frame
  prior_df <- data.frame(name=parnames, type=character(length(parnames)), lower=numeric(length(parnames)), upper=numeric(length(parnames)), stringsAsFactors=FALSE)
  
  # loop over data frame rows and fill in values
  for (i in 1:nrow(prior_df)) {
    name <- prior_df[i, 'name']
    if (grepl('sigma', name)) {
      prior_df[i, 'type'] <- 'log-normal'
      prior_df[i, 'lower'] <- qlnorm(0.05, -1, 1)
      prior_df[i, 'upper'] <- qlnorm(0.95, -1, 1)
    } else if (grepl('eps', name)) {
      prior_df[i, 'type'] <- 'log-normal'
      prior_df[i, 'lower'] <- qlnorm(0.05, -1, 1)
      prior_df[i, 'upper'] <- qlnorm(0.95, -1, 1)
    } else if (name == 'psi1') {
      prior_df[i, 'type'] <- 'normal'
      prior_df[i, 'lower'] <- 0.0001
      prior_df[i, 'upper'] <- 0.15
    } else if (name == 'psi2') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0
      prior_df[i, 'upper'] <- 50
    } else if (name == 'psi3') {
      prior_df[i, 'type'] <- 'normal'
      prior_df[i, 'lower'] <- 6.9
      prior_df[i, 'upper'] <- 14.4
    } else if (name == 'P0') {
      prior_df[i, 'type'] <- 'normal'
      prior_df[i, 'lower'] <- 0.3
      prior_df[i, 'upper'] <- 0.9
    } else if (name == 'lambda') {
      prior_df[i, 'type'] <- 'normal'
      prior_df[i, 'lower'] <- 0.6
      prior_df[i, 'upper'] <- 0.8
    } else if (name == 's') {
      prior_df[i, 'type'] <- 'normal'
      prior_df[i, 'lower'] <- 0.22
      prior_df[i, 'upper'] <- 0.26
    } else if (name == 'delta') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0.01
      prior_df[i, 'upper'] <- 0.14
    } else if (name == 'alpha') {
      prior_df[i, 'type'] <- 'normal'
      prior_df[i, 'lower'] <- 0.0007
      prior_df[i, 'upper'] <- 0.0212
    } else if (name == 'As') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 5.3
      prior_df[i, 'upper'] <- 16.11
    } else if (name == 'pi') {
      prior_df[i, 'type'] <- 'normal'
      prior_df[i, 'lower'] <- 0.62
      prior_df[i, 'upper'] <- 0.66
    } else if (name == 'A0') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0
      prior_df[i, 'upper'] <- 3
    } else if (name == 'rho2') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0
      prior_df[i, 'upper'] <- 0.75
    } else if (name == 'rho3') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0
      prior_df[i, 'upper'] <- 0.75
    } else if (name == 'tau2') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 1700
      prior_df[i, 'upper'] <- 2100
    } else if (name == 'tau3') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 1700
      prior_df[i, 'upper'] <- 2100
    } else if (name == 'tau4') {
      prior_df[i, 'type'] <- 'truncnorm'
      prior_df[i, 'lower'] <- 2050
      prior_df[i, 'upper'] <- 2150
    } else if (name == 'kappa') {
      prior_df[i, 'type'] <- 'uniform'
      prior_df[i, 'lower'] <- 0.005
      prior_df[i, 'upper'] <- 0.2
    } else if (grepl('a_', name)) {
      lab <- sub('.*_', '', name)
      if (is.na(as.numeric(lab))) {
        prior_df[i, 'type'] <- 'uniform'
        prior_df[i, 'lower'] <- 0.5
        prior_df[i, 'upper'] <- 1
      } else {
        idx <- as.numeric(unlist(strsplit(lab, '*')))
        if (idx[1] == idx[2]) {
          prior_df[i, 'type'] <- 'normal'
          prior_df[i, 'lower'] <- 0
          prior_df[i, 'upper'] <- 1
        } else {
          prior_df[i, 'type'] <- 'normal'
          prior_df[i, 'lower'] <- -1
          prior_df[i, 'upper'] <- 1
        }
      }
    }
  }
  
  prior_df
}
