##==========================================================================
## sobol_plot.R
##
## Sobol sensitvity analysis for business-as-usual emissions
## --> plotting routines
##
## Code history:
## this one by Vivek Srikrishnan, 30 April 2019, Penn State. Modified a bit
## relative to the previous versions listed below.
###################################
## Adapted from 'BRICK_Sobol_plotting.R'
## Original authored by Tony Wong
## Pennsylvania State University
## twong@psu.edu
#######################################
## Adapted from 'radialPlot_vanDantzig.R'
## Originally authored by: Perry Oddo
## Pennsylvania State University
## poddo@psu.edu
###################################
## Adapted from 'radialConvergeTest.R'
## Originally authored by: Calvin Whealton
## Cornell University
## caw324@cornell.edu
####################################
## Code for radial Sobol Analysis plot
## Original code available at:
## https://github.com/calvinwhealton/SensitivityAnalysisPlots
####################################

rm(list = ls())

n_params <- 17 # set number of parameters
# set files with sobol indices
sobol_file_1 <- paste0('output/Sobol-1-tot_temp.txt')
sobol_file_2 <- paste0('output/Sobol-2-tot_temp.txt')

## plot spider plot
library(RColorBrewer)
library(graphics)
library(plotrix)

source('R/sobol_plot_functions.R')
source('R/colorBlindPalette.R')

## Import data from sensitivity analysis
# First- and total-order indices
s1st <- read.csv(sobol_file_1,
                  sep=' ',
                  header=TRUE,
                  nrows = n_params,
                  as.is=c(TRUE,rep(FALSE,5)))
  
parnames <- s1st[,1]

# Import second-order indices
s2_table <- read.csv(sobol_file_2,
               sep=' ',
               header=TRUE,
               nrows = n_params*(n_params-1)/2,
               as.is=c(TRUE,rep(FALSE,4)))

# Convert second-order to upper-triangular matrix
s2 <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2)
s2 <- as.data.frame(s2)
colnames(s2) <- rownames(s2) <- s1st$Parameter

# Convert confidence intervals to upper-triangular matrix
s2_conf_low <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2_conf_high <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2_conf_low[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf_low)
s2_conf_high[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf_high)

s2_conf_low <- as.data.frame(s2_conf_low)
s2_conf_high <- as.data.frame(s2_conf_high)
colnames(s2_conf_low) <- rownames(s2_conf_low) <- s1st$Parameter
colnames(s2_conf_high) <- rownames(s2_conf_high) <- s1st$Parameter

# Determine which indices are statistically significant

sig.cutoff <- 0.01

# S1 & ST: using the confidence intervals
s1st1 <- stat_sig_s1st(s1st
                      ,method="congtr"
                      ,greater=sig.cutoff
                      ,sigCri='either')

# S1 & ST: using greater than a given value
#s1st1 <- stat_sig_s1st(s1st
#                      ,method="gtr"
#                      ,greater=0.01
#                      ,sigCri='either')

# S2: using the confidence intervals
s2_sig1 <- stat_sig_s2(s2
                       ,s2_conf_low
                       ,s2_conf_high
                       ,method='congtr'
                       ,greater=sig.cutoff
                       )

# S2: using greater than a given value
#s2_sig1 <- stat_sig_s2(s2
#                       ,s2_conf
#                       ,greater=0.02
#                       ,method='gtr')

## Define groups for the variables and the color schemes
# Defining lists of the variables for each group
name_list <- list('Population \n Growth' = parnames[1:4],
                  'Economic \n Output' = parnames[5:11],
                  'Emissions & \n Technology' = parnames[12:17]
                 )
                 
# add Parameter symbols to plot
name_symbols <- c(expression(psi[1]), expression(psi[2]),
                  expression(psi[3]), expression(P[0]),
                  expression(lambda), 's', expression(delta),
                  expression(alpha), expression(A[s]),
                  expression(pi), expression(A[0]),
                  expression(rho[2]), expression(rho[3]),
                  expression(tau[2]), expression(tau[3]),
                  expression(tau[4]), expression(kappa)
                 )
                 
source('R/colorBlindPalette.R')

# defining list of colors for each group
col_list <- list("Population \n Growth"     = rgb(mycol[11,1],mycol[11,2],mycol[11,3])
                  ,'Economic \n Output' = rgb(mycol[7,1],mycol[7,2],mycol[7,3]),
                   'Emissions & \n Technology'   = rgb(mycol[1,1],mycol[1,2],mycol[1,3])
                  )
                  
# using function to assign variables and colors based on group
s1st1 <- gp_name_col(name_list
                     ,col_list
                     ,s1st1)

s1st1$symbols <- name_symbols

# set filename for plot
plot_file <- paste0('figures/sobol_cum_co2')

plotRadCon(df=s1st1
           ,s2=s2
           ,scaling = .4
           ,s2_sig=s2_sig1
           ,filename = plot_file
           ,plotType = 'EPS'
           ,gpNameMult=1.6
           ,RingThick=0.1
           ,legLoc = "bottomcenter",cex = .76
           ,s1_col = rgb(mycol[3,1],mycol[3,2],mycol[3,3])
           ,st_col = rgb(mycol[6,1],mycol[6,2],mycol[6,3])
           ,line_col = rgb(mycol[10,1],mycol[10,2],mycol[10,3])
           ,STthick = 0.5
           ,legFirLabs=c(.05,.77), legTotLabs=c(.05,.83), legSecLabs=c(.02,.05)
           ,lsetback=FALSE
)

## Further analysis for the text:
##

# what are the highest first-order indices?
s1.sort <- s1st[rev(order(s1st[,'S1'])),1:4]
itmp <- which(s1.sort[,'S1'] > sig.cutoff & s1.sort[,'S1_conf_low']*s1.sort[,'S1_conf_high'] > 0)
s1.sort <- s1.sort[itmp,]
print('********************************')
print('significant first-order indices:')
print(s1.sort)
print('********************************')

# what are the highest total-order indices?
st.sort <- s1st[rev(order(s1st[,'ST'])),c(1,5:7)]
itmp <- which(st.sort[,'ST'] > sig.cutoff & st.sort[,'ST_conf_low']*st.sort[,'ST_conf_high'] > 0)
st.sort <- st.sort[itmp,]
print('********************************')
print('significant total-order indices:')
print(st.sort)
print('********************************')

# what are the highest second-order interaction indices?
s2.sort <- s2_table[rev(order(s2_table[,3])),]
itmp <- which(s2.sort[,'S2'] > sig.cutoff & s2.sort[,'S2_conf_low']*s2.sort[,'S2_conf_high'] > 0)
s2.sort <- s2.sort[itmp,]
print('********************************')
print('significant second-order indices:')
print(s2.sort)
print('********************************')

<- s2.sort[itmp,]
print('********************************')
print('significant second-order indices:')
print(s2.sort)
print('********************************')

