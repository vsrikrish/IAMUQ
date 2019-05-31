library(readxl)

data_path <- 'inst/extdata'
output_path <- 'data'
if (!dir.exists(output_path)) {
  dir.create(output_path)
}

econpop_file <- file.path(data_path, 'mpd2018.xlsx')
# read in per capita historical GDP (in 2011US$)
# we use cgadppc instead of rgdpnapc as we are interested in global GDP
cgdppc_dat <- read_excel(econpop_file, sheet = 3, skip = 1)
# read in population (in thousands)
pop_dat <- read_excel(econpop_file, sheet = 5, skip = 1)
# sum up world population figures in billions
po_dat <- data.frame(year=pop_dat[, 1], value = rowSums(pop_dat[ ,-1], na.rm = TRUE) / 1e6)
# estimate global GDP for each year (in trillions 2011US$)
pr_dat <- data.frame(year=pop_dat[, 1], value=rowSums(pop_dat[, -1] * cgdppc_dat[, -1], na.rm = TRUE) / 1e9)

# subset data based on valid years for the pop/GDP data
# we dont really trust anything before 1820, but include data going back to 1700 as a baseline and for validation
valid_yr <- 1700:2017
po_dat <- po_dat[po_dat$year %in% valid_yr, ]
pr_dat <- pr_dat[pr_dat$year %in% valid_yr, ]

# read in global emissions data (in Gt CO2)
em_dat <- read.table(file.path(data_path, 'global.1751_2014.ems'), skip = 35, colClasses = c(rep('numeric', 2), rep('NULL', 6)), fill = TRUE)
em_dat[, 2] <- em_dat[, 2]/1000
colnames(em_dat) <- c('year', 'value')

iamdata <- list(pop=po_dat, prod=pr_dat, emissions=em_dat)

devtools::use_data(iamdata, overwrite=TRUE)