library(readxl)

data_path <- 'inst/extdata'
output_path <- 'data'
if (!dir.exists(output_path)) {
  dir.create(output_path)
}

econpop_file <- file.path(data_path, 'mpd2020.xlsx')
# read in per capita historical GDP (in 2011US$)
gdppc_dat <- read_excel(econpop_file, sheet = 4, skip = 1)
# read in population (in thousands)
pop_dat <- read_excel(econpop_file, sheet = 5, skip = 1)
# sum up world population figures in billions
po_dat <- data.frame(year=pop_dat[, 1], value = rowSums(pop_dat[ ,-1], na.rm = TRUE) / 1e6)
# estimate global GDP for each year (in trillions 2011US$)
pr_dat <- data.frame(year=pop_dat[, 1], value=rowSums(pop_dat[, -1] * gdppc_dat[, -1], na.rm = TRUE) / 1e9)

### subset data based on valid years for the pop/GDP data
## we dont really trust anything before 1820, but include data going back to 1700 as a baseline and for validation
valid_yr <- 1700:2018
po_dat <- po_dat[po_dat$year %in% valid_yr, ]
pr_dat <- pr_dat[pr_dat$year %in% valid_yr, ]
### add in UN and World Bank data for 2019; these are just pulled from the relevant websites and not stored in a local file.
# there is some inconsistency here between the series, but they're pretty close
po_dat <- rbind(po_dat, c(year=2019, value=7.713468))
pr_dat <- rbind(pr_dat, c(year=2019, value=119.2541)) # adjusted for inflation from 2017 to 2011US$

### read in global emissions data (in Gt CO2)
## from Boden et al (2017)
boden_dat <- read.table(file.path(data_path, 'global.1751_2014.ems'), skip = 35, colClasses = c(rep('numeric', 2), rep('NULL', 6)), fill = TRUE)
colnames(boden_dat) <- c('year', 'value')
## from the Global Carbon Budget 2020
gcp_file <- file.path(data_path, 'Global_Carbon_Budget_2020v1.0.xlsx')
gcp_dat <- read_excel(gcp_file, sheet=3, skip=8)
gcp_dat <- as.data.frame(gcp_dat[, 1:2])
colnames(gcp_dat) <- c('year', 'value')

em_dat <- rbind(boden_dat[boden_dat$year %in% seq(1700, 1958),], gcp_dat)
em_dat[, 'value'] <- em_dat[, 'value'] * 3.667 / 1000 # convert from Mt C/yr to Gt CO2/yr

iamdata <- list(pop=po_dat, prod=pr_dat, emissions=em_dat)

usethis::use_data(iamdata, overwrite=TRUE)
