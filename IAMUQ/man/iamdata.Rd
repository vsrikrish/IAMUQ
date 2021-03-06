% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{iamdata}
\alias{iamdata}
\title{Calibration data for the IAM, consisting of global population, gross world
 product, and CO2 emissions.}
\format{
A list with three elements, one for each data component. Each list
 element is a data frame with two columns, 'year' and 'value.'
 \describe{
   \item{pop}{global mean population, in billions}
   \item{prod}{gross world product, in trillions US$2011}
   \item{emissions}{CO2 emissions, in GtC/yr}
 }
}
\source{
Population and gross world product are taken from the 2020
 Maddison Project dataset, \url{ https://www.rug.nl/ggdc
 /historicaldevelopment/maddison/releases/maddison-project-database-2020}.
These are then updated to include 2019 GDP data from the World Bank, \url{https://data.worldbank.org/indicator/NY.GDP.MKTP.PP.KD}
and population data from the UN, \url{https://population.un.org/wpp}.
 CO2 emissions are taken from Boden et al (2017), \url{https://www.osti.gov/dataexplorer
 /biblio/1389331-global-regional-national-fossil-fuel-co2-emissions} from 1751-1958, and
from the Global Carbon Project 2020, \url{https://dx.doi.org/10.18160/gcp-2020} from 1959-2019.
}
\usage{
iamdata
}
\description{
A dataset containing estimates of global population (from 1700-2019),
 gross world product (from 1700-2019), and CO2 emissions (from 1751-2019).
 The CO2 emissions neglect emissions related to land use change.
}
\keyword{datasets}
