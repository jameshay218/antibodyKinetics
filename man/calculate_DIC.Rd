% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyses.R
\name{calculate_DIC}
\alias{calculate_DIC}
\title{DIC}
\usage{
calculate_DIC(chain)
}
\arguments{
\item{chain}{the MCMC chain with a lnlike column and all columns (number of columns will be used}
}
\value{
a single DIC value
}
\description{
Calculates DIC from a given MCMC chain
}
