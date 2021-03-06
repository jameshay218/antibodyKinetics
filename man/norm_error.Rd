% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{norm_error}
\alias{norm_error}
\title{Discretised normal error}
\usage{
norm_error(actual, obs, sd, MAX_TITRE)
}
\arguments{
\item{actual}{the assumed true titre}

\item{obs}{the observed titre}

\item{sd}{standard deviation of the observation error function}

\item{MAX_TITRE}{the maximum observable titre}
}
\description{
Gives the probability of a titre observation given a true titre.
}
