% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_conditional_independence.R
\name{local_conditional_covariance}
\alias{local_conditional_covariance}
\title{Calculate the local conditional covariance between two variables}
\usage{
local_conditional_covariance(clg_object, coord = c(1, 2))
}
\arguments{
\item{clg_object}{The object produced by the clg-function}

\item{coord}{The variables for which the conditional covariance should be extracted}
}
\description{
Wrapper for the \code{clg} function that extracts the local Gaussian conditional
covariance between two variables from an object that is produced by the clg-function.
}
\details{
This function is a wrapper for the clag-function, and extracts the estimated local
conditional covariance between the first two variables in the data matrix, on the grid
specified to the clg-function.
}
