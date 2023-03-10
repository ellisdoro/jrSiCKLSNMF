% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jrSiCKLSNMF.R
\name{NormalizeCountMatrix}
\alias{NormalizeCountMatrix}
\title{NormalizeCountMatrix}
\usage{
NormalizeCountMatrix(SickleJr, frob = FALSE)
}
\arguments{
\item{SickleJr}{An object of type SickleJr with a count matrix. Users can choose to normalize using median
library size normalization or log(x+1) normalization}

\item{frob}{A Boolean. Set to TRUE if you want to perform log(x+1) normalization and FALSE for
median library size normalization as per Zhang 2017.}
}
\value{
A list of sparse, normalized matrices
}
\description{
Perform normalization for count data. If you are planning to use the Frobenius norm,
set frob=TRUE to log(x+1) normalize your count data. This step can be skipped for percentage data and
spectra data. You may also skip this if you would like to perform a different form
of normalization; however, please ensure that if using the Frobenius norm that all of the values are
divided by the sum of all of the values.
}
