% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jrSiCKLSNMF.R
\name{PlotLossvsLatentFactors}
\alias{PlotLossvsLatentFactors}
\title{Create plots to help determine the number of latent factors to use for jrSiCKLSNMF}
\usage{
PlotLossvsLatentFactors(
  SickleJr,
  rounds = 100,
  differr = 1e-04,
  d_vector = c(2:20),
  parallel = FALSE,
  nCores = detectCores() - 1,
  subsampsize = NULL,
  minibatch = FALSE,
  random = FALSE,
  random_W_updates = FALSE,
  seed = NULL,
  batchsize = dim(SickleJr@count.matrices[[1]])[2],
  lossonsubset = FALSE,
  losssubsetsize = dim(SickleJr@count.matrices[[1]])[2],
  numthreads = 1
)
}
\arguments{
\item{SickleJr}{An object of class SickleJr}

\item{rounds}{Number of rounds to use: defaults to 100; this process is time consuming,
so a high number of rounds is not recommended}

\item{differr}{Tolerance for the percentage update in the likelihood: for these plots,
this defaults to \eqn{1e-4}}

\item{d_vector}{Vector of \eqn{d} values to test: default is 2 to 20}

\item{parallel}{Boolean indicating whether to use parallel computation}

\item{nCores}{Number of desired cores; defaults to the number of cores of the current machine minus 1 for convenience}

\item{subsampsize}{Size of the random subsample (defaults to \code{NULL}, which means all cells will be used); using a random subsample decreases computation time but sacrifices accuracy}

\item{minibatch}{Boolean indicating whether to use the mini-batch algorithm: default is \code{FALSE}}

\item{random}{Boolean indicating whether to use random initialization to generate the \eqn{\mathbf{W}^v} matrices and \eqn{\mathbf{H}} matrix:
defaults to \code{FALSE}}

\item{random_W_updates}{Boolean parameter for mini-batch algorithm; if \code{TRUE}, only updates \eqn{\mathbf{W}^v} once per epoch on the
penultimate subset of \eqn{\mathbf{H}}; otherwise updates \eqn{\mathbf{W}^v} after every update of the subset of \eqn{\mathbf{H}}}

\item{seed}{Number representing the random seed}

\item{batchsize}{Desired batch size; do not use if using a subsample}

\item{lossonsubset}{Boolean indicating whether to calculate the loss on a subset rather than the full dataset; speeds up computation for larger datasets}

\item{losssubsetsize}{Number of cells to use for the loss subset; default is total number of cells}

\item{numthreads}{Number of threads to use if running in parallel}
}
\value{
An object of class SickleJr with a list of initialized \eqn{\mathbf{W}^v} matrices and an \eqn{\mathbf{H}} matrix
for each latent factor \eqn{d\in\{1,...,D\}} added to the \code{WHinitials} slot, a data frame holding relevant
values for plotting the elbow plot added to the \code{latent.factor.elbow.values} slot, diagnostic plots of the loss vs. the number of latent factors added to the \code{plots}
slot, and the cell indices used to calculate the loss on the subsample added to the \code{lossCalcSubSample} slot
}
\description{
Generate plots of the lowest achieved loss after a
pre-specified number of iterations (default 100)
of the jrSiCKLSNMF algorithm for each latent factor (defaults to 2:20). This operates similarly to a
scree plot, so please select a number of latent factors that corresponds to the
elbow of the plot. This method is not appropriate for larger sets of data (more than 1000 cells)
}
\examples{
SimSickleJrSmall@latent.factor.elbow.values<-data.frame(NULL,NULL)
SimSickleJrSmall<-PlotLossvsLatentFactors(SimSickleJrSmall,d_vector=c(2:6),rounds=10)
}
