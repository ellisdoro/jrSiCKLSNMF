% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jrSiCKLSNMF.R
\name{DetermineClusters}
\alias{DetermineClusters}
\title{Perform clustering diagnostics}
\usage{
DetermineClusters(
  SickleJr,
  numclusts = 2:20,
  clusteringmethod = "kmeans",
  diagnosticmethods = c("wss", "silhouette", "gap_stat"),
  clValidvalidation = "internal",
  createDiagnosticplots = TRUE,
  runclValidDiagnostics = TRUE,
  printPlots = TRUE,
  printclValid = TRUE,
  subset = FALSE,
  subsetsize = 1000,
  seed = NULL
)
}
\arguments{
\item{SickleJr}{An object of class SickleJr}

\item{numclusts}{A vector of integers indicating the number of clusters to test}

\item{clusteringmethod}{String holding the clustering method: defaults to k-means; since the other methods
are not implemented in jrSiCKLSNMF, it is recommended to use k-means.}

\item{diagnosticmethods}{Vector of strings indicating which methods to plot. Defaults to all three of the available: wss, silhouette, and gap_stat}

\item{clValidvalidation}{String containing validation method to use for clValid. Defaults to internal.}

\item{createDiagnosticplots}{Boolean indicating whether to create diagnostic plots for cluster size}

\item{runclValidDiagnostics}{Boolean indicating whether to calculate the diagnostics from \code{clValid}}

\item{printPlots}{Boolean indicating whether to print the diagnostic plots}

\item{printclValid}{Boolean indicating whether to print the diagnostic results from clValid}

\item{subset}{Boolean indicating whether to calculate the diagnostics on a subset of the data rather
than on the whole dataset.}

\item{subsetsize}{Numeric value indicating size of the subset}

\item{seed}{Numeric value holding the random seed}
}
\value{
An object of class SickleJr with cluster diagnostics added to its \code{clusterdiagnostics} slot
}
\description{
A wrapper for the \code{\link[clValid]{clValid}} and \code{\link[factoextra]{fviz_nbclust}} functions to perform clustering diagnostics
}
\examples{
#Since these data are too small, the clValid diagnostics do not run
#properly. See the vignette for an example with the clValid diagnostics
SimSickleJrSmall<-DetermineClusters(SimSickleJrSmall,numclusts=2:5,runclValidDiagnostics=FALSE)
}
