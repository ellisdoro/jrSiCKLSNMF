% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jrSiCKLSNMF.R
\name{PlotSickleJrUMAP}
\alias{PlotSickleJrUMAP}
\title{Generate UMAP plots for an object of class SickleJr}
\usage{
PlotSickleJrUMAP(
  SickleJr,
  umap.modality = "H",
  cluster = "kmeans",
  title = "",
  colorbymetadata = NULL,
  legendname = NULL
)
}
\arguments{
\item{SickleJr}{An object of class SickleJr}

\item{umap.modality}{String corresponding to the name of the UMAP of interest: defaults to \code{"H"}}

\item{cluster}{String input that indicates which cluster to color by: defaults to \code{"kmeans"}}

\item{title}{String input for optional plot title}

\item{colorbymetadata}{Name of metadata column if coloring by metadata}

\item{legendname}{String input that to allow specification of a different legend name}
}
\value{
An object of class SickleJr with plots added to the \code{plots} slot
}
\description{
Plot the first and second dimensions of a UMAP dimension reduction
and color either by clustering results or metadata.
}
\examples{
SimSickleJrSmall<-PlotSickleJrUMAP(SimSickleJrSmall,
title="K-Means Example")
SimSickleJrSmall<-PlotSickleJrUMAP(SimSickleJrSmall,umap.modality=1)
}
