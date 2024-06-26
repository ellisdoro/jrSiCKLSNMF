% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{jrSiCKLSNMF}
\alias{jrSiCKLSNMF}
\title{Run jrSiCKLSNMF outside of a SickleJr object}
\usage{
jrSiCKLSNMF(
  datamatL,
  WL,
  H,
  AdjL,
  DL,
  lambdaWL,
  lambdaH,
  initsamp,
  suppress_warnings,
  diffFunc,
  Hconstraint = "None",
  Wconstraint = "None",
  differr = 1e-06,
  rounds = 1000L,
  display_progress = TRUE,
  minibatch = TRUE,
  batchsize = 100L,
  random_W_updates = TRUE,
  minrounds = 100L,
  numthreads = 1L
)
}
\arguments{
\item{datamatL}{An R list where each entry contains a normalized, sparse
\eqn{\mathbf{X}^v} matrix corresponding to single-cell modality \eqn{v}}

\item{WL}{An R list containing initialized values of the \eqn{\mathbf{W}^v}
within each modality \eqn{v}}

\item{H}{A matrix containing initialized values for the shared
\eqn{\mathbf{H}}}

\item{AdjL}{An R list containing all of the adjacency matrices for the
feature-feature similarity graphs in sparse format;  note that
\eqn{\mathtt{D-Adj}} is the graph Laplacian}

\item{DL}{An R list containing all of the degree matrices of the
feature-feature similarity graphs;  note that \eqn{\mathtt{D-Adj}} is
the graph Laplacian}

\item{lambdaWL}{A list of the \eqn{\lambda_{\mathbf{W}^v}} corresponding to
modality \eqn{v}}

\item{lambdaH}{A double containing the desired value for
\eqn{\lambda_{\mathbf{H}}}}

\item{initsamp}{A vector of randomly selected rows of \eqn{\mathtt{H}} on
which to run the objective function}

\item{suppress_warnings}{A Boolean that indicates whether warnings
should be suppressed}

\item{diffFunc}{A string indicating what type of divergence to use;  set to
the Poisson Kullback-Leibler divergence
(\eqn{\mathtt{``klp"}}) by default, but the Frobenius norm
(\eqn{\mathtt{``fr"}}) is also available}

\item{Hconstraint}{A string that indicates whether you want to set an L2
norm constraint on the rows of \eqn{\mathbf{H}}. Enter `None` for
no constraints or `L2Norm` to set the L2 norm of each row of
\eqn{\mathbf{H}} to 1, or `Ortho` for an orthogonality constraint on the
rows of \eqn{\mathbf{H}}}

\item{Wconstraint}{A string that indicates whether you want to set an L2
norm constraint on the columns of \eqn{\mathbf{W}^v}. Enter `None` for
no constraints, `L2Norm` to set the L2 norm of each row of
\eqn{\mathbf{W}^v} to 1, or `Ortho` for an orthogonality constraint on the
rows of \eqn{\mathbf{W}^v}}

\item{differr}{A double containing the tolerance}

\item{rounds}{A double containing the number of rounds}

\item{display_progress}{A Boolean indicating whether to display the progress
bar}

\item{minibatch}{A Boolean indicating whether to use the mini-batch version
of the algorithm}

\item{batchsize}{Number of batches for mini-batch updates}

\item{random_W_updates}{A Boolean indicating whether to update
\eqn{\mathbf{W}^v} once per epoch (TRUE) or after every update of
the subset of \eqn{\mathbf{H}} (FALSE) for the mini-batch algorithm.}

\item{minrounds}{A minimum number of rounds for the algorithm to run: most
useful for the mini-batch algorithm}

\item{numthreads}{Number of threads to use if running in parallel}
}
\value{
An R list containing values for the objective function.
}
\description{
Perform joint non-negative matrix factorization (NMF) across
multiple modalities of single-cell data.
To measure the discrepancy between two distributions, one can use the
Poisson Kullback-Leibler divergence
(\eqn{\mathtt{diffFunc}=}\eqn{\mathtt{"klp"}}),
the Itakura-Saito divergence (\eqn{\mathtt{diffFunc}=}\eqn{\mathtt{"is"}}),
or the Frobenius norm (\eqn{\mathtt{diffFunc}=}\eqn{\mathtt{"fr"}}).
It is also possible to set graph regularization constraints on
\eqn{\mathbf{W}^v} and either a sparsity constraint on \eqn{\mathbf{H}} or
an L2 norm constraint on the rows of \eqn{\mathbf{H}}.
This function passes by reference and updates the variables
\eqn{\mathtt{WL}} and \eqn{\mathtt{H}} and does not require data to be in an
object of type SickleJr. \eqn{\mathtt{RunjrSiCKLSNMF}} calls this function.
If your data are in an object of class SickleJr,
please use the \eqn{\mathtt{RunjrSiCKLSNMF}} function instead.
}
