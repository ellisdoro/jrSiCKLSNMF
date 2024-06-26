% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jrSiCKLSNMF.R
\name{SetLambdasandRegs}
\alias{SetLambdasandRegs}
\title{Set lambda values and type of row regularization for an object of class SickleJr}
\usage{
SetLambdasandRegs(
  SickleJr,
  lambdaWlist = list(10, 50),
  lambdaH = 500,
  HrowReg = "None",
  WcolReg = "None"
)
}
\arguments{
\item{SickleJr}{An object of class SickleJr}

\item{lambdaWlist}{A list of graph regularization constraints for the \eqn{\mathbf{W}^v} matrices:
defaults to 2 modalities with the RNA modality constraint equal to 10 and the ATAC modality constraint equal to 50}

\item{lambdaH}{A numeric holding the sparsity constraint on \eqn{\mathbf{H}}: defaults to 500.}

\item{HrowReg}{A string that is equal to \code{"None"} for no constraints on the rows of \eqn{\mathbf{H}},
\code{"Ortho"} for orthogonality constraints on the rows of \eqn{\mathbf{H}}, and \code{"L2Norm"}
to set the L2 norms of the rows of \eqn{\mathbf{H}} to be equal to 1: defaults to "None"}

\item{WcolReg}{A string that is equal to \code{"None"} for no constraints on the columns of \eqn{\mathbf{W}},
\code{"Ortho"} for orthogonality constraints on the columns of \eqn{\mathbf{W}}, and
\code{"L2Norm"} to set the L2 norms of the columns of \eqn{\mathbf{W}} to be equal to 1: defaults to "None"}
}
\value{
An object of class SickleJr with the lambda hyperparameter values added to its \code{lambdaWlist} and \code{lambdaH} slots
}
\description{
Provide the values for the graph regularization \eqn{\lambda_{\textbf{W}^v}}
for each modality as a list and provide a
}
\examples{
SimSickleJrSmall<-SetLambdasandRegs(SimSickleJrSmall,
lambdaWlist=list(10,50),lambdaH=500,HrowReg="None",WcolReg="None")
SimSickleJrSmall<-SetLambdasandRegs(SimSickleJrSmall,
lambdaWlist=list(3,15),lambdaH=0,HrowReg="L2Norm",WcolReg="None")
}
