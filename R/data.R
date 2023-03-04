#' A simulated dataset with U(1,1.25) multiplicative noise for the
#' scRNA-seq variability parameter in SPARSim for the simulated scRNA-seq data and
#' with N(-0.25,0.25) additive noise to the expression levels of the scATAC-seq data
#' for data simulated via SimATAC. The simulated matrices are located in SimData$Xmatrices
#' and the identities for the cell types are contained in SimData$cell_type. This corresponds
#' to the Xmatrix data found in both XandLmatrices25/XandindividLKNNLmatrices1Sparsity5.RData
#' and XandBulkLmatrix25/XandBulkLKNNmatrices1Sparsity5.RData.on our Github
#' ellisdoro/jrSiCKLSNMF_Simulations
#' @name SimData
#'
#' @docType data
#'
#' @usage data(SimData)
#'
#' @format A list made up of a list of 2 simulated matrices and a vector containing cell identities.
#'
#' @keywords datasets
#' @source \href{https://github.com/ellisdoro/jrSiCKLSNMF_Simulations}{jrSicKLSNMF Simulations}

NULL
