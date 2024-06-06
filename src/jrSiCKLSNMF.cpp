// -*- mode: C++;  c-indent-level: 4;  c-basic-offset: 4;
// indent-tabs-mode: nil;  -*-
// Copyright 2024 Dorothy Ellis
#include <omp.h>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <vector>
#include <math.h>
#include "Utils.h"
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

arma::sp_mat col_sp(const arma::sp_mat& x, const arma::uvec& index) {
  int n_cols = index.n_elem;
  arma::sp_mat x_subset(x.n_rows, index.n_elem);
  for (int i = 0;  i < n_cols;  i++) {
    x_subset.col(i) = x.col(index(i));
  }
  return x_subset;
  }
arma::sp_mat row_sp(const arma::sp_mat& x, const arma::uvec& index) {
  int n_rows = index.n_elem;
  arma::sp_mat x_subset(index.n_elem, x.n_cols);
  for (int i=0;  i < n_rows;  i++) {
    x_subset.row(i) = x.row(index(i));
  }
  return x_subset;
}
void ErrorCheck(const std::vector<std::string>& diffFunc,
                const std::string& Hconstraint,
                const std::string& Wconstraint) {
  for (int i = 0; i < static_cast<int>(diffFunc.size()); i++) {
    if ((diffFunc[i] != "klp") & (diffFunc[i] != "fr") &
        (diffFunc[i] != "is"))
      throw std::invalid_argument("Please enter 'klp' for Kullback-Leibler "
                                    "divergence for count data, 'is' for the "
                                    "Itakura-Saito divergence for spectra data"
                                    ", or 'fr' to calculate the Frobenius norm "
                                    "of the difference between X and WH.");
    }
  if ((Hconstraint != "None") & (Hconstraint != "L2Norm") &
      (Hconstraint != "Ortho"))
    throw std::invalid_argument("Please enter \n 'None' for no row constraints "
                                  "on H, \n 'L2Norm' for an L2 Norm constraint "
                                  "(i.e. the square root of the sum of squares "
                                  "of each row entry equals 1) on the rows of "
                                  "H, or \n 'Ortho' for an orthogonality "
                                  "constraint on the rows of H");
  if ((Wconstraint != "None") & (Wconstraint != "Ortho") &
      (Wconstraint != "L2Norm"))
    throw std::invalid_argument("Please enter \n 'None' for no constraints on "
                                  "the columns of each W matrix, \n 'L2Norm' "
                                  "for an L2 Norm constraint () or \n 'Ortho' "
                                  "for an orthogonality constraint on the "
                                  "columns of each W matrix");
}
double lossmatcalc(const arma::sp_mat& datamat, const arma::mat& W,
                   const arma::mat& H, const arma::sp_mat& Adj,
                   const arma::sp_mat& D, const double lambdaW,
                   const double lambdaH, const std::string& diffFunc = "klp",
                   int numv = 2) {
  double lik = 0.0;
  arma::mat Ht = trans(H);
  arma::mat WH = W*Ht;
  if (diffFunc == "klp") {
    arma::mat tmpWH = WH;
    arma::mat tmpWHlog = tmpWH.transform([](double val) {
      return std::max(val, 1e-16); });
    arma::sp_mat logxtmp(datamat);
    logxtmp.for_each([](arma::sp_mat::elem_type& val) {val = log(val); });
    lik = lik+arma::as_scalar(arma::accu(datamat%logxtmp))-
      arma::as_scalar(arma::accu(datamat%log(tmpWHlog)))-
      arma::as_scalar(arma::accu(datamat))+arma::as_scalar(arma::accu(tmpWH));
  } else if (diffFunc == "fr") {
    arma::mat tmp = datamat-WH;
    double frobnorm2 = pow(arma::norm(tmp, "fro"), 2);
    lik = lik+0.5*frobnorm2;
  } else if (diffFunc == "is") {
    arma::mat tmpWH = WH;
    tmpWH = tmpWH.transform([](double val) {
      return std::max(val, 1e-16); });
    arma::mat invtmpWH = tmpWH.transform([](double val) {
      return 1/val; });
    arma::mat xtmp(datamat);
    xtmp = xtmp.transform([](double val) {
      return std::max(val, 1e-16); });
    lik = lik+arma::as_scalar(arma::accu(datamat%invtmpWH+1))-
      arma::as_scalar(arma::accu(log(xtmp)))+
      arma::as_scalar(arma::accu(log(tmpWH)));
  } else {
    perror("Please enter 'klp' for the Kullback-Leibler divergence "
             "for count data, 'is' for the Itakura-Saito divergence, "
             "or 'fr' for the Frobenius norm loss");
    return -1;
  }
  if (lambdaW > 0) {
    arma::mat Wt = trans(W);
    double traceWAW = arma::trace(Wt*(D-Adj)*W);
    lik = lik+0.5*lambdaW*traceWAW;
  }
  if (lambdaH > 0) {
    double numvterm = static_cast<double>(numv);
    numvterm = 1/numvterm;
    double hfrob = arma::as_scalar(arma::accu(arma::trace(H*Ht)));
    lik = lik+0.5*lambdaH*numvterm*hfrob;
  }
  return lik;
}

double losscalc(const arma::field<arma::sp_mat>& datamatF,
                const arma::field<arma::mat>& WF,
                const arma::mat& H, const arma::field<arma::sp_mat>& AdjF,
                const arma::field<arma::sp_mat>& DF,
                arma::vec lambdaWV, const double& lambdaH,
                const std::vector<std::string>& diffFunc) {
  double lik = 0.0;
  int numv(datamatF.n_rows);
  for (int i = 0; i < numv; i++) {
    double liki = lossmatcalc(datamatF[i], WF[i], H, AdjF[i], DF[i],
                              lambdaWV[i], lambdaH, diffFunc[i], numv);
    if (liki == -1) {
      break;
    } else {
      lik = lik+liki;
    }
  }
  return lik;
}

//' @title Calculate the loss function in R
//' @description This is a wrapper function that allows for the calculation of
//' the loss function in R code.
//' @name lossCalcRWrapper
//' @param datamatL An R list where each entry contains a normalized, sparse
//' \eqn{\mathbf{X}^v} matrix corresponding to single-cell modality \eqn{v}
//' @param WL An R list containing initialized values of the \eqn{\mathbf{W}^v}
//' within each modality \eqn{v}
//' @param H A matrix containing initialized values for the shared
//' \eqn{\mathbf{H}}
//' @param AdjL An R list containing all of the adjacency matrices for the
//' feature-feature similarity graphs in sparse format;  note that
//' \eqn{\mathtt{D-Adj}} is the graph Laplacian
//' @param DL An R list containing all of the degree matrices of the
//' feature-feature similarity graphs;  note that \eqn{\mathtt{D-Adj}} is the
//' graph Laplacian
//' @param lambdaWL A list of the \eqn{\lambda_{\mathbf{W}^v}} corresponding to
//' modality \eqn{v}
//' @param lambdaH A double containing the desired value for
//' \eqn{\lambda_{\mathbf{H}}}
//' @param diffFunc A vector of strings indicating what type of divergence to
//' use;  set to the Poisson Kullback-Leibler divergence
//' (\eqn{\mathtt{``klp"}}) by default, but the Frobenius norm
//' (\eqn{\mathtt{``fr"}}) is also available
//' @returns A double containing the value of the jrSiCKLSNMF loss function
//' @export
// [[Rcpp::export]]

double lossCalcRWrapper(const Rcpp::List& datamatL, const Rcpp::List& WL,
                        const arma::mat& H, const Rcpp::List& AdjL,
                        const Rcpp::List& DL, const Rcpp::List& lambdaWL,
                        const double& lambdaH,
                        const std::vector<std::string>& diffFunc) {
  int numv(datamatL.size());
  arma::field<arma::sp_mat> datamatF(numv, 1);
  arma::field<arma::mat> WF(numv, 1);
  arma::vec lambdaWV(numv);
  arma::field<arma::sp_mat> DF(numv, 1);
  arma::field<arma::sp_mat> AdjF(numv, 1);
  arma::field<arma::sp_mat> datamatFsub(numv, 1);
  int numrowsH = H.n_rows;
  arma::uvec numbercells = arma::linspace<arma::uvec>(0, numrowsH-1, numrowsH);
  arma::field<arma::uvec> numberfeats(numv, 1);
  for (int i = 0;  i < numv; i++) {
    datamatF[i] = Rcpp::as<arma::sp_mat>(datamatL[i]);
    WF[i] = Rcpp::as<arma::mat>(WL[i]);
    int numrowsW = WF[i].n_rows;
    numberfeats[i] = arma::linspace<arma::uvec>(0, numrowsW-1, numrowsW);
    lambdaWV[i] = Rcpp::as<double>(lambdaWL[i]);
    AdjF[i] = Rcpp::as<arma::sp_mat>(AdjL[i]);
    DF[i] = Rcpp::as<arma::sp_mat>(DL[i]);
    if ((sum(diagvec(DF[i])) == 0)&(lambdaWV[i] > 0)) {
      DF[i] = arma::speye(WF[i].n_rows, WF[i].n_rows);
      // Adjacency matrix is already 0s so no need to worry about it
    }
  }
  double lik = losscalc(datamatF, WF, H, AdjF, DF, lambdaWV, lambdaH, diffFunc);
  return lik;
}

void NMFinview(arma::mat& old, arma::mat denom, arma::mat numer) {
  numer = numer.transform([](double val) {return std::max(val, 1e-16); });
  denom = denom.transform([](double val) {return std::max(val, 1e-16); });
  denom = denom.transform([](double val) {return 1/val; });
  old = old%numer%denom;
}
// Perform tangent gradient as in Le Roux et al. 2015
void normalizeMat(arma::mat& Mat, const std::string& Matconstraint) {
  if (Matconstraint == "L2Norm") {
    Mat = normalise(Mat);
    return;
  } else {
    return;
  }
}
arma::mat regFunc(const arma::mat& denomnumer, const arma::mat& Matrix,
                  const std::string& Matconstraint) {
  arma::mat regularizernumerdenom(denomnumer.n_rows, denomnumer.n_cols,
                                  arma::fill::value(0.0));
  if (Matconstraint == "L2Norm") {
    arma::mat sumdenomnumer = arma::sum(Matrix%denomnumer);
    arma::mat regularizermat = arma::repmat(sumdenomnumer, Matrix.n_rows, 1);
    regularizernumerdenom = Matrix%regularizermat;
    return regularizernumerdenom;
  } else if (Matconstraint == "Ortho") {
    regularizernumerdenom = Matrix*(denomnumer.t()*Matrix);
    return regularizernumerdenom;
  } else {
    return regularizernumerdenom;
  }
}

void perviewNMFMUR(const arma::field<arma::sp_mat>& datamatF,
                   arma::field<arma::mat>& WF, arma::mat& H,
                   const arma::field<arma::sp_mat>& AdjF,
                   const arma::field<arma::sp_mat>& DF,
                   const arma::vec lambdaWV, const double lambdaH,
                   const std::vector<std::string>& diffFunc,
                   arma::mat& Hfull, const std::string& Hconstraint,
                   const std::string& Wconstraint,
                   const arma::uvec batchsub, bool random_W_updates = false,
                   int batchiter = -1, int batchsize = -1) {
  int numv(datamatF.n_rows);
  arma::mat numerH(H.n_rows, H.n_cols, arma::fill::value(0.0));
  arma::mat denomH(H.n_rows, H.n_cols, arma::fill::value(0.0));

  arma::mat InumerH(H.n_rows, H.n_cols, arma::fill::value(0.0));
  arma::mat IdenomH(H.n_rows, H.n_cols, arma::fill::value(0.0));
  // Adding in column regularization to improve identifiability
  // normalize H to have l2 norm equal to 1
  // new additions are all multiplied by HHones for reference later
  // Since H is non-negative, |H|=H and is differentiable for all H
  // arma::mat Hones= H.ones();

  // Now we need to divide each column by the squared L2 norm of the
  // corresponding column entries of HHones.Since each entry is 1, this is the
  // same as dividing by the the number of rows of HHones.
  // See Douglas, 2000;  Fu 2018.
  for (int i = 0; i < numv; i++) {
    if (diffFunc[i] == "fr") {
      arma::mat WW = WF[i].t()*WF[i];
      // Take gradient with normalized parameters
      IdenomH = H*WW;
      InumerH = datamatF[i].t()*WF[i];
    } else if (diffFunc[i] == "klp") {
      arma::mat WHinv = WF[i]*H.t();
      WHinv = WHinv.transform([](double val) {return std::max(val, 1e-16); });
      WHinv = WHinv.transform([](double val) {return 1/val; });
      arma::sp_mat XWHinv(datamatF[i]%WHinv);
      InumerH = XWHinv.t()*WF[i];
      IdenomH = repmat(sum(WF[i]), H.n_rows, 1);
    } else if (diffFunc[i] == "is") {
      arma::mat WHinv = WF[i]*H.t();
      WHinv = WHinv.transform([](double val) {return std::max(val, 1e-16); });
      WHinv = WHinv.transform([](double val) {return 1/val; });
      arma::mat WHinv2 = WHinv;
      WHinv2 = WHinv.transform([](double val) {return pow(val, 2); });
      arma::sp_mat XWHinv2(datamatF[i]%WHinv2);
      InumerH = XWHinv2.t()*WF[i];
      IdenomH = WHinv.t()*WF[i];
    }
    numerH = numerH+InumerH;
    denomH = denomH+IdenomH;
  }
  if (lambdaH > 0) {
    denomH = denomH+lambdaH*H;
  }
  if (Hconstraint != "None") {
    arma::mat numerHold = numerH;
    numerH = numerH+regFunc(denomH, H, Hconstraint);
    denomH = denomH+regFunc(numerHold, H, Hconstraint);
  }
  NMFinview(H, denomH, numerH);
  H = H.transform([](double val) {
    return(val < 1e-10)? static_cast<double>(0) : static_cast<double>(val); });
  if (!batchsub.empty()) {
    Hfull.rows(batchsub) = H;
    normalizeMat(Hfull, Hconstraint);
    H = Hfull.rows(batchsub);
  } else {
    normalizeMat(H, Hconstraint);
  }
  if (!random_W_updates || (batchiter%batchsize == batchsize-1)) {
    for (int i = 0; i < numv; i++) {
      arma::mat numerW(WF[i].n_rows, WF[i].n_cols, arma::fill::value(0.0));
      arma::mat denomW(WF[i].n_rows, WF[i].n_cols, arma::fill::value(0.0));
      arma::mat regnumerW(WF[i].n_rows, WF[i].n_cols, arma::fill::value(0.0));
      arma::mat regdenomW(WF[i].n_rows, WF[i].n_cols, arma::fill::value(0.0));
      if (diffFunc[i] == "fr") {
        arma::mat HH = H.t()*H;
        numerW = numerW+datamatF(i, 0)*H;
        denomW = denomW+WF[i]*HH;
      } else if (diffFunc[i] == "klp") {
        arma::mat WHinv = WF[i]*H.t();
        WHinv = WHinv.transform([](double val) {return std::max(val, 1e-16); });
        WHinv = WHinv.transform([](double val) {return 1/val; });
        arma::mat XWHinv(datamatF[i]%WHinv);
        numerW = numerW+XWHinv*H;
        denomW = denomW+repmat(sum(H), WF[i].n_rows, 1);
      } else if (diffFunc[i] == "is") {
        arma::mat WHinv = WF[i]*H.t();
        WHinv = WHinv.transform([](double val) {return std::max(val, 1e-16); });
        WHinv = WHinv.transform([](double val) {return 1/val; });
        arma::mat WHinv2 = WHinv;
        WHinv2 = WHinv2.transform([](double val) {return pow(val, 2); });
        arma::mat XWHinv2(datamatF[i]%WHinv2);
        numerW = numerW+XWHinv2*H;
        denomW = denomW+WHinv*H;
      }
      if (lambdaWV[i] > 0) {
        denomW = denomW+0.5*lambdaWV[i]*(DF[i]*WF[i]+DF[i].t()*WF[i]);
        numerW = numerW+0.5*lambdaWV[i]*(AdjF[i]*WF[i]+AdjF[i].t()*WF[i]);
      }
      if (Wconstraint != "None") {
        arma::mat numerWold = numerW;
        arma::mat denomWold = denomW;
        numerW = numerW+regFunc(denomWold, WF[i], Wconstraint);
        denomW = denomW+regFunc(numerWold, WF[i], Wconstraint);
      }
      NMFinview(WF[i], denomW, numerW);
      WF[i] = WF[i].transform([](double val) {
        return(val < 1e-10)? static_cast<double>(0) : static_cast<double>(val);
        });
      normalizeMat(WF[i], Wconstraint);
    }
  }
}
//' @title Run jrSiCKLSNMF outside of a SickleJr object
//' @description Perform joint non-negative matrix factorization (NMF) across
//' multiple modalities of single-cell data.
//' To measure the discrepancy between two distributions, one can use the
//' Poisson Kullback-Leibler divergence
//' (\eqn{\mathtt{diffFunc}=}\eqn{\mathtt{"klp"}}),
//' the Itakura-Saito divergence (\eqn{\mathtt{diffFunc}=}\eqn{\mathtt{"is"}}),
//' or the Frobenius norm (\eqn{\mathtt{diffFunc}=}\eqn{\mathtt{"fr"}}).
//' It is also possible to set graph regularization constraints on
//' \eqn{\mathbf{W}^v} and either a sparsity constraint on \eqn{\mathbf{H}} or
//' an L2 norm constraint on the rows of \eqn{\mathbf{H}}.
//' This function passes by reference and updates the variables
//' \eqn{\mathtt{WL}} and \eqn{\mathtt{H}} and does not require data to be in an
//' object of type SickleJr. \eqn{\mathtt{RunjrSiCKLSNMF}} calls this function.
//' If your data are in an object of class SickleJr,
//' please use the \eqn{\mathtt{RunjrSiCKLSNMF}} function instead.
//' @name jrSiCKLSNMF
//' @param datamatL An R list where each entry contains a normalized, sparse
//' \eqn{\mathbf{X}^v} matrix corresponding to single-cell modality \eqn{v}
//' @param WL An R list containing initialized values of the \eqn{\mathbf{W}^v}
//' within each modality \eqn{v}
//' @param H A matrix containing initialized values for the shared
//' \eqn{\mathbf{H}}
//' @param AdjL An R list containing all of the adjacency matrices for the
//' feature-feature similarity graphs in sparse format;  note that
//' \eqn{\mathtt{D-Adj}} is the graph Laplacian
//' @param DL An R list containing all of the degree matrices of the
//' feature-feature similarity graphs;  note that \eqn{\mathtt{D-Adj}} is
//' the graph Laplacian
//' @param lambdaWL A list of the \eqn{\lambda_{\mathbf{W}^v}} corresponding to
//' modality \eqn{v}
//' @param lambdaH A double containing the desired value for
//' \eqn{\lambda_{\mathbf{H}}}
//' @param initsamp A vector of randomly selected rows of \eqn{\mathtt{H}} on
//' which to run the objective function
//' @param suppress_warnings A Boolean that indicates whether warnings
//' should be suppressed
//' @param diffFunc A string indicating what type of divergence to use;  set to
//' the Poisson Kullback-Leibler divergence
//' (\eqn{\mathtt{``klp"}}) by default, but the Frobenius norm
//' (\eqn{\mathtt{``fr"}}) is also available
//' @param Hconstraint A string that indicates whether you want to set an L2
//' norm constraint on the rows of \eqn{\mathbf{H}}. Enter `None` for
//' no constraints or `L2Norm` to set the L2 norm of each row of
//' \eqn{\mathbf{H}} to 1, or `Ortho` for an orthogonality constraint on the
//' rows of \eqn{\mathbf{H}}
//' @param Wconstraint A string that indicates whether you want to set an L2
//' norm constraint on the columns of \eqn{\mathbf{W}^v}. Enter `None` for
//' no constraints, `L2Norm` to set the L2 norm of each row of
//' \eqn{\mathbf{W}^v} to 1, or `Ortho` for an orthogonality constraint on the
//' rows of \eqn{\mathbf{W}^v}
//' @param differr A double containing the tolerance
//' @param rounds A double containing the number of rounds
//' @param display_progress A Boolean indicating whether to display the progress
//' bar
//' @param minibatch A Boolean indicating whether to use the mini-batch version
//' of the algorithm
//' @param batchsize Number of batches for mini-batch updates
//' @param random_W_updates A Boolean indicating whether to update
//' \eqn{\mathbf{W}^v} once per epoch (TRUE) or after every update of
//' the subset of \eqn{\mathbf{H}} (FALSE) for the mini-batch algorithm.
//' @param minrounds A minimum number of rounds for the algorithm to run: most
//' useful for the mini-batch algorithm
//' @param numthreads Number of threads to use if running in parallel
//' @returns An R list containing values for the objective function.
//' @export
// [[Rcpp::export]]

Rcpp::List jrSiCKLSNMF(const Rcpp::List& datamatL, Rcpp::List& WL, arma::mat& H,
                        const Rcpp::List& AdjL, const Rcpp::List& DL,
                        const Rcpp::List& lambdaWL, const double& lambdaH,
                        const arma::uvec& initsamp, bool suppress_warnings,
                        const std::vector<std::string>& diffFunc,
                        const std::string& Hconstraint = "None",
                        const std::string& Wconstraint = "None",
                        const double differr = 1e-6,
                        const int rounds = 1000, bool display_progress = true,
                        bool minibatch = true, int batchsize = 100,
                        bool random_W_updates = true, int minrounds = 100,
                        int numthreads = 1) {
  int numv(datamatL.size());
  try {
  ErrorCheck(diffFunc, Hconstraint, Wconstraint);
  } catch(std::invalid_argument& e) {
    Rcpp::Rcerr << e.what() << std::endl;
    return NULL;
  }
  // First convert all Rcpp lists to fields. We will return WL at the end
  arma::field<arma::sp_mat> datamatF(numv, 1);
  arma::field<arma::mat> WF(numv, 1);
  arma::vec lambdaWV(numv);
  arma::field<arma::sp_mat> DF(numv, 1);
  arma::field<arma::sp_mat> AdjF(numv, 1);
  arma::field<arma::sp_mat> datamatFsub(numv, 1);
  int numrowsH = H.n_rows;
  arma::uvec numbercells = arma::linspace<arma::uvec>(0, numrowsH-1, numrowsH);
  arma::field<arma::uvec> numberfeats(numv, 1);
  for (int i = 0;  i < numv; i++) {
    datamatF[i] = Rcpp::as<arma::sp_mat>(datamatL[i]);
    WF[i] = Rcpp::as<arma::mat>(WL[i]);
    int numrowsW = WF[i].n_rows;
    numberfeats[i] = arma::linspace<arma::uvec>(0, numrowsW-1, numrowsW);
    lambdaWV[i] = Rcpp::as<double>(lambdaWL[i]);
    AdjF[i] = Rcpp::as<arma::sp_mat>(AdjL[i]);
    DF[i] = Rcpp::as<arma::sp_mat>(DL[i]);
    if ((sum(diagvec(DF[i])) == 0)&(lambdaWV[i] > 0)) {
      DF[i] = arma::speye(WF[i].n_rows, WF[i].n_rows);
      // Adjacency matrix is already 0s so no need to worry about it
      }
    if (initsamp.n_elem != H.n_rows) {
      arma::sp_mat subsetmati = col_sp(datamatF[i], initsamp);
      datamatFsub[i] = subsetmati;
      }
    }
  arma::field<arma::mat> WFmin(WF);
  arma::mat Hmin(H);
  double minlik;
  std::vector<double> LL;
  double initlikelihood;
  if ((initsamp.n_elem < H.n_rows)) {
    initlikelihood = losscalc(datamatFsub, WF, H.rows(initsamp), AdjF, DF,
                               lambdaWV, lambdaH, diffFunc);
     LL.push_back(initlikelihood);
  } else {
  initlikelihood = losscalc(datamatF, WF, H, AdjF, DF, lambdaWV, lambdaH,
                             diffFunc);
  LL.push_back(initlikelihood);
  }
  Progress p(rounds, display_progress);
  minlik = initlikelihood;
  // Ensure that each column of H is normalized if necessary (since all
  // positive, the derivative of the norm is differentiable everywhere)
  normalizeMat(H, Hconstraint);
  for (int i = 0; i < rounds; i++) {
    p.increment();
    if (minibatch) {
      numbercells = Rcpp::RcppArmadillo::sample(numbercells,
                                                 numbercells.size(), false);
      auto numbersPartition = Utils::partition(numbercells.begin(),
                                                numbercells.end(), batchsize);
      if (Progress::check_abort()) {
        return Rcpp::List::create(Rcpp::Named("Loss") = LL);
       }
#pragma omp parallel for num_threads(numthreads)
       for (unsigned batchCount = 0; batchCount < numbersPartition.size();
       batchCount++) {
         auto batch = numbersPartition[batchCount];
         std::vector<arma::uword> batchsubset;
         std::copy(batch.begin(), batch.end(), std::back_inserter(batchsubset));
         // for (const auto& item:batch) {
            // batchsubset.push_back(item);
            // }
         arma::field<arma::sp_mat> datamatFmini(numv, 1);
         arma::uvec batchsub(batchsubset);
         // Take a subset of datamatF
         for (int j = 0;  j < numv; j++) {
           arma::sp_mat minimatj = col_sp(datamatF[j], batchsub);
           datamatFmini[j] = minimatj;
         }
         arma::mat Hmini = H.rows(batchsub);
         int size = numbersPartition.size();
         perviewNMFMUR(datamatFmini, WF, Hmini, AdjF, DF, lambdaWV, lambdaH,
                       diffFunc, H, Hconstraint, Wconstraint, batchsub,
                       random_W_updates, batchCount, size);
         H.rows(batchsub) = Hmini;
       }} else {
         p.increment();
         arma::uvec batchsub;
         arma::mat emptyH;
         perviewNMFMUR(datamatF, WF, H, AdjF, DF, lambdaWV, lambdaH, diffFunc,
                       emptyH, Hconstraint, Wconstraint, batchsub,
                       random_W_updates);
       }
       if ((initsamp.n_elem < H.n_rows)) {
         double liki = losscalc(datamatFsub, WF, H.rows(initsamp), AdjF, DF,
                                lambdaWV, lambdaH, diffFunc);
         if (liki < minlik) {
           minlik = liki;
           WFmin = WF;
           Hmin = H;
         }
         LL.push_back(liki);
       } else {
          double liki = losscalc(datamatF, WF, H, AdjF, DF, lambdaWV, lambdaH,
                                 diffFunc);
         if (liki < minlik) {
           minlik = liki;
           WFmin = WF;
           Hmin = H;
         }
         LL.push_back(liki);
       }
       if ((((LL.rbegin()[1]-LL.rbegin()[0])) < differr*LL.rbegin()[1])&
           (i > minrounds)) {
         if (LL.rbegin()[1] < LL.rbegin()[0]) {
           if (!suppress_warnings) {
             Rcpp::Rcout << "\n Value for the loss has increased after round "
             << i+1 << ". \n If this is an minibatch algorithm, this behavior "
             "is expected. \n Please note that the final update is: "
             << (abs((LL.rbegin()[1]-LL.rbegin()[0]))/LL.rbegin()[1])*100 <<
             "%. \n Returning W and H matrices corresponding to lowest "
             "achieved loss";
           }
           WF = WFmin;
           H = Hmin;
         }
         break;
       }
       if ((i == rounds-1)&(!suppress_warnings)) {
         Rcpp::Rcout << "Algorithm not converged. Maximum number of rounds "
         "reached. \n Final update is: "
         << (abs((LL.rbegin()[1]-LL.rbegin()[0]))/LL.rbegin()[1])*100 << "%.";
         WF = WFmin;
         H = Hmin;
       }
  }
  for (int i = 0;  i < numv; i++) {
    WL[i] = WF[i];
  }
  return Rcpp::List::create(Rcpp::Named("Loss") = LL);
}
