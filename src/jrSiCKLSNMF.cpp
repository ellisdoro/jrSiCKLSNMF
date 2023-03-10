// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <vector>
#include <math.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

void ErrorCheck(std::string diffFunc, std::string Hconstraint){
  if((diffFunc!="klp")&(diffFunc!="klm")&(diffFunc!="fr")&(diffFunc!="is"))
    throw invalid_argument("Please enter 'klp' for Kullback-Leibler divergence for count data, 'klm' for Kullback-Leibler divergence for proportional data, 'is' for Itakura-Saito divergence (spectra data), or 'fr' for Frobenius norm.");
  if((Hconstraint!="None")&(Hconstraint!="L1Norm")&(Hconstraint!="L2Norm")&(Hconstraint!="L1NormCol"))
    throw invalid_argument("Please enter 'None' for no row constraints on H, 'L1Norm' for an L1 Norm constraint (i.e. all entries in a row sum to 1), 'L1NormCol' for an L1 Norm constraint (i.e. all entries in each column sum to 1), or 'L2Norm' for an L2 Norm constraint (i.e. the square root of the sum of squares of each row entry equals 1)");
}


double lossmatcalc(const arma::sp_mat& datamat, const arma::mat& W,
                         const arma::mat& H, const arma::sp_mat& A,
                         const double lambdaW, const double lambdaH,
                         const std::string diffFunc="klp",int numviews=2){
  double lik=0.0;
  arma::mat Ht=trans(H);
  arma::mat WH=W*Ht;
  if(diffFunc=="klp"){
    arma::mat tmpWH=WH;
    arma::mat tmpWHlog=tmpWH.transform([](double val){
      return std::max(val,1e-16);});
    arma::mat xtmp(datamat);
    xtmp=xtmp.transform([](double val){return std::max(val,1e-16);});
    lik=lik+arma::as_scalar(arma::accu(datamat%log(xtmp)))-arma::as_scalar(arma::accu(datamat%log(tmpWHlog)))-
      arma::as_scalar(arma::accu(datamat))+arma::as_scalar(arma::accu(tmpWH));
  } else if (diffFunc=="klm"){
    arma::mat tmpWH=WH;
    tmpWH=tmpWH.transform([](double val){return std::max(val,1e-16);});
    arma::mat xtmp(datamat);
    xtmp=xtmp.transform([](double val){return std::max(val,1e-16);});
    lik=lik+arma::as_scalar(arma::accu(datamat%log(xtmp)))-arma::as_scalar(arma::accu(datamat%log(tmpWH)));
  } else if (diffFunc=="fr"){
    arma::mat tmp=datamat-WH;
    arma::mat tmpt=trans(tmp);
    double frobnorm2=arma::as_scalar(arma::accu(tmp*tmpt));
    lik=lik+0.5*frobnorm2;
  } else if (diffFunc=="is"){
    arma::mat WHinv=WH;
    arma::mat xtmp(datamat);
    xtmp=xtmp.transform([](double val){return std::max(val,1e-16);});
    WHinv=WHinv.transform([](double val){return std::max(val,1e-16);});
    WHinv=WHinv.transform([](double val){return 1/val;});
    arma::mat WHtmp=WH;
    WHtmp=WHtmp.transform([](double val){return std::max(val,1e-16);});
    lik=lik+arma::as_scalar(arma::accu(datamat%WHinv))-arma::as_scalar(arma::accu(log(xtmp)))+
      arma::as_scalar(arma::accu(log(WHtmp)-1));
  } else{
    perror("Please enter 'klp' for the Kullback-Leibler divergence for count data, 'klm' for the Kullback-Leibler divergence for proportional data, 'is' for the Itakura-Saito divergence, or 'fr' for the Frobenius norm loss");
    return -1;
  }
  if(lambdaW>0){
    arma::mat Wt=trans(W);
    double traceWAW=arma::trace(Wt*A*W);
    lik=lik+0.5*lambdaW*traceWAW;
  }
  if(lambdaH>0){
    double numviewsterm=(double)numviews;
    numviewsterm=1/numviewsterm;
    arma::mat Ht=trans(H);
    double hfrob=arma::as_scalar(arma::accu(arma::trace(H*Ht)));
    lik=lik+0.5*lambdaH*numviewsterm*hfrob;
  }
  return lik;
}

double losscalc(const arma::field<arma::sp_mat>& datamatF, const arma::field<arma::mat>& WF,
                      const arma::mat& H,const arma::field<arma::sp_mat>& AF,
                      arma::vec lambdaWV, const double& lambdaH,
                      const std::string& diffFunc="klp"){
  double lik=0.0;
  int viewnum(datamatF.n_rows);
  for(int i=0;i<viewnum;i++){
    double liki=lossmatcalc(datamatF[i],WF[i],H,AF[i],lambdaWV[i],lambdaH,diffFunc,viewnum);
    if (liki==-1){
      break;
    }else{
      lik=lik+liki;
    }
  }
  return lik;
}
void NMFinview(arma::mat& old,arma::mat denom,arma::mat numer){
  numer=numer.transform([](double val){return std::max(val,1e-16);});
  denom=denom.transform([](double val){return std::max(val,1e-16);});
  denom=denom.transform([](double val){return 1/val;});
  old=old%numer%denom;
}
//Perform tangent gradient as in Le Roux et al. 2015
void normalizeH(arma::mat& H, const std::string Hconstraint){
  if(Hconstraint=="L1Norm"){
    arma::rowvec Hcolsums=sum(H);
    H.each_row()/=Hcolsums;
    return;
  }else if (Hconstraint=="L2Norm"){
    H=normalise(H);
    return;
  }else if (Hconstraint=="L1NormCol"){
    arma::colvec Hrowsums=sum(H,1);
    H.each_col()/=Hrowsums;
  }else{
    return;
  }
}
arma::mat regFunc(const arma::mat& denomnumer, const arma::mat& H, const std::string Hconstraint){
  arma::mat regularizernumerdenom(denomnumer.n_rows,denomnumer.n_cols,arma::fill::value(0.0));
  if(Hconstraint=="L1Norm"){
    for(int i=0;i<(int)H.n_cols;i++){
      arma::mat Hitrans=H.col(i).t();
      arma::mat Hones(H.n_rows,1,arma::fill::value(1.0));
      regularizernumerdenom.col(i)=Hones*Hitrans*denomnumer.col(i);
    }
    return regularizernumerdenom;
  } else if (Hconstraint=="L2Norm"){
    arma::mat squareones(H.n_rows,H.n_rows);
    squareones.fill(1);
    regularizernumerdenom=H%(squareones*(denomnumer));
    return regularizernumerdenom;
  }else if(Hconstraint=="L1NormCol"){
    for(int i=0;i<(int)H.n_rows;i++){
      arma::mat Hitrans=H.row(i).t();
      arma::mat Hones(H.n_cols,1,arma::fill::value(1.0));
      regularizernumerdenom.row(i)=Hones*(Hitrans*denomnumer.row(i));
    }
    return regularizernumerdenom;
  }else{
    return regularizernumerdenom;
  }
}

void perviewNMFMUR(const arma::field<arma::sp_mat>& datamatF, arma::field<arma::mat>& WF,
                   arma::mat& H,const arma::field<arma::sp_mat>& AF,
                   const arma::vec lambdaWV, const double lambdaH,
                   const std::string diffFunc, const std::string Hconstraint){
  int viewnum(datamatF.n_rows);
  for(int i=0;i<viewnum;i++){
    arma::sp_mat D=AF[i];
    arma::sp_mat Adj=AF[i];
    D=D.transform([](double val){return std::max(val,0.0);});
    Adj=Adj.transform([](double val){return std::min(val,0.0);});
    Adj=-Adj;
    if((sum(diagvec(D))==0)&(lambdaWV[i]>0)){
      D=arma::speye(WF[i].n_rows,WF[i].n_rows);
      //Adjacency matrix is already 0s so no need to worry about it
    }
    arma::mat numerW(WF[i].n_rows,WF[i].n_cols,arma::fill::value(0.0));
    arma::mat denomW(WF[i].n_rows,WF[i].n_cols,arma::fill::value(0.0));
    if(diffFunc=="fr"){
      arma::mat HH=H.t()*H;
      numerW=numerW+WF[i]*HH;
      denomW=denomW+datamatF(i,0)*H;
    }else if(diffFunc=="klp"){
      arma::mat WHinv=WF[i]*H.t();
      WHinv=WHinv.transform([](double val){return std::max(val,1e-16);});
      WHinv=WHinv.transform([](double val){return 1/val;});
      arma::mat XWHinv(datamatF[i]%WHinv);
      numerW=numerW+XWHinv*H;
      denomW=denomW+repmat(sum(H),WF[i].n_rows,1);
    }else if (diffFunc=="klm"){
      arma::mat WHinv=WF[i]*H.t();
      WHinv=WHinv.transform([](double val){return std::max(val,1e-16);});
      WHinv=WHinv.transform([](double val){return 1/val;});
      arma::sp_mat XWHinv(datamatF[i]%WHinv);
      numerW=numerW+XWHinv*H;
    }else if (diffFunc=="is"){
      arma::mat WHinv=WF[i]*H.t();
      WHinv=WHinv.transform([](double val){return std::max(val,1e-16);});
      WHinv=WHinv.transform([](double val){return 1/val;});
      arma::sp_mat XWHinv(datamatF[i]%WHinv);
      arma::mat WHinvsquare=WHinv.transform([](double val){return pow(val,2);});
      WHinvsquare=WHinvsquare.transform([](double val){return std::max(val,1e-16);});
      arma::sp_mat XWHinvsquare(datamatF[i]%WHinvsquare);
      numerW=numerW+XWHinvsquare*H;
      denomW=denomW+XWHinv*H;
    }
    if(lambdaWV[i]>0){
      denomW=denomW+0.5*lambdaWV[i]*(D*WF[i]+D.t()*WF[i]);
      numerW=numerW+0.5*lambdaWV[i]*(Adj*WF[i]+Adj.t()*WF[i]);
    }
    NMFinview(WF[i],denomW,numerW);
    WF[i]=WF[i].transform([](double val){return(val<1e-10)? double(0) : double(val);});
  }
  arma::mat numerH(H.n_rows,H.n_cols,arma::fill::value(0.0));
  arma::mat denomH(H.n_rows,H.n_cols,arma::fill::value(0.0));

  arma::mat InumerH(H.n_rows,H.n_cols,arma::fill::value(0.0));
  arma::mat IdenomH(H.n_rows,H.n_cols,arma::fill::value(0.0));

  arma::mat regnumerH(H.n_rows,H.n_cols,arma::fill::value(0.0));
  arma::mat regdenomH(H.n_rows,H.n_cols,arma::fill::value(0.0));
  //Adding in column regularization to improve identifiability
  //normalize H to have l2 norm equal to 1
  //new additions are all multiplied by HHones for reference later
  //Since H is non-negative, |H|=H and is differentiable for all H
  //arma::mat Hones= H.ones();

  //Now we need to divide each column by the squared L2 norm of the corresponding
  //column entries of HHones.Since each entry is 1, this is the same as
  //dividing by the the number of rows of HHones.
  //See Douglas, 2000; Fu 2018.
  for(int i=0;i<viewnum;i++){
    if(diffFunc=="fr"){
      arma::mat WW=WF[i].t()*WF[i];

      //Take gradient with normalized parameters
      InumerH=H*WW;
      IdenomH=datamatF[i]*H;
      if(Hconstraint!="None"){
        regnumerH=regFunc(IdenomH,H,Hconstraint);
        regdenomH=regFunc(InumerH,H,Hconstraint);
      }

    } else if(diffFunc=="klp"){
      arma::mat WHinv=WF[i]*H.t();
      WHinv=WHinv.transform([](double val){return std::max(val,1e-16);});
      WHinv=WHinv.transform([](double val){return 1/val;});
      arma::sp_mat XWHinv(datamatF[i]%WHinv);

      InumerH=XWHinv.t()*WF[i];
      IdenomH=repmat(sum(WF[i]),H.n_rows,1);
      if(Hconstraint!="None"){
        regnumerH=regFunc(IdenomH,H,Hconstraint);
        regdenomH=regFunc(InumerH,H,Hconstraint);
      }
    } else if(diffFunc=="klm"){

      arma::mat WHinv=WF[i]*H.t();
      WHinv=WHinv.transform([](double val){return std::max(val,1e-16);});
      WHinv=WHinv.transform([](double val){return 1/val;});
      arma::sp_mat XWHinv(datamatF[i]%WHinv);
      InumerH=XWHinv.t()*WF[i];
      if(Hconstraint!="None"){
        //Here regnumerH and IdenomH are not updated since the gradient is all negative
        regdenomH=regFunc(InumerH,H,Hconstraint);
      }
    } else if(diffFunc=="is"){
      arma::mat WHinv=WF[i]*H.t();
      WHinv=WHinv.transform([](double val){return std::max(val,1e-16);});
      WHinv=WHinv.transform([](double val){return 1/val;});
      arma::sp_mat XWHinv(datamatF[i]%WHinv);
      arma::mat WHinvsquare=WHinv.transform([](double val){return pow(val,2);});
      WHinvsquare=WHinvsquare.transform([](double val){return std::max(val,1e-16);});
      arma::sp_mat XWHinvsquare(datamatF[i]%WHinvsquare);
      InumerH=XWHinvsquare.t()*WF[i];
      IdenomH=XWHinv.t()*WF[i];
      if(Hconstraint!="None"){
        regnumerH=regFunc(IdenomH,H,Hconstraint);
        regdenomH=regFunc(InumerH,H,Hconstraint);
      }
    }
    numerH=numerH+InumerH+regnumerH;
    denomH=denomH+IdenomH+regdenomH;
  }
  if(lambdaH>0){
    IdenomH=H;
    if(Hconstraint!="None"){
      regnumerH=regFunc(IdenomH,H,Hconstraint);
    }

    numerH=numerH+lambdaH*regnumerH;
    denomH=denomH+lambdaH*IdenomH;
  }
  NMFinview(H,denomH,numerH);
  H=H.transform([](double val){return(val<1e-10)? double(0) : double(val);});

  normalizeH(H,Hconstraint);
}
//' @title jrSiCKLSNMF
//' @description Perform joint non-negative matrix factorization (NMF) across multiple views of single cell data.
//' Users can choose to use the Poisson Kullback-Leibler divergence or the Frobenius norm.
//' Users can also set graph regularization constraints on W and sparsity constraints on H.
//' This function updates WL and H.
//' @name jrSiCKLSNMF
//' @param datamatL An R list where each entry contains a sparse X matrix corresponding to a single cell data view
//' Each X is m^v features by n cells. Features can differ across matrices; however, n must be the same for
//' each view. All data are measured on the same set of cells
//' @param WL An R list containing initialized values for the W within each view. These are passed by reference
//' @param H A matrix containing initialized values for the shared H
//' @param AL An R list containing all of the graph laplacians in sparse format
//' @param lambdaWL A list of each lambdaW for each view
//' @param lambdaH A double containing the desired value for H
//' @param diffFunc A string indicating what type of divergence to use. It is Poisson Kulback Leibler by default
//' @param Hconstraint A string that indicates whether you want to force constraints on the rows of H. Enter "None" for
//' no constraints, enter "L1Norm" to ensure all rows of H sum to 1, and "L2Norm" to ensure that the the
//' L2 norm of each row of H equals 1. Please note that jrSiCKLSNMF stores H as the transpose of H, so the code
//' will perform regularization on the columns of the transpose of H.
//' @param differr A double containing the tolerance
//' @param rounds A double containing the number of rounds
//' @param display_progress A boolean to indicate whether the user wants to display the progress
//' @returns An R list containing values for the objective function.
//' @export
// [[Rcpp::export]]
 Rcpp::List jrSiCKLSNMF(const Rcpp::List& datamatL, Rcpp::List& WL, arma::mat& H,
                        const Rcpp::List& AL,const Rcpp::List& lambdaWL,
                        const double& lambdaH=0.0,const std::string diffFunc="klp",
                        const std::string Hconstraint="L2Norm",const double differr=1e-6,
                        const double rounds=10000, bool display_progress=true){
   int viewnum(datamatL.size());
   try{
     ErrorCheck(diffFunc,Hconstraint);
   }catch(invalid_argument& e){
     Rcerr <<e.what()<<endl;
     return NULL;
   }
   //First convert all Rcpp lists to fields. We will return WL at the end
   arma::field<arma::sp_mat> datamatF(viewnum,1);
   arma::field<arma::mat> WF(viewnum,1);
   arma::vec lambdaWV(viewnum);
   arma::field<arma::sp_mat> AF(viewnum,1);
   for(int i=0; i<viewnum;i++){
     datamatF[i]=as<arma::sp_mat>(datamatL[i]);
     WF[i]=as<arma::mat>(WL[i]);
     AF[i]=as<arma::sp_mat>(AL[i]);
     lambdaWV[i]=as<double>(lambdaWL[i]);
   }
   std::vector<double> LL;
   double initlikelihood;
   initlikelihood=losscalc(datamatF,WF,H,AF,lambdaWV,lambdaH,diffFunc);
   LL.push_back(initlikelihood);
   Progress p(rounds,display_progress);
   //Ensure that each column of H is L2 normalized (since all positive, the derivative
   //of the norm is differentiable everywhere)
   H=normalise(H);
   for(int i=0;i<rounds;i++){
     if(Progress::check_abort()){
       return Rcpp::List::create(Rcpp::Named("W")=WL,
                                 Rcpp::Named("H")=H,
                                 Rcpp::Named("Loss")=LL);
     }
     p.increment();
     perviewNMFMUR(datamatF,WF,H,AF,lambdaWV,lambdaH,diffFunc,Hconstraint);
     double liki=losscalc(datamatF,WF,H,AF,lambdaWV,lambdaH,diffFunc);
     LL.push_back(liki);
     if(abs((LL.end()[-2]-LL.end()[-1]))<differr*LL.end()[-2]){
       break;
     }
     if(i==rounds-1){
       Rcout<<"Algorithm not converged. Maximum number of rounds reached. Final update is: "<<(abs((LL.end()[-2]-LL.end()[-1]))/LL.end()[-2])*100<<"%.";
     }
   }
   for(int i=0; i<viewnum;i++){
     WL[i]=WF[i];
   }
   return Rcpp::List::create(Rcpp::Named("Loss")=LL);

 }
