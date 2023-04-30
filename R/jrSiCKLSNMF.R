NULL

#' @title A SickleJr object
#' @description Defines the SickleJr class for use with jrSiCKLSNMF
#' @slot count.matrices list. List of all of the QC'd count matrices.
#' Note that these count matrices should not use all features and should
#' only include features that appear in at a minimum 10 cells.
#' @slot normalization.type character. Holds the type of normalization
#' @slot normalized.count.matrices list. Holds the normalized count matrices
#' @slot graph.laplacian.list list. A list of the graph laplacians used for graph regularizations
#' @slot rowRegularization character. A string that indicates the type of row regularization to
#' use. Types include "None" and "L2Norm"
#' @slot diffFunc character. Holds the name of the function used to measure the "distance" between
#' data matrix X and WH for each view. Can be "klp" for Kullback-Leibler based on the Poisson
#' distribution or "fr" for the Frobenius norm.
#' @slot lambdaWlist list. List of lambda values used in the W view
#' @slot lambdaH numeric. Value for the constraint on H
#' @slot Wlist list. A list of the generated W matrices, one for each view
#' @slot H matrix. The shared H matrix
#' @slot WHinitials If in PlotLossvsLatentFactors, you use all of the cells to calculate,
#' you can store these in a WHinitials vector and use it when performing jrSiCKLSNMF to save
#' time
#' @slot lossCalcSubsample holds the cell indices on which plotlossvslatentfactors was calculated
#' @slot latent.factor.elbow.values Values for the latent factor elbow plot
#' @slot online Indicator variable that states whether the algorithm is online
#' @slot clusterdiagnostics list. List of the cluster diagnostic results for the SickleJr object
#' @slot clusters list. List of different clusters performed on the SickleJr object
#' @slot metadata list. List of metadata
#' @slot loss vector. Vector of the value for the loss function
#' @slot umap list. List of different UMAP based dimension reductions
#' @return an object of type SickleJr
#' @export
SickleJr<-methods::setClass(
  "SickleJr",
  slots=c(
    count.matrices="list",
    normalization.type="character",
    normalized.count.matrices="list",
    graph.laplacian.list="list",
    rowRegularization="character",
    diffFunc="character",
    lambdaWlist="list",
    lambdaH="numeric",
    Wlist="list",
    H="matrix",
    WHinitials="list",
    lossCalcSubsample="vector",
    latent.factor.elbow.values="data.frame",
    online="logical",
    clusterdiagnostics="list",
    clusters="list",
    metadata="list",
    loss="vector",
    umap="list"
  )
)

#' @title Create an object of type SickleJr
#' @name CreateSickleJr
#' @description Creates an object of type SickleJr that includes lists of sparse
#' count matrices and allows users to specify the names of those count matrices
#' @param count.matrices A list of count matrices; each view is one entry in the
#' list. These count matrices should already be QC'd and the features filtered
#' @param names Optional parameter with names for the count matrices
#' @return A SickleJr object with sparse count matrices
#' @examples
#' SickleJrSim<-CreateSickleJr(SimData$Xmatrices)
CreateSickleJr<-function(count.matrices,names=NULL){
  if(length(names)!=length(count.matrices)){
    warning("\n Length of names does not match length of list of count matrices. Names of count matrices will remain unchanged.\n")
    names<-NULL
  }
  if(!is.null(names)){
    names(count.matrices)<-names
  }
  #Force matrix to be stored in sparse format.
  count.matrices<-lapply(count.matrices, function(x){if(!.is.sparseMatrix(x)){
    x<-.sparsify(x)
  }else{x}})
  object<-methods::new(Class="SickleJr",count.matrices=count.matrices)
  return(object)
}

#' @title Build KNN graphs and generate their graph Laplacians
#' @description Generate graph Laplacians for graph regularization of
#' jrSiCKLSNMFNMF from the list of normalized count matrices
#'@name BuildKNNGraphLaplacians
#'@param SickleJr A SickleJr object
#'@param k Number of KNN neighbors to calculate. By default, is set to 20
#'@returns A list of graph Laplacians in sparse matrix format
#'@export
BuildKNNGraphLaplacians<-function(SickleJr,k=20){
    counts<-SickleJr@count.matrices
    SickleJr@graph.laplacian.list<-lapply(counts,function(x){
    laplacian_matrix(buildKNNGraph(x,transposed=TRUE,k=k))})
  return(SickleJr)
}

#' @title Build SNN graphs and generate their graph Laplacians
#' @description Generate graph Laplacians for graph regularization of
#' jrSiCKLSNMFNMF from the list of normalized count matrices
#'@name BuildSNNGraphLaplacians
#'@param SickleJr A SickleJr object
#'@param k Number of KNN neighbors to calculate. By default, is set to 20
#'@returns A list of graph Laplacians in sparse matrix format
#'@export
BuildSNNGraphLaplacians<-function(SickleJr,k=20){
    counts<-SickleJr@count.matrices
    SickleJr@graph.laplacian.list<-lapply(counts,function(x){
    laplacian_matrix(buildSNNGraph(x,transposed=TRUE,k=k))})
  return(SickleJr)
}



#' @title Normalize the count matrices and set whether to use the KL divergence
#' based on the Poisson distribution or the Frobenius norm
#' @description Perform normalization for count data. If you are planning to use the Frobenius norm,
#'set frob=TRUE to log(x+1) normalize your count data. This step can be skipped for percentage data and
#'spectra data. You may also skip this if you would like to perform a different form
#'of normalization; however, please ensure that if using the Frobenius norm that all of the values are
#'divided by the sum of all of the values.
#'@name NormalizeCountMatrices
#'@param SickleJr An object of type SickleJr with a count matrix. Users can choose to normalize using median
#'library size normalization or log(x+1) normalization
#'@param diffFunc A string. Set to "klp" to use the KL divergence based on the Poisson distribution
#'or to "fr" to use the Frobenius norm. Defaults to KL divergence. This also determines
#'the type of normalization to use.
#'@param scaleFactor A list of numeric values to use as scale factors for the log(x+1)
#'normalization when using the Frobenius diffFunction
#'@returns A list of sparse, normalized matrices
#'@export

NormalizeCountMatrices<-function(SickleJr,diffFunc="klp",scaleFactor=NULL){
  Xmatrixlist<-SickleJr@count.matrices
  SickleJr@diffFunc<-diffFunc
  frob=FALSE
  if(diffFunc=="fr"){
    frob=TRUE
  }
  NormalizedXmatrices<-Xmatrixlist
  medianlibsize<-list()
  scaleFactorlist<-list()
  if(!frob){
    medianlibsize<-lapply(Xmatrixlist,function(x) {median(unlist(
      lapply(.listCols(x),function(y) sum(y))))})
    for(i in 1:length(Xmatrixlist)){
      normalizeXmatrixlist<-.listCols(Xmatrixlist[[i]])
      normalizeXmatrixlist<-lapply(normalizeXmatrixlist,function(x){x/sum(x)*medianlibsize[[i]]})
      NormalizedXmatrices[[i]]@x<-unlist(normalizeXmatrixlist)
    }
  } else{
    if(is.list(scaleFactor)){
      if(length(scaleFactor)==length(Xmatrixlist)){
        scaleFactorlist=scaleFactor
      } else{
        if(length(scaleFactor)==1){
          warning(paste0("\n List contains only 1 scale factor. Setting all scale factors to ",scaleFactor[[1]],".\n"))
          for(i in 1:length(Xmatrixlist)){
            scaleFactorlist[[i]]=scaleFactor[[1]]
          }
        } else{
          stop("\n Length of scale factor list not equal to number of assays\n")
        }
      }
    } else if(is.null(scaleFactor)){
      warning("\n No scale factor given. Setting to 1e6 \n")
      for(i in 1:length(Xmatrixlist)){
        scaleFactorlist[[i]]=1e6
      }
    } else if (is.numeric(scaleFactor)){
      warning(paste0("\n Only one scale factor specified. Setting all to ",scaleFactor,".\n"))
      for(i in 1:length(Xmatrixlist)){
        scaleFactorlist[[i]]=scaleFactor
      }
    } else{
      stop("\n Scale factor is not numeric or a list of proper length.\n")
    }
    for(i in 1:length(Xmatrixlist)){
      normalizeXmatrixlist<-.listCols(Xmatrixlist[[i]])
      scaleFactori=scaleFactorlist[[i]]
      normalizeXmatrixlist<-lapply(normalizeXmatrixlist,function(x){x/sum(x)*scaleFactori})
      normalizeXmatrixlist<-lapply(normalizeXmatrixlist,function(x){log(x+1)})
      NormalizedXmatrices[[i]]@x<-unlist(normalizeXmatrixlist)
    }
  }
  SickleJr@normalized.count.matrices<-NormalizedXmatrices
  return(SickleJr)
}


#' @title Set lambda values, type of row regularization, for a SickleJr object
#' @description Provide the values for the graph regularization \eqn{\lambda_{\textbf{W}^v}}
#' for each view as a list and provide a
#' @name SetLambdasandRowReg
#' @param SickleJr a SickleJr object
#' @param lambdaWlist A list of graph regularization constraints for the W matrices.
#' Defaults to 2 views with the RNA view equal to 10 and the ATAC view equal to 50
#' @param lambdaH lambda H value to add. Defaults to 500.
#' @param rowReg Choose whether or not to enforce constraints on the rows. Enter
#' "None" for no constraints and "L2Norm" for an L2 Norm constraint on the rows
#' of H. Defaults to "None"
#' @return SickleJr object with added lambda values
#' @export
#'
SetLambdasandRowReg<-function(SickleJr,lambdaWlist=list(10,50),lambdaH=500,rowReg="None"){
  SickleJr@lambdaWlist<-lambdaWlist
  SickleJr@lambdaH<-lambdaH
  SickleJr@rowRegularization<-rowReg
  return(SickleJr)
}

#' @title Generate the W and H matrices
#' @description Create the W and H matrices via non-negative double singular
#' value decomposition (nndsvd) or randomization. For randomization, we will run
#' the algorithm for 10 steps at the chosen value and pick the W and H values with
#' the lowest achieved loss
#' @name GenerateWmatricesandHmatrix
#' @param SickleJr a SickleJr object
#' @param d number of latent factors to use. Defaults to 10.
#' @param random Choose whether to use nndsvd or random initialization. Defaults to nndsvd
#' @param numberReps Number of random initializations to use. Defaults to 5.
#' @param seed Random seed for reproducibility with random initializations
#' @param online Indicates whether or not to use an online algorithm
#' @param batchsize Indicates size of batches for online NMF
#' @param random_W_updates Indicates whether to only update W once per round of
#' H updates. Only appropriate for online algorithms and for randomly generated
#' algorithms
#' @param subsample A vector of values to use for subsampling. Only appropriate
#' for determining proper values for d.
#' @returns SickleJr A SickleJr object with W and H matrices added.
#' @export
#'
GenerateWmatricesandHmatrix<-function(SickleJr,d=10,random=FALSE,
                                      numberReps=100,seed=5,online=FALSE,batchsize=-1,
                                      random_W_updates=FALSE,subsample=1:dim(SickleJr@count.matrices[[1]])[2]){
  Hconcatmat<-NULL
  NormalizedXmatrices<-SickleJr@normalized.count.matrices
  if(length(subsample)<dim(SickleJr@normalized.count.matrices[[1]])[2]){
    for(i in 1:length(SickleJr@normalized.count.matrices)){
      NormalizedXmatrices[[i]]<-SickleJr@normalized.count.matrices[[i]][,subsample]
    }
  }
  rowReg=SickleJr@rowRegularization

  if(!random){
    svdlist<-lapply(NormalizedXmatrices,function(x) .nndsvd(A=x,k=d,flag=1))
    WandMeanlist<-lapply(svdlist,function(x) list(W=x$W,meanW=mean(apply(x$W,2,function(y) sum(y)))))
    SickleJr@Wlist<-lapply(WandMeanlist,function(x) apply(x$W,2,function(y) y/sum(y))*x$meanW)
    bigXmatrix<-do.call(rbind,NormalizedXmatrices)
    svdH<-.nndsvd(A=bigXmatrix,k=d,flag=1)
    Hconcatmat<-t(svdH$H)
    if(rowReg=="L2Norm"){
      norms<-apply(Hconcatmat,2,function(x)sqrt(sum(x^2)))
      for(i in 1:dim(Hconcatmat)[[2]]){
        Hconcatmat[,i]=Hconcatmat[,i]/norms[i]
      }
    } else if(rowReg=="L1Norm"){
      Hconcatmat<-apply(Hconcatmat,2,function(x)x/sum(x))
    } else{
      meanHconcat<-mean(apply(Hconcatmat,2,function(x) sum(x)))
      Hconcatmat<-apply(Hconcatmat,2,function(x) x/sum(x))*meanHconcat
    }
    SickleJr@H<-Hconcatmat
  }else{
    output<-rep(0,numberReps)
    for(i in 1:numberReps){

    set.seed(seed+i)
    Unormvals<-lapply(NormalizedXmatrices,
                      function(x) {list(normvals=median(unlist(.listCols(x))),
                      dimW=dim(x)[1])})
    Hvals<-mean(unlist(lapply(Unormvals,function(x) x$normvals)))
    Wvals<-lapply(Unormvals,function(x) {x$normvals<-x$normvals/(Hvals*sqrt(d))
                                          return(x)})
    Hvals<-Hvals/sqrt(d)
    Hdim<-dim(NormalizedXmatrices[[1]])[2]
    Wnew<-lapply(Wvals,function(x) matrix(runif(d*x$dimW,min=0,max=x$normvals),nrow = x$dimW,ncol=d))
    Hconcatmat<-matrix(runif(Hdim*d,min=0,max=Hvals),nrow=Hdim,ncol=d)
    if(rowReg=="L2Norm"){
      norms<-apply(Hconcatmat,2,function(x)sqrt(sum(x^2)))
      for(i in 1:dim(Hconcatmat)[[2]]){
        Hconcatmat[,i]=Hconcatmat[,i]/norms[i]
      }
    } else if(rowReg=="L1Norm"){
      Hconcatmat<-apply(Hconcatmat,2,function(x)x/sum(x))
    } else{
      meanHconcat<-mean(apply(Hconcatmat,2,function(x) sum(x)))
      Hconcatmat<-apply(Hconcatmat,2,function(x) x/sum(x))*meanHconcat
    }
    if(i==1){
      old.W=copy(Wnew)
      old.H=copy(Hconcatmat)
    }
    SickleJr@H<-Hconcatmat
    SickleJr@Wlist<-Wnew
    rm(Hconcatmat)
    rm(Wnew)
    SickleJr<-RunjrSiCKLSNMF(SickleJr,rounds=5,online=online,random_W_updates=random_W_updates,
                   batchsize=batchsize,seed=seed,display_progress = FALSE,suppress_warnings = TRUE)
    output[i]=tail(SickleJr@loss$Loss,1)
    if(which.min(output[which(output>0)])!=i){
      SickleJr@H<-copy(old.H)
      SickleJr@Wlist=copy(old.W)
    }else{
      old.H=copy(SickleJr@H)
      old.W=copy(SickleJr@Wlist)
    }
    }
  }
  return(SickleJr)
}

#' @title Creating plots to help determine the number of latent factors to use for jrSiCKLSNMF
#' @description This generates plots of the lowest achieved loss after a
#' user-specified number of iterations (default 200)
#' of the jrSiCKLSNMF algorithm for each latent factor (defaults to 2:20). This operates similarly to a
#' scree plot, so please select a number of latent factors that corresponds to the
#' elbow of the plot. This function will not uniformly decrease, so an increase in loss
#' may be expected. Select the minimum loss or a point that appears to be an elbow.
#' @name PlotLossvsLatentFactors
#' @param SickleJr A SickleJr object
#' @param rounds Number of rounds to use. Defaults to 20. This process will take some time,
#' so a high number of rounds is not recommended
#' @param differr A tolerance of updates after which to stop updating. For these plots,
#' this is set to 1e-4 by default
#' @param d_vector Vector of D values to test over. Defaults from 2 to 20.
#' @param parallel Indicates whether or not you want to perform this using
#' parallel computing
#' @param nCores Number of desired cores. If null, defaults to the number of cores
#' minus 1 for user convenience.
#' @param subsampsize Indicates whether you want to perform this on a random subsample rather
#' than on the whole dataset. Will speed up the process but will have lower accuracy on desired
#' number of d
#' @param online Indicates whether you want an online algorithm. Defaults to FALSE
#' @param random Indicates whether to use random initialization to generate the W and H matrices
#' Defaults to FALSE as this in general is less accurate.
#' @param random_W_updates Only used for online algorithms. If TRUE, only updates W using the second
#' to last H subset on the online algorithm
#' @param seed Set a random seed
#' @param batchsize Desired batch size. Do not use if you are using a subsample.
#' @param lossonsubset Indicates whether to calculate the loss on a subset rather than the full dataset.
#' Speeds up computation for larger datasets.
#' @param losssubsetsize Number of cells to use for the loss subset. Defaults to full number of cells
#' @return A SickleJr with an added ggplot object in slot LatentFactorDeterminationPlot
#' @export

PlotLossvsLatentFactors<-function(SickleJr,rounds=200,differr=1e-5,d_vector=c(2:20),
                                  parallel=FALSE,nCores=detectCores()-1,subsampsize=NULL,
                                  online=FALSE,random=FALSE,random_W_updates=FALSE,seed=NULL,batchsize=-1,
                                  lossonsubset=FALSE,losssubsetsize=dim(SickleJr@count.matrices[[1]])[2]){
  latentfactors<-SickleJr@latent.factor.elbow.values$latent_factor_number
  if(any(latentfactors%in%d_vector)){
    alreadycalced<-dput(latentfactors[which(latentfactors%in%d_vector)])
      print(paste(alreadycalced,
            "has already been calculated. Skipping calculation for this value"))
    d_vector<-d_vector[-which(d_vector%in%latentfactors)]
    emptyvec<-is_empty(d_vector)
    if(emptyvec){
      latent_factor_number<-SickleJr@latent.factor.elbow.values$latent_factor_number
      Loss<-SickleJr@latent.factor.elbow.values$Loss
      elbowvals<-data.frame(latent_factor_number=latent_factor_number,Loss=Loss)
      s=""
      if(rounds>1){
        s="s"
      }
      plotvals<-ggplot(elbowvals,aes(x=latent_factor_number,y=Loss))+geom_line()+geom_point()+theme_bw()+ggtitle(paste0("Plot of Number of Latent Factors vs Loss after ",rounds," Iteration",s))+xlab("Number of Latent Factors")
      print(plotvals)
      stop("No uncalculated values for d. Plotting loss values for latent factors and exiting...\n")
    }
  }
  if(is.numeric(seed)){
    set.seed(seed)
  }
  initsamp=NULL
  if(!lossonsubset){
    initsamp=1:losssubsetsize
  }else{
    initsamp<-sample(1:dim(SickleJr@normalized.count.matrices[[1]])[2],losssubsetsize,replace=FALSE)
  }
  if(is.numeric(subsampsize)&(online)){
    stop("It is not appropriate to use both an online algorithm and a subsample.
         Please select only one.\n")
  }
  samp<-NULL
  if(is.null(subsampsize)){
    subsampsize=dim(SickleJr@count.matrices[[1]])[2]
    message("\n Subsample not specified. Setting to number of cells. \n")
    samp<-1:subsampsize
  } else  if(subsampsize>=dim(SickleJr@count.matrices[[1]])[2]){
    warning("\n Subsample size is greater than or equal to number of cells. Using all cells.\n")
    subsampsize=dim(SickleJr@count.matrices[[1]])[2]
    samp=1:dim(SickleJr@count.matrices[[1]])[2]
  } else if (!is.numeric(subsampsize)){
    stop("\n Subsample size is not numeric. Exiting...")
  }else {
    samp=sample(1:dim(SickleJr@count.matrices[[1]])[2],subsampsize,replace=FALSE)
  }
  if(random_W_updates & !online){
    warning("\n Random W updates are only appropriate for online algorithms. Setting random_W_updates to FALSE.\n" )
    random_W_updates=FALSE
  }
  if(batchsize>dim(SickleJr@count.matrices[[1]])[2]){
    warning("\n Batch size larger than number of cells. Will not use an online algorithm \n")
    batchsize=-1
    online=FALSE
  }
  AdjL=lapply(SickleJr@graph.laplacian.list,function(x) {x@x[which(x@x>0)]=0
  return(-x)})
  DL=lapply(SickleJr@graph.laplacian.list,function(x) {x@x[which(x@x<0)]=0
  return(x)})
  SickleJrSub<-list()
  for(i in 1:length(SickleJr@normalized.count.matrices)){
    SickleJrSub[[i]]<-SickleJr@normalized.count.matrices[[i]][,samp]
  }
  if(!online){
    initsamp=1:dim(SickleJrSub[[1]])[2]
  }
  if(losssubsetsize<dim(SickleJrSub[[1]])[2]){
    initsamp<-sample(1:dim(SickleJrSub[[1]])[2],batchsize,replace=FALSE)
    #prepare initsamp for c++
    initsamp<-initsamp-1
  }else{initsamp=1:dim(SickleJrSub[[1]])[2]-1}
  WHinitials<-NULL
  lossvals<-NULL
  if(!parallel){
    cl=NULL} else{
      cl<<-makeCluster(nCores)
      clusterExport(cl,varlist=c("SickleJr","d_vector","random","samp","AdjL","DL","differr","rounds","initsamp","SickleJrSub"),envir = environment())
      clusterEvalQ(cl,library("jrSiCKLSNMF"))
    }
  message("\n Calculating initial WH matrices... \n")
  WHinitials<-pblapply(X=d_vector,FUN=function(x) {newvals=GenerateWmatricesandHmatrix(SickleJr,d=x,random=random,subsample=samp)
  list(d=x,Wlist=newvals@Wlist,H=newvals@H)},cl=cl)
  message("\n Running jrSiCKLSNMF... \n")
  lossvals<-pblapply(X=WHinitials,FUN=function(x){
      vals<-jrSiCKLSNMF(datamatL=SickleJrSub,
                      WL=x[["Wlist"]],
                      H=x[["H"]],
                      AdjL=AdjL,
                      DL=DL,
                      lambdaWL=SickleJr@lambdaWlist,
                      lambdaH=SickleJr@lambdaH,
                      Hconstraint = toString(SickleJr@rowRegularization),
                      differr=differr,
                      display_progress=FALSE,
                      rounds=rounds,
                      diffFunc = toString(SickleJr@diffFunc),
                      minrounds=200,
                      suppress_warnings=TRUE,
                      initsamp = 1:dim(SickleJrSub[[1]])[2])
    return(vals[length(vals)])},cl=cl)
   if(parallel) {
     stopCluster(cl)
   }
    newlikelihoodvalues<-NULL
    for(i in 1:length(lossvals)){
      newlikelihoodvalues<-rbind(newlikelihoodvalues,min(lossvals[[i]]$Loss))
    }
    oldvals<-SickleJr@latent.factor.elbow.values
    SickleJr@latent.factor.elbow.values<-rbind(oldvals,data.frame(latent_factor_number=d_vector,Loss=newlikelihoodvalues))
    elbowvals=SickleJr@latent.factor.elbow.values
    plotvals<-ggplot(elbowvals,aes(x=latent_factor_number,y=Loss))+geom_line()+geom_point()+theme_bw()+ggtitle(paste0("Plot of Number of Latent Factors vs Loss after ",rounds," Iterations"))+xlab("Number of Latent Factors")
    print(plotvals)
    SickleJr@WHinitials<-WHinitials
    SickleJr@lossCalcSubsample<-samp
  return(SickleJr)
}


#' @title Set W matrices and H matrix from pre-calculated values
#' @description
#' Use values calculated in the step to determine number of latent factors in the initial
#' steps for the jrSiCKLSNMF algorithm. If only a subset was calculated, this uses
#' nndsvd to calculate missing H values
#'
#' @param SickleJr A SickleJr object
#' @param d The number of desired latent factors
#' @param usebatchupdate Whether to update in series of batches or not. Use for large numbers of cells
#' @param batchsize Desired size of batches for initialization of large H values
#' @param parallel Whether you want to run this operation in parallel
#' @param seed Desired seed.
#' @param nCores Number of desired cores
#'
#' @return a SickleJr object with a filled in W and H matrices.
#' @export
#'
SetWandHfromWHinitials<-function(SickleJr,d,usebatchupdate=TRUE,batchsize=1000,parallel=TRUE,seed=NULL,nCores=detectCores()-1){
  if(is.numeric(seed)){
    set.seed(seed)
  }
  if(!usebatchupdate&parallel){
    warning("\n When not using batch updates, parallel computation is not feasible. \n")
  }
  numcells<-dim(SickleJr@normalized.count.matrices[[1]])[2]
  numcellsinit<-dim(SickleJr@WHinitials[[1]][["H"]])[1]
  latent_factor_number=SickleJr@latent.factor.elbow.values$latent_factor_number
  Hconcatmat<-matrix(data=0,nrow=numcells,ncol=d)
  #SickleJr@H<-Hconcatmat
  if(!(d%in%latent_factor_number)){
    stop("\n Initial values for W and H and have not been calculated for the specified number of latent factors. Exiting...\n")
  }
  whichd<-which(unlist(lapply(SickleJr@WHinitials,function(x) x[["d"]]==d)))
  listwh<-copy(SickleJr@WHinitials[[whichd]])

  if(numcells!=numcellsinit){
    message("\n Number of total cells differs from number of cells used in calculation for the initial plot. Using values previously
            calculated and using nndsvd initialization for the rest. \n")
    #Set values to previously calculated ones
    Hconcatmat[SickleJr@lossCalcSubsample,]=listwh$H
    #Subset for nndsvd

    bigXmatrix<-do.call(rbind,SickleJr@normalized.count.matrices)
    if(usebatchupdate){
      allcells<-1:numcells
      shuffleindices<-sample(allcells[-SickleJr@lossCalcSubsample],(numcells-numcellsinit),replace=FALSE)
      batchlist<-split(shuffleindices,ceiling(seq_along(shuffleindices)/batchsize))
      if(!parallel){cl=NULL
        }else{
          cl<<-makeCluster(nCores)
          clusterExport(cl,varlist=c("bigXmatrix","Hconcatmat","batchlist","d"),envir = environment())
          clusterEvalQ(cl,library("jrSiCKLSNMF"))
        }
      svdlist<-pblapply(batchlist, function(x) {svdH=jrSiCKLSNMF:::.nndsvd(A=bigXmatrix[,x],k=d,flag=1)
        t(svdH$H)},cl=cl)
      stopCluster(cl)
      for(i in 1:length(svdlist)){
        Hconcatmat[batchlist[[i]],]<-svdlist[[i]]
      }
    }else{
      svdH<-.nndsvd(A=bigXmatrix[,-SickleJr@lossCalcSubsample],k=d,flag=1)
      Hconcatmat[-SickleJr@lossCalcSubsample,]<-t(svdH$H)
    }
    rowReg=SickleJr@rowRegularization
    if(rowReg=="L2Norm"){
      norms<-apply(Hconcatmat,2,function(x)sqrt(sum(x^2)))
      for(i in 1:dim(Hconcatmat)[[2]]){
        Hconcatmat[,i]=Hconcatmat[,i]/norms[i]
      }
    } else if(rowReg=="L1Norm"){
      Hconcatmat<-apply(Hconcatmat,2,function(x)x/sum(x))
    } else{
      meanHconcat<-mean(apply(Hconcatmat,2,function(x) sum(x)))
      Hconcatmat<-apply(Hconcatmat,2,function(x) x/sum(x))*meanHconcat
    }
    SickleJr@H=Hconcatmat
    SickleJr@Wlist=SickleJr@WHinitials[[whichd]][["Wlist"]]
  } else{
    SickleJr@Wlist<-SickleJr@WHinitials[[whichd]][["Wlist"]]
    SickleJr@H<-SickleJr@WHinitials[[whichd]][["H"]]}
  return(SickleJr)
}

#' @title RunjrSiCKLSNMF
#' @description Wrapper function to run jrSiCKLSNMF on a SickleJr object. Performs jrSiCKLSNMF on
#' the given SickleJr
#' @name RunjrSiCKLSNMF
#' @param SickleJr a SickleJr object
#' @param rounds Number of rounds. Defaults to 300.
#' @param differr Tolerance for updating. Defaults to 1e-6.
#' @param display_progress whether to display the progress bar for jrSiCKLSNMF
#' @param lossonsubset Indicates whether you will be using a subset to calculate the loss function
#' rather than the whole dataset
#' @param losssubsetsize Size of the subset of data to calculate the loss on
#' @param online Indicates whether the algorithm will be online or not
#' @param batchsize indicates size of batch for large matrices.
#' @param random_W_updates indicates whether or not to use random_W_updates updates (i.e. only update
#' W once per online iteration)
#' @param seed Number to specify seed desired.
#' @param minrounds A minimum number of rounds to use. Most helpful for the online algorithm
#' @param suppress_warnings Indicates whether or not you want warnings suppressed
#' @param subsample Used primarily when finding an appropriate number of latent factors. Defaults
#' to all cells, but you can run jrSiCKLSNMF on just a subset of cells.
#' @return a SickleJr object with updated W matrices, updated H matrices, and a vector of values for
#' the loss function
#' @export
RunjrSiCKLSNMF<-function(SickleJr,rounds=30000,differr=1e-6,
                         display_progress=TRUE,lossonsubset=FALSE,losssubsetsize=dim(SickleJr@H)[1],
                         online=FALSE,batchsize=1000,random_W_updates=FALSE,seed=NULL,minrounds=200,
                         suppress_warnings=FALSE,subsample=1:dim(SickleJr@normalized.count.matrices[[1]])[2]){
  SickleJr@online<-FALSE
  if(online){
    SickleJr@online<-TRUE
  }
  if(is.numeric(seed)){
    set.seed(seed)
  }
  initsamp=NULL
  if(!lossonsubset){
    if(losssubsetsize!=dim(SickleJr@H)[1]){
      warning("\n Using a subset to calculate the loss was set to FALSE while the size of the subset for the loss was not equal to the number of cells. Setting the loss subset to the total number of cells... \n")
    }
    losssubsetsize=dim(SickleJr@H)[1]
  }
  if(losssubsetsize<dim(SickleJr@H)[1]){
    initsamp<-sample(1:dim(SickleJr@normalized.count.matrices[[1]])[2],losssubsetsize,replace=FALSE)
    #prepare initsamp for c++
    initsamp<-initsamp-1
  }else{initsamp=1:dim(SickleJr@H)[1]-1}
  AdjL=lapply(SickleJr@graph.laplacian.list,function(x) {x@x[which(x@x>0)]=0
                                                          return(-x)})
  DL=lapply(SickleJr@graph.laplacian.list,function(x) {x@x[which(x@x<0)]=0
                                                        return(x)})
  usemats<-lapply(SickleJr@normalized.count.matrices,function(x) x[,subsample])
  outputloss<-jrSiCKLSNMF(datamatL=usemats,
                          WL=SickleJr@Wlist,H=SickleJr@H,
                          AdjL=AdjL,
                          DL=DL,
                          lambdaWL=SickleJr@lambdaWlist,
                          lambdaH=SickleJr@lambdaH,
                          Hconstraint=toString(SickleJr@rowRegularization),
                          differr=differr,
                          display_progress = display_progress,
                          rounds=rounds,
                          diffFunc=SickleJr@diffFunc,
                          online=online,
                          batchsize=batchsize,
                          initsamp=initsamp,
                          random_W_updates=random_W_updates,
                          minrounds=minrounds,
                          suppress_warnings=suppress_warnings)
  SickleJr@loss<-outputloss
  return(SickleJr)
}

#' @title Plot a diagnostic plot for the online algorithm
#' @description
#' Plot the loss vs the number of iterations for the online algorithm. After a certain number of iterations,
#' the loss should appear to oscillate around a value.
#'
#'
#' @param SickleJr A SickleJr Object
#'
#' @return A diagnostic plot. SickleJr object remains unchanged.
#' @export

OnlineDiagnosticPlot<-function(SickleJr){
  if(!SickleJr@online){
    warning("\n This is not an online algorithm, so this plot won't give you particularly helpful diagnostics.")
  }
  niterations<-4:length(SickleJr@loss$Loss)
  loss<-SickleJr@loss$Loss[-c(1:3)]
  lossdata<-data.frame(niterations=niterations,loss=loss)
  g<-ggplot(lossdata,aes(x=niterations,y=loss))+geom_line()+theme_bw()+ggtitle("Plot of Loss vs. Number of Iterations")+xlab("Number of Iterations")+ylab("Loss")
  print(g)
}

#' @title Perform clustering diagnostics
#' @description
#' A wrapper for clValid and factoextra functions to perform clustering diagnostics
#' @param SickleJr A SickleJr object
#' @param numclusts a vector of clusters to test
#' @param clusteringmethod Clustering method. Defaults to k-means. Since other methods
#' are not implemented in jrSiCKLSNMF, it is recommended to use k-means.
#' @param diagnosticmethods Which methods to plot. Defaults to all three of the available: wss, silhouette, and gap_stat
#' @param clValidvalidation validation method to use for ClValid. Defaults to internal.
#' @param printDiagnosticplots Whether to print the diagnostic plots
#' @param printClValidDiagnostics Whether to print the diagnostic results from clValid
#'
#' @return a SickleJr object with cluster diagnostics
#' @export
determineClusters<-function(SickleJr,numclusts=2:20,clusteringmethod="kmeans",
                            diagnosticmethods=c("wss","silhouette","gap_stat"),
                            clValidvalidation="internal",printDiagnosticplots=TRUE,printClValidDiagnostics=TRUE,subset=TRUE,subsetsize=1000,seed=15){

  ingraph<-FALSE
  inclValid<-FALSE
  if((clusteringmethod%in% c("kmeans","clara","fanny","dbscan","Mclust","HCPC","hkmeans"))&
     !(clusteringmethod%in% c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes"))){
    warning("\n This clustering method is not available in clValid. Only graphs will be printed and saved in the clusterdiagnostics node.\n")
    ingraph<-TRUE
    }else if(!(clusteringmethod%in% c("kmeans","clara","fanny","dbscan","Mclust","HCPC","hkmeans"))&
           (clusteringmethod%in% c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes"))){
      warning("\n This clustering method is not available in factoextra. Only information from clValid will be printed and saved in the clusterdiagnostics node.\n")
      inclValid<-TRUE
  }else if ((clusteringmethod%in% c("kmeans","clara","fanny","dbscan","Mclust","HCPC","hkmeans"))&
            (clusteringmethod%in% c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes"))){
    ingraph<-TRUE
    inclValid<-TRUE
  }else{
    stop("\n Valid clustering method not selected. Please select a valid clustering method. Options include:\n 'kmeans','clara','fanny','dbscan','Mclust',HCPC','hkmeans','hierarchical','diana','som','model','pam',and 'agnes'.")
  }
  ggplotlist<-list()
  rownames(SickleJr@H)<-colnames(SickleJr@count.matrices[[1]])
  clValidobj<-NULL
  Hrun<-NULL
  if(subset&subsetsize>dim(SickleJr@H)[1]){
    warning("Size of subset is greater than the number of cells. Using full dataset")
    subset=FALSE
  }
  if(subset){
    if(is.numeric(seed)){
      set.seed(seed)
    }
    subset<-sample(1:dim(SickleJr@H)[1],subsetsize,replace = FALSE)
    Hrun<-SickleJr@H[subset,]
  }else{
    Hrun<-SickleJr@H
  }
  if(ingraph){
    ggplotlist<-lapply(diagnosticmethods,function(x) fviz_nbclust(Hrun,match.fun(clusteringmethod),method=x,k.max=max(numclusts)))
    if(printDiagnosticplots){
      lapply(ggplotlist,function(x) print(x))
    }
  }
  if(inclValid){

    clValidobj<-clValid(Hrun,clMethods=clusteringmethod,validation=clValidvalidation,nClust=numclusts)
    if(printClValidDiagnostics){
      print(summary(clValidobj))
    }
  }
  SickleJr@clusterdiagnostics<-list(DiagnosticPlots=ggplotlist,clValidResults=clValidobj)
  return(SickleJr)
}

#' @title ClusterSickleJr
#' @description Perform either K-means or spectral clustering on the H matrix.
#' Defaults to kmeans.
#' @name ClusterSickleJr
#' @param SickleJr an object of class SickleJr
#' @param numclusts Number of clusters
#' @param method Clustering method. Can choose between kmeans and spectral clustering
#' @param neighbors parameter for spectral clustering
#'
#' @return a SickleJr object with added clustering information
#' @export
ClusterSickleJr<-function(SickleJr,numclusts,method="kmeans",neighbors=NULL){
  clust<-NULL
  if(method=="kmeans"){
   clust<-kmeans(SickleJr@H,centers=numclusts,nstart=1000)$cluster
  }else if(method=="spectral"){
    clust<-specClust(SickleJr@H,centers=numclusts,nn=neighbors)$cluster
  }else {print("Please enter 'kmeans' for k-means clustering or 'spectral' for spectral clustering")}
  SickleJr@clusters[[method]]<-clust
  return(SickleJr)
}


#' @title Calculate the UMAP for a SickleJr object.
#' @description Perform UMAP on the H matrix alone (default) or within a view by
#' using UMAP on the \eqn{\textbf{W}^v}H corresponding to view \eqn{v}
#' @name CalculateUMAPSickleJr
#' @param SickleJr An object of class SickleJr
#' @param umap.settings Optional settings for UMAP. Set to default.
#' @param view If this parameter is set, this will perform UMAP on W_view%*%t(H) rather than
#' on H alone. Warning that this takes a very long time for large datasets and is
#' not recommended for them.
#'
#' @return a SickleJr with added UMAP calculations based on the H matrix alone or within a view
#' @export
CalculateUMAPSickleJr<-function(SickleJr,umap.settings=umap::umap.defaults,view=NULL){
  if(!is.null(view)){
    WH<-t(SickleJr@Wlist[[view]]%*%t(SickleJr@H))
    WHname<-paste0("W",view,"H")
    SickleJr@umap[[WHname]]<-umap::umap(WH,config=umap.settings)
  }else{
    SickleJr@umap[["H"]]<-umap::umap(SickleJr@H,config=umap.settings)
  }
  return(SickleJr)
}

#' @title AddSickleJrMetadata
#' @description Add metadata to a SickleJr object
#' @name AddSickleJrMetadata
#' @param SickleJr a SickleJr object
#' @param metadata The metadata you wish to add to the object
#' @param metadataname The name that you wish to call the added metadata
#' @return SickleJr with added metadata
#' @export
AddSickleJrMetaData<-function(SickleJr,metadata,metadataname){
  SickleJr@metadata[[metadataname]]<-metadata
  return(SickleJr)
}

#' @title PlotSickleJrUMAP
#' @description Plot the first and second dimensions of a UMAP dimension reduction
#' and color either by cluster or metadata
#' @name PlotSickleJrUMAP
#' @param SickleJr An object of type SickleJr
#' @param umap.view String corresponding to the name of the UMAP you are interested in. Defaults to H
#' @param cluster Cluster you wish to color by
#' @param title Optional title for your plot
#' @param colorbymetadata Name of metadata column if you wish to color by metadata
#' @param legendname Option if you want to give a different name for your legend
#'
#' @return A ggplot2 object of the plots of UMAP1 and UMAP2 for the SickleJr object
#' @export
PlotSickleJrUMAP<-function(SickleJr,umap.view="H",cluster="kmeans",title="",colorbymetadata=NULL,legendname=NULL){
  if(is.null(colorbymetadata)){
    color=SickleJr@clusters[[cluster]]
    if(is.null(legendname)){
      legendname=paste(cluster, "cluster")
      }

  }else{
    color=SickleJr@metadata[[colorbymetadata]]
    if(is.null(legendname)){legendname=colorbymetadata}
  }
  UMAP1=SickleJr@umap[[umap.view]]$layout[,1]
  UMAP2=SickleJr@umap[[umap.view]]$layout[,2]
  umapvals<-data.frame(UMAP1=UMAP1,UMAP2=UMAP2,cluster=color)
  ggplot(umapvals,aes(x=UMAP1,y=UMAP2,color=factor(color)))+geom_point()+theme_bw()+
    ggtitle(title)+guides(color=guide_legend(title=legendname))
}
