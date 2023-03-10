NULL

#' @title SickleJr
#' @description defines the SickleJr class
#'
#' @slot count.matrices list. List of all of the QC'd count matrices.
#' Note that these count matrices should not use
#' @slot normalization.type character. Holds the type of normalization
#' @slot normalized.count.matrices list. Holds the normalized count matrices
#' @slot graph.laplacian.list list. A list of the graph laplacians used for graph regularizations
#' @slot rowRegularization character. A string that indicates the type of row regularization to
#' use. Types include "None" and "L2Norm"
#' @slot Wlist list. A list of the generated W matrices, one for each view
#' @slot H matrix. The shared H matrix
#' @slot diffFunc character. Holds the name of the function used to measure the "distance" between
#' data matrix X and WH for each view. Can be "klp" for Kullback-Leibler based on the Poisson
#' distribution or "frob" for the Frobenius norm.
#' @slot lambdaWlist list. List of lambda values used, one for H view
#' @slot lambdaH numeric. Value for the constraint on H
#' @slot clusters list. List of different clusters performed on the SickleJr object
#' @slot metadata list. List of metadata
#' @slot loss vector. Vector of the value for the loss function
#' @slot umap list. List of different UMAP based dimension reductions
#'
#' @return an object of type SicklJr
#' @export
SickleJr<-methods::setClass(
  "SickleJr",
  slots=c(
    count.matrices="list",
    normalization.type="character",
    normalized.count.matrices="list",
    graph.laplacian.list="list",
    rowRegularization="character",
    Wlist="list",
    H="matrix",
    diffFunc="character",
    lambdaWlist="list",
    lambdaH="numeric",
    clusters="list",
    metadata="list",
    loss="vector",
    umap="list"
  )
)

#' @title createSickleJr
#' @name createSickleJr
#' @description creates an object of type SickleJr
#' @param count.matrices A count matrix that has already been filtered such that
#' all genes and peaks appear in at least 10 cells.
#' @param dense Indicates whether you want the count matrices stored densely or
#' sparsely
#' @return A SickleJr object with sparse or dense count matrices
createSickleJr<-function(count.matrices,dense=FALSE){
  if(!dense){
    count.matrices<-.sparsify(count.matrices)
  } else{
    count.matrices<-as.matrix(count.matrices)
  }
  object<-methods::new(Class="SickleJr",count.matrices=count.matrices)
}


#' @title BuildKNNGraphLaplacians
#' @description Generate graph Laplacians for graph regularization of
#' jrSiCKLSNMFNMF from the list of count matrices
#'@name BuildKNNGraphLaplacians
#'@param SickleJr A SickleJr object
#'@param k Number of KNN neighbors to calculate. By default, is set to 20
#'@param dense Indicates whether or not you want a dense matrix. Defaults to
#'FALSE
#'@returns A list of graph Laplacians in sparse matrix format
#'@export
BuildKNNGraphLaplacians<-function(SickleJr,k=20,dense=FALSE){
  if(!dense){
    SickleJr@graph.laplacian.list<-lapply(SickleJr@count.matrices,function(x){
      laplacian_matrix(buildKNNGraph(x,transposed=TRUE,k=k))})} else{
        SickleJr@graph.laplacian.list<-lapply(SickleJr@count.matrices,function(x){
          as.matrix(laplacian_matrix(buildKNNGraph(x,transposed=TRUE,k=k)))})
      }
  SickleJr@graph.laplacian.list<-.sparsify(SickleJr@graph.laplacian.list)
  return(SickleJr)
}

#' @title NormalizeCountMatrix
#' @description Perform normalization for count data. If you are planning to use the Frobenius norm,
#'set frob=TRUE to log(x+1) normalize your count data. This step can be skipped for percentage data and
#'spectra data. You may also skip this if you would like to perform a different form
#'of normalization; however, please ensure that if using the Frobenius norm that all of the values are
#'divided by the sum of all of the values.
#'@name NormalizeCountMatrix
#'@param SickleJr An object of type SickleJr with a count matrix. Users can choose to normalize using median
#'library size normalization or log(x+1) normalization
#'@param frob A Boolean. Set to TRUE if you want to perform log(x+1) normalization and FALSE for
#'median library size normalization as per Zhang 2017.
#'@returns A list of sparse, normalized matrices
#'@export
NormalizeCountMatrix<-function(SickleJr, frob=FALSE){
  XmatrixList<-SickleJr@count.matrices
  NormalizedXMatrices<-list()
  if(frob==TRUE){
    for(i in 1:length(XmatrixList)){
      NormalizedXMatrices[[i]]<-.sparsify(log(XmatrixList[[i]]+1))
      NormalizedXMatrices[[i]]<-NormalizedXMatrices[[i]]/sum(NormalizedXMatrices[[i]])
    }

    } else{
      medianLibSize<-list()
      for(i in 1:length(XmatrixList)){
        medianLibSize[[i]]<-median(apply(XmatrixList[[i]],2,sum))
        NormalizedXMatrices[[i]]<-apply(XmatrixList[[i]],2,
                                        function(x) x/sum(x))*medianLibSize[[i]]
    }
    }
  normalizationtype<-"Median Library Size"
  if(frob==TRUE){
    normalizationtype<-"Log(x+1) normalization"
  }
  SickleJr@normalized.count.matrices<-.sparsify(NormalizedXMatrices)
  SickleJr@normalization.type<-normalizationtype
  return(SickleJr)
}

#' @title GenerateWmatricesandHmatrix
#' @description Perform normalization for count data. If you are planning to use the Frobenius norm,
#'set frob=TRUE
#'@name GenerateWmatricesandHmatrix
#' @param SickleJr a SickleJr object
#' @param d number of latent factors to use. Defaults to 10.
#' @param rowReg Choose whether or not to enforce constraints on the rows. Enter "None" for no constraints
#'  Defaults to L2Norm.
#' @returns SickleJr A SickleJr object with W and H matrices added.
#' @export
#'
GenerateWmatricesandHmatrix<-function(SickleJr,d=10,rowReg="L2Norm"){
  WandHlist<-list()
  NormalizedXmatrices<-SickleJr@normalized.count.matrices
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
  SickleJr@rowRegularization=rowReg
  return(SickleJr)
}

#' @title setLambdas
#' @description Set the lambda values for a SickleJr object
#' @name setLambdas
#' @param SickleJr a SickleJr object
#' @param lambdaWlist A list of graph regularization constraints for the W matrices.
#' @param lambdaH lambda H value to add. Defaults to 0.
#' @return SickleJr object with added lambda values
#' @export
#'
setLambdas<-function(SickleJr,lambdaWlist,lambdaH=0.0){
  SickleJr@lambdaWlist<-lambdaWlist
  SickleJr@lambdaH<-lambdaH
  return(SickleJr)
}

#' @title RunjrSiCKLSNMF
#' @description Wrapper function to run jrSiCKLSNMF on a SickleJr object. Performs jrSiCKLSNMF on
#' the given SickleJr
#' @name RunjrSiCKLSNMF
#' @param SickleJr a SickleJr object
#' @param rounds Number of rounds. Defaults to 300.
#' @param differr Tolerance for updating. Defaults to 1e-6.
#' @param diffFunc choose between "klp" for Poisson-based kullback leibler divergence
#' @param display_progress whether to display the progress bar for jrSiCKLSNMF
#'
#' @return a SickleJr object with updated W matrices, updated H matrices, and a vector of values for
#' the loss function
#' @export
RunjrSiCKLSNMF<-function(SickleJr,rounds=30000,differr=1e-6,diffFunc="klp",display_progress=TRUE){
  outputloss<-jrSiCKLSNMF(datamatL=SickleJr@normalized.count.matrices,WL=SickleJr@Wlist,H=SickleJr@H,
                                AL=SickleJr@graph.laplacian.list,lambdaWL=SickleJr@lambdaWlist,
                                lambdaH=SickleJr@lambdaH, Hconstraint=toString(SickleJr@rowRegularization),
                                differr=differr,display_progress = display_progress,diffFunc=diffFunc)
  SickleJr@loss<-outputloss
  return(SickleJr)
}


#' @title ClusterSickleJr
#' @description perform either kmeans or spectral clustering on the H matrix.
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
    clust<-kknn::specClust(SickleJr@H,centers=numclusts,nn=neighbors)$cluster
  }else {print("Please enter 'kmeans' for k-means clustering or 'spectral' for spectral clustering")}
  SickleJr@clusters[[method]]<-clust
  return(SickleJr)
}


#' @title CalculateUMAPSickleJr
#' @description Perform UMAP on the H matrix alone (default) or within a view
#' @name CalculateUMAPSickleJr
#' @param SickleJr An object of class SickleJr
#' @param umap.settings Optional settings for UMAP. Set to default.
#' @param view If this parameter is set, this will perform UMAP on W_view%*%t(H) rather than
#' on H alone.
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
#'
#' @return SickleJr with added metadata
#' @export
AddSickleJrMetaData<-function(SickleJr,metadata,metadataname){
  SickleJr@metadata[[metadataname]]<-metadata
  return(SickleJr)
}

#' @title plotSickleJrUMAP
#' @description Plot first and second dimensions of UMAP and color either by cluster or metadata
#' @name plotSickleJrUMAP
#' @param SickleJr An object of type SickleJr
#' @param umap.view String corresponding to the name of the UMAP you are interested in. Defaults to H
#' @param cluster Cluster you wish to color by
#' @param title Optional title
#' @param colorbymetadata Name of metadata column if you wish to color by metadata
#' @param legendname Option if you want to give a different name for your legend
#'
#' @return A ggplot2 object of the plots of UMAP1 and UMAP2 for the SickleJr object
#' @export
plotSickleJrUMAP<-function(SickleJr,umap.view="H",cluster="kmeans",title="",colorbymetadata=NULL,legendname=NULL){
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
