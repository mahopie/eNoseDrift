# Author: Pierre Maho
# Lab: Gipsa-Lab
# PhD Supervisors: Simon Barthelme and Pierre Comon
# Date: 2020/01/12


###########################
# Principal Component Analysis
PCACC_tune <- function(X,ncomp=NULL,varExp=0.95){
  # Principal Component Analysis - Component Correction from Artursson et al.
  # Tom Artursson, Tomas Eklöv, Ingemar Lundström, et al. “Drift correction for gas sensors using multivariate methods”. In: Journal of Chemometrics (Sept. 2000)
  #
  # Parameters
  # X: Calibrant data matrix (n_obs x n_sens)
  # ncomp: Number components to return
  # varExp: minimal number of components to explain varExp variance in the calibrant matrix
  #
  # Return
  # ncomp drift directions
  pca <- prcomp(X)
  if(is.null(ncomp)){
    var_comp = cumsum(pca$sdev**2)/sum(pca$sdev**2)
    ncomp = which(var_comp > varExp)[1]
    print(paste(ncomp,'components have been chosen'))
  }
  return(matrix(pca$rotation[,1:ncomp],ncol=ncomp))
}

PCACC_cor <- function(X, loadings){
  # Principal Component Analysis - Component Correction from Artursson et al.
  # Tom Artursson, Tomas Eklöv, Ingemar Lundström, et al. “Drift correction for gas sensors using multivariate methods”. In: Journal of Chemometrics (Sept. 2000)
  #
  # Parameters
  # X: drifted data matrix (n_obs x n_sens)
  # loadings: ncomp PCs for correction
  #
  # Return
  # Drift corrected data matrix
  return(X%*%(diag(rep(1,ncol(X))) - loadings %*% t(loadings)))
}
###########################
