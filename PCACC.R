# Author: Pierre Maho
# Lab: Gipsa-Lab
# PhD Supervisors: Simon Barthelme and Pierre Comon
# Date: 2020/01/12


###########################
# Principal Component Analysis
PCACC_tune <- function(X,ncomp=2){
  # Principal Component Analysis - Component Correction from Artursson et al.
  # Tom Artursson, Tomas Eklöv, Ingemar Lundström, et al. “Drift correction for gas sensors using multivariate methods”. In: Journal of Chemometrics (Sept. 2000)
  #
  # Parameters
  # X: Calibrant data matrix (n_obs x n_sens)
  # ncomp: Number components to return
  #
  # Return
  # ncomp drift directions
  pca <- prcomp(X)
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
  Xcor <- X
  ncomp <- ncol(loadings)
  for(i in 1:ncomp){
    t <- Xcor %*% loadings[,i]
    Xcor <- Xcor - t%*% t(loadings[,i])
  }
  return(Xcor)
}
###########################