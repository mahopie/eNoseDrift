# Author: Pierre Maho
# Lab: Gipsa-Lab
# PhD Supervisors: Simon Barthelme and Pierre Comon
# Date: 2020/01/12


###########################
indicator_matrix <- function(y){
  # Transform a label vector into a label matrix
  M <- matrix(0, ncol = length(unique(y)), nrow=length(y))
  cl <- unique(y)
  for(i in 1:length(cl))
    M[which(y==cl[i]),i] <- 1
  return(M)
}

OSC_DO_tune <- function(X,y,ncomp=2){
  # Orthogonal Signal Correction from M. Padilla et al.
  # M. Padilla, A. Perera, I. Montoliu, et al. “Drift compensation of gas sensor array data by Orthogonal Signal Correction”. In: Chemometrics and Intelligent Laboratory Systems (Jan. 2010)
  #
  # There are multiple ways for performing OSC, here we implement the solution of C. A. Andersson, called Direct Orthogonalization (DO)
  # C. A. Andersson, "Direct orthogonalization". In: Chemometrics and Intelligent Laboratory Systems (1999)
  #
  # Parameters
  # X: Multi-Session data matrix ((n_obs1 + ... + n_obsS) x n_sens)
  # y: Multi-Session label vector ((n_obs1 + ... + n_obsS) x 1)
  # ncomp: Number components to return
  #
  # Return
  # ncomp drift latent directions
  #
  # Useful references
  # O. Svensson, T. Kourti and J. F. MacGregor, "An investigation of orthogonal signal correction algorithms and their characteristics". In: Journal of Chemometrics (2002)
  
  Y <- indicator_matrix(y)
  Z <- X - Y %*% ginv(t(Y) %*% Y) %*% t(Y) %*% X
  pca <- prcomp(Z)
  return(as.matrix(pca$rotation[,1:ncomp], ncol=ncomp))
}

OSC_DO_cor <- function(X,loadings){
  # Orthogonal Signal Correction (OSC) from M. Padilla et al.
  # M. Padilla, A. Perera, I. Montoliu, et al. “Drift compensation of gas sensor array data by Orthogonal Signal Correction”. In: Chemometrics and Intelligent Laboratory Systems (Jan. 2010)
  #
  # Parameters
  # X: drift data matrix (n_obs x n_sens)
  # loadings: n components for correction (n_sens x ncomp)
  #
  # Return
  # Drift corrected data matrix 
  return(X%*%(diag(rep(1,ncol(X))) - loadings %*% t(loadings)))
}
###########################
