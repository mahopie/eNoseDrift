# Author: Pierre Maho
# Lab: Gipsa-Lab
# PhD Supervisors: Simon Barthelme and Pierre Comon
# Date: 2020/01/12


###########################
# Common Principal Component Analysis
CPCA_tune <- function(X,y, ncomp = 2){
  # Common Principal Component Analysis - Component Correction from A. Ziyatdinov et al.
  # A. Ziyatdinov, S. Marco, A. Chaudry, et al. “Drift compensation of gas sensor array data by common principal component analysis”. In: Sensors and Actuators B: Chemical (Apr. 2010)
  #
  # Parameters
  # X: Multi-Session data matrix ((n_obs1 + ... + n_obsS) x n_sens)
  # y: Multi-Session label vector ((n_obs1 + ... + n_obsS) x 1)
  # ncomp: Number components to return
  #
  # Return
  # ncomp drift directions
  
  # Compute second order moment for each class
  labels = unique(y)
  covMats = array(0,c(ncol(X),ncol(X),length(labels)))
  for(i in 1:length(labels)){
    covMats[,,i] = cov(X[y==labels[i],])
  }
  
  # Joint diagonalization 
  jd <- JADE::cjd(covMats)
  
  if(sum(Im(jd$V))!=0)
    print('Warning: imaginary part of loadings is non-null')
  loadings <- Re(jd$V)
  eig <- Re(jd$D)
  
  # Reorder the eigenvalues (decreasing)
  SSE = numeric(ncol(X))
  for(i in 1:ncol(X)){
    for(j in 1:length(labels)){
      SSE[i] = SSE[i] + sum((covMats[,,j] - eig[i,i,j] * loadings[,i] %*% t(loadings[,i]))**2)
    }
  }
  SSE_ord = sort(SSE, index.return=T)$ix
  
  return(matrix(loadings[,SSE_ord[1:ncomp]],ncol=ncomp))
}
CPCA_cor <- function(X, loadings){
  # Common Principal Component Analysis - Component Correction from A. Ziyatdinov et al.
  # A. Ziyatdinov, S. Marco, A. Chaudry, et al. “Drift compensation of gas sensor array data by common principal component analysis”. In: Sensors and Actuators B: Chemical (Apr. 2010)
  #
  # Parameters
  # X: drift data matrix (n_obs x n_sens)
  # loadings: CPCA drift directions (n_sens x ncomp)
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