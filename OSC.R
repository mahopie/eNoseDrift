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

OSC_tune <- function(X,y,ncomp=2){
  # Orthogonal Signal Correction from M. Padilla et al.
  # M. Padilla, A. Perera, I. Montoliu, et al. “Drift compensation of gas sensor array data by Orthogonal Signal Correction”. In: Chemometrics and Intelligent Laboratory Systems (Jan. 2010)
  #
  # Parameters
  # X: Multi-Session data matrix ((n_obs1 + ... + n_obsS) x n_sens)
  # y: Multi-Session label vector ((n_obs1 + ... + n_obsS) x 1)
  # ncomp: Number components to return
  #
  # Return
  # ncomp drift latent directions
  Y <- indicator_matrix(y)
  W <- matrix(0, ncol = ncomp, nrow = ncol(X))
  loadings <- matrix(0, ncol = ncomp, nrow = ncol(X))
  nruns <- 100
  
  for(n in 1:ncomp){
    pca <- prcomp(X)
    t <- pca$x[,1]
    i <- 0
    err <- 1
    while(i < nruns & err > 1e-6){
      t_old <- t
      nt = (diag(rep(1,nrow(Y))) - Y %*% ginv(t(Y) %*% Y) %*% t(Y)) %*% t
      w <- ginv(X) %*% nt
      t <- X %*% w 
      err <- (sqrt(sum((t-t_old)**2)) / sqrt(sum((t)**2)))
      i <- i+1
    }
    p <- t(t(t) %*% X / as.numeric(t(t)%*%t))
    loadings[,n] <- p
    W[,n] <- w
    X <- X - t %*% t(loadings[,n])
  }
  return(list(weights=W,loadings=loadings))
}

OSC_cor <- function(X,osc){
  # Orthogonal Signal Correction from M. Padilla et al.
  # M. Padilla, A. Perera, I. Montoliu, et al. “Drift compensation of gas sensor array data by Orthogonal Signal Correction”. In: Chemometrics and Intelligent Laboratory Systems (Jan. 2010)
  #
  # Parameters
  # X: drift data matrix (n_obs x n_sens)
  # osc: OSC object containing drift latent directions (n_sens x ncomp)
  #
  # Return
  # Drift corrected data matrix 
  Xcor <- X
  for(i in 1:ncol(osc$weights)){
    scores <- Xcor %*% osc$weights[,i] 
    Xcor <- Xcor - scores %*% t(osc$loadings[,i])
  }
  return(Xcor)
}
###########################
