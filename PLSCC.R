# Author: Pierre Maho
# Lab: Gipsa-Lab
# PhD Supervisors: Simon Barthelme and Pierre Comon
# Date: 2020/01/12


###########################
# Partial least squares
PLSCC_tune <- function(X,y,ncomp=2){
  # Partial Least Squares - Component Correction from Artursson et al.
  # Tom Artursson, Tomas Eklöv, Ingemar Lundström, et al. “Drift correction for gas sensors using multivariate methods”. In: Journal of Chemometrics (Sept. 2000)
  #
  # Parameters
  # X: Calibrant data matrix (n_obs x n_sens)
  # y: time (n_obs x 1)
  # ncomp: Number components to return
  #
  # Return
  # ncomp drift latent directions
  res_pls <- pls::plsr(y~X)
  return(list(weights = matrix(res_pls$loading.weights[,1:ncomp],ncol = ncomp),
              loadings = matrix(res_pls$loadings[,1:ncomp],ncol = ncomp)))
}

PLSCC_cor <- function(X,pls,ncomp=2){
  # Partial Least Squares - Component Correction from Artursson et al.
  # Tom Artursson, Tomas Eklöv, Ingemar Lundström, et al. “Drift correction for gas sensors using multivariate methods”. In: Journal of Chemometrics (Sept. 2000)
  #
  # Parameters
  # X: Drifted data matrix (n_obs x n_sens)
  # pls: PLS object containing the latent directions to remove
  #
  # Return
  # Corrected data matrix
  Xcor <- X
  ncomp <- ncol(pls$loadings)
  for(i in 1:ncomp){
    wd <- matrix(pls$weights[,i])
    pd <- matrix(pls$loadings[,i])
    td <- Xcor %*% wd %*% ginv(t(pd) %*% wd)
    Xcor <- Xcor - td %*% t(pd)
  }
  return(Xcor)
}
###########################