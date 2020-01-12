# Author: Pierre Maho
# Lab: Gipsa-Lab
# PhD Supervisors: Simon Barthelme and Pierre Comon
# Date: 2020/01/12


cost_function_mahalanobis <- function(M, X, pred, mus, covMats){
  # Compute Mahalanobis cost function
  #
  # Parameters
  # M: drift matrix (n_sens x n_sens)
  # X: data matrix (n_obs x n_sens)
  # pred: labels vector with integer from 1 to n_class (n_obs x 1)
  # mus: list of class centroids (n_class elements of size n_sens x 1)
  # covMats: list of class covariance matrices (n_class elements of size n_sens x n_sens)
  #
  # Note
  # If M is a matrix, then M must be vectorized col-by-col (you can call: as.vector(M))
  res <- 0
  M = matrix(M, ncol = ncol(X), nrow = ncol(X))
  for(i in 1:nrow(X)){
    tmp =  t(matrix(X[i,])) %*%  M - mus[[pred[i]]] 
    tmp = matrix(tmp)
    res = res + t(tmp) %*% ginv(covMats[[pred[i]]]) %*% tmp
  }
  return(res)
}

evolutionaryBased_DiCarlo <- function(X1, y1, X2, Minit=NULL, optimMethod =c('optim', 'cmaes')){
  # Drift correction method of Di Carlo et al.
  # S. Di Carlo, M. Falasconi, E. Sanchez, et al. “Increasing pattern recognition accuracy for chemical sensing by evolutionary based drift compensation”. In: Pattern Recognition Letters (Oct. 2011)
  #
  # Parameters
  # X1: training data matrix (n_obs1 x n_sens)
  # y1: training label vector with integer from 1 to n_class (n_obs1 x 1)
  # X2: drifted data matrix (n_obs2 x n_sens)
  # optimMethod: optimization method used to find M. In their article, they used CMAES but it can be slow (several minutes), so a quicker optim is based on BFGS.
  # Minit: initial drift matrix (n_sens x n_sens)
  #
  # Return
  # X2cor: corrected version of X2 (n_obs2 x n_sens)
  # M: drift matrix (n_sens x n_sens)
  # predCor: predicted labels after correction
  
  if(is.null(Minit)){
    Minit <- diag(rep(1,ncol(X1)))
  }
  M <- Minit
  
  # Extract first and second moment of each class
  mus = covMats = list()
  labels = unique(y1)
  for(i in 1:length(labels)){
    mus[[i]] = colMeans(X1[y1==labels[i],])
    covMats[[i]] = cov(X1[y1==labels[i],])
  }
  
  # Train linear-SVM classifier (for other classifier, modify this part)
  svm_tune <- tune(svm, train.x=X1, train.y=factor(y1), kernel='linear', scale = F, ranges=list(cost=c(10^(-5:5))))
  
  # First classification step based on Raw drifted data X2
  predRaw <- predict(svm_tune$best.model, X2 %*% M) %>% as.numeric
  
  # Adapt the drift matrix M
  if(optimMethod == "optim"){
    optim_res = optim(as.vector(M), cost_function_mahalanobis, X = X2, pred = predRaw, mus = mus, covMats = covMats,method="BFGS")
  }else{
    optim_res = cmaes::cma_es(par = as.vector(M), cost_function_mahalanobis, X = X2, pred = predRaw, mus = mus, covMats = covMats)
  }
  
  # Correct and pred with the estimated M
  M = matrix(optim_res$par, ncol = ncol(X1), nrow = ncol(X1))
  X2cor = X2 %*% M
  predCor <- predict(svm_tune$best.model,X2cor) %>% as.numeric
  
  return(list(X2cor=X2cor, M=M, pred=predCor))
}