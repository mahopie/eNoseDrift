####################################################################
# Author : Pierre Maho
# 24/07/2020
# GIPSA-Lab
# Expectation-Maximization Component Correction

within_class_covariance <- function(X,y){
  # Compute within-class scatter matrix
  Sw <- matrix(0, nrow=ncol(X), ncol=ncol(X))
  ynum <- as.numeric(as.factor(y))
  
  for(i in 1:max(ynum)){
    tmp <- which(ynum==i)
    mui <- as.matrix(colMeans(X[tmp,]))
    for(j in 1:length(tmp)){
      xj <- as.matrix(X[tmp[j],])
      Sw <- Sw + (xj-mui) %*% t(xj-mui)
    }
  }
  return(Sw)
}

init_EMC2 <- function(X1, y1, X, N_ptX=10){
  # Initialize d, alpha and pi for EM algorithm (drift compensation under the Blind Scenario)
  #
  # Take N_ptX random points from X
  # Compute directions from X1 centroids to these random points
  # For each direction :
  #        Project X1 and X in the orthogonal subspace
  #        Learn a kNN using (X1, y1) and predict X
  #        Compute the trace of the within-class covariance matrix
  # Choose the direction with the minimal trace
  # According to the labels, compute pi and alpha
  #
  # return a list containing alpha, d, pi and the labels.
  
  labels <- unique(y1)
  R <- length(labels)
  P <- ncol(X)
  N <- nrow(X)
  
  ##### ESTIMATE THE CENTROIDS OF X1
  mus <- list()
  for(p in 1:R){
    mus[[p]] <- colMeans(X1[y1==labels[p],])
  }
  
  ##### SELECT RANDOM DIRECTIONS
  p <- 1
  ds <- matrix(0, nrow = P, ncol = N_ptX*R)
  ptX <- sample(1:N)[1:N_ptX]
  for(j in 1:R){
    for(i in 1:N_ptX){
      ds[,p] <- X[ptX[i],] - mus[[j]]
      ds[,p] <- ds[,p]/sqrt(sum(ds[,p]**2)) # normalize the directions
      p <- p+1
    }
  }
  
  # add "trivial" init for direction 
  ds <- cbind(ds, (colMeans(X) - colMeans(X1)) / sqrt(sum((colMeans(X) - colMeans(X1))**2)))
  
  ##### ORTHOGONALIZE AND SELECT THE BEST CANDIDATE
  WVs <- numeric(ncol(ds))
  for(p in 1:(ncol(ds))){
    d <- matrix(ds[,p])
    d_ortho <- diag(rep(1,P)) - d %*% t(d)
    X1p <- X1 %*% d_ortho
    Xp <- X %*% d_ortho
    
    pred_knn <- class::knn(X1p, Xp, y1,k=3)
    WVs[p] <- sum(diag(within_class_covariance(rbind(X1p,Xp),c(y1,pred_knn))))
  }
  
  # Estimate d as the one minimizing the sum of the distances in the orthogonal subspace
  d <- matrix(ds[,which.min(WVs)])
  
  #####  ESTIMATE alpha and pi
  d_ortho <- diag(rep(1,P)) - d %*% t(d)
  X1p <- X1 %*% d_ortho
  Xp <- X %*% d_ortho
  pred_knn <- class::knn(X1p, Xp, y1,k=3)
  
  pi <- numeric(R)
  alpha <- numeric(R)
  for(r in 1:R){
    which_r = which(pred_knn == labels[r])
    if(length(which_r)==0){
      pi[r] <- 0
      alpha[r] <- 0
    }else{
      pi[r] <- length(which_r)/N
      muXr <- colMeans(matrix(X[which_r,],ncol = P)) %>% matrix
      alpha[r] <- t(d) %*% (muXr - mus[[r]])
    }
  }
  return(list(d=d, alpha=alpha, pi=pi, class=pred_knn))
}


plot_EMC2 <- function(X1, y1, X, y, mus, covMats,i){
  # plot function for EMC2 algorithm
  R = length(covMats)
  plot(X1[,1], X1[,2], col = factor(y1), pch = 1, main=paste("Iteration ",i), xlim = c(min(rbind(X1,X)[,1]),max(rbind(X1,X)[,1])), ylim = c(min(rbind(X1,X)[,2]),max(rbind(X1,X)[,2])))
  points(X[,1], X[,2],  col = factor(y),pch = 21, lwd = 2)  
  for(k in 1:R)
    ellipse(mus[[k]][1:2], covMats[[k]][1:2,1:2], alpha = .05, npoints = 250, newplot = FALSE, draw = TRUE)
}

EMC2 <- function(X1, y1, X, N_ptX = 10, covDiag = F, covSame = F, init_em = NULL, plot=F, y=NULL){
  # Expectation-Maximization Component Correction from Maho et al. (2021)
  # EMC2 is a drift correction method which assumes that each class drifts along a common direction weighted by factor depending on the class. 
  # Input:
  #   - X1: data matrix (N1 x P) of Session 1 (training set)
  #   - y1: labels vector (N1 x 1) of Session 1 
  #   - X: data matrix (N x P) containing another session 
  #   - N_ptX: number of points used for the initialization
  #   - covDiag: use diagonal covariance matrices (suitable for high dimension)
  #   - covSame: use the same covariance matrix for all the classes (suitable for high dimension)
  #   - init_em: EM list containing initial values for alpha, d, pi and the classes of X
  #   - plot: each iteration is plotted
  #   - y: if each iteration is plotted, function needs the labels of X
  # Output:
  #   - a list containing the corrected data, the predicted classes, the final em object and all the em objects estimated during the maximization
  labels <- unique(y1)
  R <- length(labels) # number of compounds in Session 1
  P <- ncol(X) # number of chemical sensors
  N <- nrow(X) # number of samples in Session 2
  
  ##### ESTIMATION OF GAUSSIAN PARS
  mus <- covMats <- covMats_inv <- list()
  for(p in 1:R){
    mus[[p]] <- colMeans(X1[y1==labels[p],])
    if(covDiag){ # Diagonal covariance matrices
      covMats[[p]] <- diag(colSd(X1[y1==labels[p],])**2)
    }else{
      covMats[[p]] <- cov(matrix(X1[y1==labels[p],], ncol=P))
    }
  }
  if(covSame){ # Same covariance matrix for all classes
    tmp = 0
    for(p in 1:R)
      tmp <- tmp + covMats[[p]]
    covMats <- lapply(covMats, function(x) x = tmp / R)
  }
  # Compute the inverse of the covariance matrices
  for(r in 1:R){
    covMats_inv[[r]] <- MASS::ginv(covMats[[r]])
  }
  
  
  ##### INITIALIZATION OF EM ALGO
  if(!is.null(init_em)){
    d_next <- init_em$d
    alpha <- init_em$alpha
    pi <- init_em$pi
    class <- as.numeric(init_em$class)
  }else{
    initEM <- init_EMC2(X1, y1, X, N_ptX)
    d_next <- initEM$d
    pi <- initEM$pi
    alpha <- initEM$alpha
    class <- as.numeric(initEM$class)
  }

  pi_min <- 1/nrow(X) 
  post <- matrix(0, nrow = N, ncol = R)
  eps <- 1e-5 # convergence threshold (between two successive iterations)
  imax <- 100 # max iterations
  i <- 1
  crit <- 1
  alpha_all <- matrix(NA, nrow = imax, ncol = R)
  d_all <- matrix(NA, nrow = imax, ncol = P)
  pi_all <- matrix(NA, nrow = imax, ncol = R)
  class_all <- matrix(NA, nrow = imax, ncol = N)
  
  if(plot)
    plot_EMC2(X1, y1, X, y, mus, covMats,0)
  
  ##### EMC2 algorithm
  while(i < imax & crit > eps){
    d_cur <- d_next
    
    ##### EXPECTATION STEP
    for(p in 1:N){
      for(k in 1:R){
        post[p,k] <- pi[k] * mixtools::dmvnorm(X[p,], mus[[k]] + alpha[k]*d_cur, covMats[[k]])
      }
      sig_tmp <- covMats
      # while(sum(post[p,])==0){ # increase ellipse x2
      #   print('Increase')
      #   sig_tmp <- lapply(sig_tmp, function(x) x%*%diag(rep(2,nrow(sig_tmp[[k]]))))
      #   for(k in 1:R){
      #     post[p,k] <- mixtools::dmvnorm(X[p,], mus[[k]] + alpha[k]*d_cur, sig_tmp[[k]])
      #   }
      # }
      post[p,] <- post[p,] / sum(post[p,])
    }
    
    if(plot){
      ##### APPLY THE CORRECTION
      Xc <- X
      class <- apply(post, 1, which.max)
      for(p in 1:nrow(X)){
        Xc[p,] <- Xc[p,] - alpha[class[p]] * d_next
      }
      plot_EMC2(X1, y1, Xc, y, mus, covMats,i)
    }
    
    ##### MAXIMIZATION STEP
    # save before updating
    d_all[i,] <- d_next
    alpha_all[i,] <- alpha
    pi_all[i,] <- pi
    class_all[i,] <- apply(post, 1, which.max)
    
    # Compute pi
    pi <- colMeans(post)
    
    # Compute drift direction
    d1 <- d2 <- 0 
    for(n in 1:N){
      for(k in 1:R){
        d1 <- d1 + post[n,k] * alpha[k]**2 * (covMats_inv[[k]])
        d2 <- d2 + post[n,k] * alpha[k] * (covMats_inv[[k]]) %*% matrix(X[n,] - mus[[k]])
      }
    }
    d_next <- MASS::ginv(d1) %*% d2
    d_next <- d_next / sqrt(sum(d_next**2))
    
    # Compute alpha
    for(k in 1:R){
      if(pi[k] > pi_min){
        alpha[k] <- (t(d_next) %*% covMats_inv[[k]] %*% (t(X) %*% matrix(post[,k]) - matrix(mus[[k]]) * sum(post[,k]))) / as.numeric(t(d_next) %*% covMats_inv[[k]] %*% d_next * sum(post[,k]))
      }
    }
    
    ##### CONVERGENCE ?
    crit <- sqrt(sum((d_next-d_cur)**2))
    i <- i+1
  }
  
  ##### APPLY THE CORRECTION
  Xc <- X
  class <- apply(post, 1, which.max)
  for(p in 1:nrow(X)){
    Xc[p,] <- Xc[p,] - alpha[class[p]] * d_next
  }
  
  ##### SAVE EM OBJECT
  em <- list(post = post, d = d_next, alpha = alpha, pi = pi)
  em_all <- list(d = d_all, alpha = alpha_all, pi = pi_all, class = class_all)
  
  return(list(Xc = Xc, class = class, em = em, em_all = em_all))
}
