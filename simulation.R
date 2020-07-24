library(mixtools)
library(tidyr)
library(e1071)
library(MASS)
source('PCACC.R')
source('EMC2.R')
source('OSC.R')
source('CPCACC.R')
source('EvolutionaryMethod.R')

P = 4 # number of chemical sensors
R = 6 # number of compounds
N1 = 60 # number of samples in Session 1 (in total)

deltat = 8 # window size
ts = 0:251 # experiment duration
beta = 12 # hyperparemeter
h = 30 # hyperparemeter
Ntimes = 10 # number of simulations
# Ntimes = 100
cr_raw  = cr_emc2 = cr_diCarlo  = cr_pcacc = matrix(0,nrow=Ntimes,ncol=R-1) # classification rates
Nwindow = length(ts)/(deltat+1) # number of windows (sessions)
m = 1
while(m<=Ntimes){
  print(paste("SIMULATION",m))
  centroids = rmvnorm(R,rep(0,P), diag(rep(1,P)) * beta**2/(2*P))
  drift_direction = runif(P,-1,1)
  drift_direction = drift_direction / sqrt(sum(drift_direction**2))
  alpha = runif(R,0.5,3)
  
  # check the separation
  dist_centroid = dist(centroids)
  if(any(dist_centroid < beta/2)){
    print("Simulation discarded")
  }else{

    for(K in 2:R){
      print(K)
      molToKeep = 1:K
      
      #### GENERATE SESSION 1
      X1 <- NULL
      for(r in 1:R){
        X1 = rbind(X1,rmvnorm(N1/R, centroids[r,], diag(rep(1,P))))
      }
      y1 <- rep(1:R, each = N1/R)
      
      #### GENERATE DRIFTING SESSIONS
      Xcont <- NULL
      for(t in ts){
        for(r in molToKeep){
          Xcont = rbind(Xcont,rmvnorm(1, centroids[r,] + ((t*drift_direction)/(h) + t/(h)*sin(t/(h)*drift_direction)) * alpha[r]  , diag(rep(1,P))))
        }
      }
      ycont <- rep(molToKeep, length(ts))
      tall <- rep(ts, each = length(molToKeep))
      
      #### DRIFT CORRECTION WINDOW-BY-WINDOW
      Ws = ts %>% matrix(ncol=deltat+1, byrow=T)
      Xcont_cor <- Xcont_pcacc <- Xcont
      cr_cont = cr_cont_diCarlo = numeric(nrow(Ws))
      Minit = NULL
      for(nw in 1:nrow(Ws)){
        # print(paste('Di Carlo ', nw))
        tmp = which(tall%in%Ws[nw,])

        # EMC2
        res_emc2 = EMC2(X1, y1, Xcont[tmp,])
        Xcont_cor[tmp,] = res_emc2$Xc
        cr_cont[nw] = mean(ycont[tmp] == res_emc2$class)
        
        # # Di Carlo method
        res_diCarlo = evolutionaryBased_DiCarlo(X1, y1, Xcont[tmp,], Minit=Minit, optimMethod='optim')
        cr_cont_diCarlo[nw] =  mean(ycont[tmp] == res_diCarlo$pred)
        Minit = res_diCarlo$M
      }
      
      #### DRIFT CORRECTION PCA-CC
      Xcal <- rbind(X1[y1==1,],Xcont[ycont==1,])
      pcacc <- PCACC_tune(Xcal, ncomp=NULL, varExp=0.8)
      Xcont_pcacc = PCACC_cor(Xcont, pcacc)
      
      # Raw method
      svm_tune <- tune(svm, train.x=X1, train.y=factor(y1), kernel='linear', scale = F, ranges=list(cost=c(10^(-5:5))))
      predRaw <- predict(svm_tune$best.model, Xcont) %>% as.numeric
      cr_raw[m,K-1] = mean(predRaw == ycont)
      
      predPCACC <- predict(svm_tune$best.model, Xcont_pcacc) %>% as.numeric
      cr_pcacc[m,K-1] = mean(predPCACC == ycont)
      
      cr_emc2[m,K-1] = mean(cr_cont)
      cr_diCarlo[m,K-1] = mean(cr_cont_diCarlo)
    }
    
    # save results just in case
    # saveRDS(cr_raw, "classRate_raw.rds") 
    # saveRDS(cr_emc2, "classRate_EMC2.rds")
    # saveRDS(cr_diCarlo, "classRate_diCarlo.rds")
    # saveRDS(cr_pcacc, "classRate_PCACC.rds")
    m = m + 1
  }
}
