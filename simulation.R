# Author: Pierre Maho
# Lab: Gipsa-Lab
# PhD Supervisors: Simon Barthelme and Pierre Comon
# Date: 2020/01/12


suppressMessages(library(ggplot2))
suppressMessages(library(extrafont))
suppressMessages(library(MASS))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(e1071))
suppressMessages(library(grid))
suppressMessages(library(ggthemes))
suppressMessages(loadfonts())
source('PCACC.R')
source('PLSCC.R')
source('OSC.R')
source('CPCACC.R')
source('EvolutionaryMethod.R')

# A function for nice plots
theme_Publication <- function(base_size=10, base_family="LM Roman 10") {
  # taken from https://rstudio-pubs-static.s3.amazonaws.com/71792_1acccaaf1b5b4fac88cb986da852c2df.html
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = "black", fill=NA),
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(1, "cm"),
            legend.margin = unit(0.5, "cm"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="black",fill="#f0f0f0")))
}

##### GENERATE ARTIFICIAL DATA SET
# parameters
P  <- 5# number of chemical sensors
R  <- 5 # number of VOCs
N  <- 30 # number of repetitions/session
Ns <- 5 # number of session


# Generate Artificial data set
# Drift is a single translation vector
# Each class is drifting along this direction but with a different weight
mus <- replicate(R,runif(P,-15,15)) # class centroids
Sig <- diag(rep(1,P)) # covariance matrix

# Generate Artificial data set
drift_direction  <- runif(P, 0, 20)
drift_class      <- replicate(R,drift_direction)
weight_class     <-  replicate(R,runif(Ns, 0, 1)) # coefficient alpha in PhD thesis (class weight of the drift)
weight_class[1,] <- 0 # first session is not drifting

Xs <- NULL # all session data
for(i in 1:Ns){
  Xi <- NULL
  for(j in 1:R){
    # drift_class[,j] <- drift_class[,j] + drift_direction * weight_class[i,j] # cumulate the drift
    drift_class[,j] <- drift_class[,j] + drift_direction * 1 # cumulate the drift
    Xi              <- rbind(Xi,mvrnorm(N, mus[,j] + drift_class[,j], Sig))  
  }
  Xs <- rbind(Xs,Xi)
}
mols     <- rep(rep(1:R, each = N), Ns)
sessions <-  rep(1:Ns, each = N*R)

# Plot the artificial data set
pca_session1 <- prcomp(Xs[sessions==1,])
df           <- data.frame(Xs, mol = mols, session = sessions) %>% cbind(Xs%*%pca_session1$rotation[,1:5])
ggplot(df,aes(PC1, PC2,col=as.factor(mol),shape=as.factor(session))) + 
  geom_point() + 
  xlab("Principal direction 1 of Session 1")+
  ylab("Principal direction 2 of Session 1")+
  ggtitle("Raw data") +
  scale_color_discrete("VOC")+
  scale_shape_manual("Session",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)


#### CALIBRANT SCENARIO
# Correct with PCA-CC
Xcal         <- Xs[mols==1,] # calibrant matrix
pcacc        <- PCACC_tune(Xcal, ncomp=1) 
Xcor         <- PCACC_cor(Xs, pcacc)
pca_session1 <- prcomp(Xcor[sessions==1,])
df           <- data.frame(Xcor, mol = mols, session = sessions) %>% cbind(Xcor%*%pca_session1$rotation[,1:5])
ggplot(df,aes(PC1, PC2,col=as.factor(mol),shape=as.factor(session))) + 
  geom_point() + 
  xlab("Principal direction 1 of Session 1")+
  ylab("Principal direction 2 of Session 1")+
  ggtitle("PCA-CC") +
  scale_color_discrete("VOC")+
  scale_shape_manual("Session",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)

# Correct with PLS-CC
Xcal         <- Xs[mols==1,] # calibrant matrix
ycal         <- 1:nrow(Xcal) 
plscc        <- PLSCC_tune(Xcal, ycal, ncomp=1) 
Xcor         <- PLSCC_cor(Xs, plscc)
pca_session1 <- prcomp(Xcor[sessions==1,])
df           <- data.frame(Xcor, mol = mols, session = sessions) %>% cbind(Xcor%*%pca_session1$rotation[,1:5])
ggplot(df,aes(PC1, PC2,col=as.factor(mol),shape=as.factor(session))) + 
  geom_point() + 
  xlab("Principal direction 1 of Session 1")+
  ylab("Principal direction 2 of Session 1")+
  ggtitle("PLS-CC") +
  scale_color_discrete("VOC")+
  scale_shape_manual("Session",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)


#### MULTI SESSION SCENARIO
# Correct with OSC
Xtrain       <- Xs[sessions %in% 1:2,]
ytrain       <- mols[sessions %in% 1:2]
# osc          <- OSC_tune(Xtrain,ytrain,ncomp=1)
# Xcor         <- OSC_cor(Xs,osc)
osc          <- OSC_DO_tune(Xtrain,ytrain,ncomp=1)
Xcor         <- OSC_DO_cor(Xs,osc)
pca_session1 <- prcomp(Xcor[sessions==1,])
df           <- data.frame(Xcor, mol = mols, session = sessions) %>% cbind(Xcor%*%pca_session1$rotation[,1:5])
ggplot(df,aes(PC1, PC2,col=as.factor(mol),shape=as.factor(session))) + 
  geom_point() + 
  xlab("Principal direction 1 of Session 1")+
  ylab("Principal direction 2 of Session 1")+
  ggtitle("OSC") +
  scale_color_discrete("VOC")+
  scale_shape_manual("Session",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)

# Correct with CPCA
Xtrain       <- Xs[sessions %in% 1:2,]
ytrain       <- mols[sessions %in% 1:2]
cpca         <- CPCA_tune(Xtrain,ytrain,ncomp=1)
Xcor         <- CPCA_cor(Xs,cpca)
pca_session1 <- prcomp(Xcor[sessions==1,])
df           <- data.frame(Xcor, mol = mols, session = sessions) %>% cbind(Xcor%*%pca_session1$rotation[,1:5])
ggplot(df,aes(PC1, PC2,col=as.factor(mol),shape=as.factor(session))) + 
  geom_point() + 
  xlab("Principal direction 1 of Session 1")+
  ylab("Principal direction 2 of Session 1")+
  ggtitle("CPCA") +
  scale_color_discrete("VOC")+
  scale_shape_manual("Session",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)




#### BLIND SCENARIO
# Correct with the method of Di Carlo et al (CMA-ES)
# The method is a bit long so just test on Session 1 and 2, and with a different kind of drift (appropriate for the drift correction method
# We do not use CMA-ES as proposed by the authors. This method is quite long while BFGS method is much faster and gives similar results. However, we can try CMA-ES bty setting optimMethod to 'cmaes'
Xtrain <- Xs[sessions == 1,]
ytrain <- mols[sessions == 1]
Mgt <- runif(P*P,0,0.1) %>% matrix(ncol = ncol(Xs), nrow = ncol(Xs))
diag(Mgt) <- 1
Xtest <- Xtrain %*% Mgt

Xnew <- rbind(Xtrain, Xtest)
pca_session1 <- prcomp(Xnew[which(sessions==1),])
df           <- data.frame(Xnew, mol = mols[sessions%in%1:2], session = sessions[sessions%in%1:2]) %>%
  cbind(Xnew%*%pca_session1$rotation[,1:5])
ggplot(df,aes(PC1, PC2,col=as.factor(mol),shape=as.factor(session))) + 
  geom_point() + 
  xlab("Principal direction 1 of Session 1")+
  ylab("Principal direction 2 of Session 1")+
  ggtitle("New raw data (for Di Carlo's method)") +
  scale_color_discrete("VOC")+
  scale_shape_manual("Session",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)

res_diCarlo <- evolutionaryBased_DiCarlo(Xtrain, ytrain, Xtest, optimMethod='optim')
Xnew_cor <- rbind(Xtrain, res_diCarlo$X2cor)
pca_session1 <- prcomp(Xnew_cor[which(sessions==1),])
df           <- data.frame(Xnew_cor, mol = mols[sessions%in%1:2], session = sessions[sessions%in%1:2]) %>%
  cbind(Xnew_cor%*%pca_session1$rotation[,1:5])
ggplot(df,aes(PC1, PC2,col=as.factor(mol),shape=as.factor(session))) + 
  geom_point() + 
  xlab("Principal direction 1 of Session 1")+
  ylab("Principal direction 2 of Session 1")+
  ggtitle("Di Carlo's method") +
  scale_color_discrete("VOC")+
  scale_shape_manual("Session",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)


