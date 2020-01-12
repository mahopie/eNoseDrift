suppressMessages(library(ggplot2))
suppressMessages(library(extrafont))
suppressMessages(loadfonts())
source('PCACC.R')
source('PLSCC.R')
source('OSC.R')
source('CPCACC.R')
source('EvolutionaryMethod.R')

# A function for nice plots
theme_Publication <- function(base_size=10, base_family="LM Roman 10") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
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
            # legend.text = element_text(size=base_size-5),
            # legend.title = element_text(size=base_size-5),
            legend.margin = unit(0.5, "cm"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="black",fill="#f0f0f0")
    ))
  
}

##### GENERATE ARTIFICIAL DATA SET
# parameters
P <- 5# number of chemical sensors
R <- 3 # number of VOCs
N <- 30 # number of repetitions/session
Ns <- 5 # number of session


# Generate Artificial data set
# Drift is a single translation vector for each session
# Each class is drifting along this direction bu with a different weight
mus <- replicate(R,runif(P,-15,15)) # class centroids
Sigs <- list() # class covariance
for(i in 1:R){
  A <- matrix(runif(P^2,-1,1)*2-1, ncol=P)
  Sigs[[i]] <- t(A) %*% A * 2/P
}

# Generate Artificial data set
drift_direction <- runif(P, 0, 40)
drift_class <- replicate(R,drift_direction)
coef <-  replicate(R,runif(Ns, 0, 10))
coef[1,] <- 0 # first session is not drifting

Xs <- NULL # all session data
for(i in 1:Ns){
  Xi <- NULL
  for(j in 1:R){
    drift_class[,j] <- drift_class[,j] + drift_direction * coef[i,j] # cumulate the drift
    Xi <- rbind(Xi,mvrnorm(N, mus[,j] + drift_class[,j], Sigs[[j]]))  
  }
  Xs <- rbind(Xs,Xi)
}
mol <- rep(rep(1:R, each = N), Ns)
session <- rep(1:Ns, each = N*R)

# Plot the data set
pca <- prcomp(Xs)
df <- as.data.frame(Xs) %>% mutate(mol=mol, session=session,pca1=pca$x[,1], pca2=pca$x[,2],pca3=pca$x[,3])
ggplot(df,aes(pca1, pca2,col=as.factor(mol),shape=as.factor(session))) + geom_point() + xlab("PCA 1")+
  ylab("PCA 2")+
  scale_color_discrete("VOC")+
  scale_shape_manual("Data from",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)


######## DRIFT CORRECTION METHODS ######## 

#### CALIBRANT HYPOTHESIS
# Correct with PCA-CC
Xcal <- Xs[mol==1,]
pcacc <- PCACC_tune(Xcal, ncomp=1) 
Xcor = PCACC_cor(Xs, pcacc)
pca <- prcomp(Xcor)
df <- as.data.frame(Xcor) %>% mutate(mol=mol, session=session,pca1=pca$x[,1], pca2=pca$x[,2],pca3=pca$x[,3])
ggplot(df,aes(pca1, pca2,col=as.factor(mol),shape=as.factor(session))) + geom_point() + xlab("PCA 1")+
  ylab("PCA 2")+
  scale_color_discrete("VOC")+
  scale_shape_manual("Data from",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)

# Correct with PLS-CC
Xcal = Xs[mol==1,]
ycal = 1:nrow(Xcal)
plscc <- PLSCC_tune(Xcal, ycal, ncomp=1) 
Xcor = PLSCC_cor(Xs, plscc)
pca <- prcomp(Xcor)
df <- as.data.frame(Xcor) %>% mutate(mol=mol, session=session,pca1=pca$x[,1], pca2=pca$x[,2],pca3=pca$x[,3])
ggplot(df,aes(pca1, pca2,col=as.factor(mol),shape=as.factor(session))) + geom_point() + xlab("PCA 1")+
  ylab("PCA 2")+
  scale_color_discrete("VOC")+
  scale_shape_manual("Data from",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)


#### MULTI SESSION HYPOTHESIS
# Correct with OSC
Xtrain = Xs[session %in% 1:2,]
ytrain = mol[session %in% 1:2]
Xtest = Xs
osc = OSC_tune(Xtrain,ytrain,ncomp=1)
Xcor = OSC_cor(Xtest,osc)
pca <- prcomp(Xcor)
df <- as.data.frame(Xcor) %>% mutate(mol=mol, session=session,pca1=pca$x[,1], pca2=pca$x[,2],pca3=pca$x[,3])
ggplot(df,aes(pca1, pca2,col=as.factor(mol),shape=as.factor(session))) + geom_point() + xlab("PCA 1")+
  ylab("PCA 2")+
  scale_color_discrete("VOC")+
  scale_shape_manual("Data from",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)


# Correct with CPCA
Xtrain = Xs[session %in% 1:2,]
Xtest = Xs
ytrain = mol[session %in% 1:2]
cpca = CPCA_tune(Xtrain,ytrain,ncomp=1)
Xcor = CPCA_cor(Xtest,cpca)
pca <- prcomp(Xcor)
df <- as.data.frame(Xcor) %>% mutate(mol=mol, session=session,pca1=pca$x[,1], pca2=pca$x[,2],pca3=pca$x[,3])
ggplot(df,aes(pca1, pca2,col=as.factor(mol),shape=as.factor(session))) + geom_point() + xlab("PCA 1")+
  ylab("PCA 2")+
  scale_color_discrete("VOC")+
  scale_shape_manual("Data from",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)



#### BLIND HYPOTHESIS
# Correct with the method of Di Carlo et al (CMA-ES)
# The method is abit long so just test on Session 1 and 2, and with a different kind of drift
Xtrain = Xs[session == 1,]
ytrain = mol[session == 1]
Mgt = runif(P*P,0,0.25) %>% matrix(ncol = ncol(Xs), nrow = ncol(Xs))
diag(Mgt) <- 1
Xtest = Xtrain %*% Mgt

Xnew = rbind(Xtrain, Xtest)
pca <- prcomp(Xnew)
df <- as.data.frame(Xnew) %>% mutate(mol=mol[session %in%1:2], session=session[session %in%1:2],pca1=pca$x[,1], pca2=pca$x[,2],pca3=pca$x[,3])
ggplot(df,aes(pca1, pca2,col=as.factor(mol),shape=as.factor(session))) + geom_point() + xlab("PCA 1")+
  ylab("PCA 2")+
  scale_color_discrete("VOC")+
  scale_shape_manual("Data from",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)

res_diCarlo = evolutionaryBased_DiCarlo(Xtrain, ytrain, Xtest, optimMethod='optim')
Xnew = rbind(Xtrain, res_diCarlo$X2cor)
pca <- prcomp(Xnew)
df <- as.data.frame(Xnew) %>% mutate(mol=mol[session %in%1:2], session=session[session %in%1:2],pca1=pca$x[,1], pca2=pca$x[,2],pca3=pca$x[,3])
ggplot(df,aes(pca1, pca2,col=as.factor(mol),shape=as.factor(session))) + geom_point() + xlab("PCA 1")+
  ylab("PCA 2")+
  scale_color_discrete("VOC")+
  scale_shape_manual("Data from",values=c(1, 17,4,5,8, 7, 13,16,3,9,15))+
  theme_Publication(12)


