#################################################
saemixDir <- "C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension"
workDir <- file.path(saemixDir, "alexandra","joint_alex")
setwd(workDir)

library(ggplot2)
library(Cairo)
library("viridis")
library(rlang)

# Chargement des fonctions originelles de la librairie
progDir<-file.path(saemixDir, "R")
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"func_plots.R")) # for saemix.plot.setoptions

################################################# Data and model (original files)
# Creating data and model objects

# data
data_joint <- read.csv("C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/joint_comp.csv", header=TRUE)
dataJM<-saemixData(name.data=data_joint, name.group=c("id"), name.predictors=c("time"),
                   name.response="obs",name.ytype = "ytype")

#model JM lin + 2 competing risks
JMmodel<-function(psi,id,xidep) {
  ytype<-xidep$ytype  
  
  p1 <- psi[id,1]  ## baseline risks 
  p2 <- psi[id,2]
  g1 <- psi[id,3]  ## effet cov sur risques 1 et 2 
  g2 <- psi[id,4]
  alpha1 <- psi[id,5]
  alpha2 <- psi[id,6]
  b0 <- psi[id,7]
  b1 <- psi[id,8]
  
  T1<-xidep[ytype==2,1] # vector of times partie survie ev 1
  Nj1 <- length(T1)
  cens1<-which(T1==max(T1))  # censoring time=30
  init1 <- which(T1==0)
  ind1 <- setdiff(1:Nj1, append(init1,cens1)) # indices of events
  
  T2<-xidep[ytype==3,1] # vector of times partie survie ev 2
  Nj2 <- length(T2)
  cens2<-which(T2==max(T2))  # censoring time=30
  init2 <- which(T2==0)
  ind2 <- setdiff(1:Nj2, append(init2,cens2))
  
  p1b <- unique(p1)
  p2b <- unique(p2)
  g1b <- unique(g1)  
  g2b <- unique(g2)
  alpha1b <- unique(alpha1)
  alpha2b <- unique(alpha2)
  b0b <- unique(b0)
  b1b <- unique(b1)
  
  f=function(x) seq(0,x,length.out=1000)
  tab1 = mapply(f,T1[T1!=0])
  tab1 = t(tab1)
  pas1 = tab1[,2]-tab1[,1]
  
  haz1 = p1b*g1b*exp(-g1b*tab1)/(1-p1b*(1-exp(-g1b*tab1)))*exp(alpha1b*(b0b+b1b*tab1))
  H1 = apply(haz1,1,sum)*pas1
  hazt1 = haz1[,1000]
  
  tab2 = mapply(f,T2[T2!=0])
  tab2 = t(tab2)
  pas2 = tab2[,2]-tab2[,1]
  
  haz2 = p2b*g2b*exp(-g2b*tab2)/(1-p2b*(1-exp(-g2b*tab2)))*exp(alpha2b*(b0b+b1b*tab2))
  H2 = apply(haz2,1,sum)*pas2
  hazt2 = haz2[,1000]
  
  logpdf1 <- rep(0,Nj1)
  logpdf1[cens1] <- -H1[cens1/2] #+ H1[cens1-1] à généraliser par la suite (facile)
  logpdf1[ind1] <- -H1[ind1/2] + log(hazt1[ind1/2]) #+ H1[ind1-1]
  
  logpdf2 <- rep(0,Nj2)
  logpdf2[cens2] <- -H2[cens2/2] #+ H2[cens2-1]
  logpdf2[ind2] <- -H2[ind2/2] + log(hazt2[ind2/2]) #+ H2[ind2-1]
  
  ypred = b0+b1*xidep[,1]
  
  ypred[ytype==2] = logpdf1
  ypred[ytype==3] = logpdf2
  
  return(ypred)
}


# joint TTE  
param<-c(0.01,0.02,0.2,0.1,0,0,15,0.1) # p1, g1, p2, g2, alpha1, alpha2, b0, b1 
jointTTE<-saemixModel(model=JMmodel,description="JM lin+competing risks",modeltype=c("structural","likelihood","likelihood"),
                      psi0=matrix(param,ncol=8,byrow=TRUE,dimnames=list(NULL, c("p1", "g1", "p2", "g2", "alpha1", "alpha2", "b0", "b1"))),
                      transform.par=c(1,1,1,1,0,0,0,0), covariance.model=diag(c(0,0,0,0,0,0,1,1)),
                      fixed.estim = c(1,1,1,1,1,1,1,1), omega.init = diag(c(0,0,0,0,0,0,1,0.1)))



################################################# Running
# Computational function
source(file.path(workDir,"multi_aux2.R"))
source(file.path(workDir,"multi_initializeMainAlgo.R"))
source(file.path(workDir,"multi_estep.R"))
source(file.path(workDir,"multi_mstep.R"))
source(file.path(workDir,"multi_main.R"))
source(file.path(workDir,"multi_map.R"))
saemix.data<-dataJM
saemix.model<-jointTTE
saemix.options<-saemixControl(seed=12345, map=FALSE, fim=FALSE, ll.is=FALSE)

yfit <- saemix.multi(saemix.model, saemix.data, saemix.options)

# Population estimates
yfit
param
sigma.sim
sqrt(diag(yfit@results@omega))
omega.sim

# Individual estimates
yfit1 <- map.saemix(yfit)
summary(yfit1@results@map.psi)

ipar<-data_pkpd.prop[!duplicated(data_pkpd.prop$id),6:10]
plot(yfit1@results@map.psi)
par(mfrow=c(2,3))
for(i in 1:5) {
  plot(ipar[,i], yfit1@results@map.psi[,i], xlab=paste("Simulated",colnames(ipar)[i]), ylab=paste("Estimated",colnames(ipar)[i]), pch=20)
  abline(0,1)
}
################################################# NOT DONE (need to be adjusted to multiple response models)
# compute FIM (and LL by linearisation)
# compute LL by IS
# conditional distributions
# plots
# all the other functions

# Fonctions originelles (non ajustées)
if(FALSE) {
  source(file.path(progDir,"func_FIM.R"))
  source(file.path(progDir,"func_distcond.R"))
  source(file.path(progDir,"func_simulations.R"))
  source(file.path(progDir,"compute_LL.R"))
  source(file.path(progDir,"func_estimParam.R"))
  source(file.path(progDir,"func_npde.R"))
  source(file.path(progDir,"backward.R"))
  source(file.path(progDir,"forward.R"))
  source(file.path(progDir,"stepwise.R"))
  source(file.path(progDir,"func_stepwise.R"))
  source(file.path(progDir,"func_compare.R"))
  source(file.path(progDir,"func_bootstrap.R"))
  source(file.path(progDir,"func_exploreData.R"))
  source(file.path(progDir,"func_discreteVPC.R"))
  
}
