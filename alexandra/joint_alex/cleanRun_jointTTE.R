#################################################
saemixDir <- "C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension"
workDir <- file.path(saemixDir, "alexandra","joint_alex")
setwd(workDir)

library(ggplot2)
library(Cairo)
library("viridis")  

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
data_joint <- read.csv("C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/joint_tte.csv", header=TRUE)
dataJM<-saemixData(name.data=data_joint, name.group=c("id"), name.predictors=c("time"), 
                   name.response="obs",name.ytype = "ytype")

#model JM longi+TTE
JMmodel<-function(psi,id,xidep) {
  ytype<-xidep$ytype  # type of response (1: continuous, 2: event)
  b0 <- psi[id,1]  ## coeffs
  b1 <- psi[id,2]  ## longi
  h0 <- psi[id,3]
  alpha <- psi[id,4]  ## coeff de lien
  
  ypred <- b0+b1*xidep[,1]  ## pred longi 
  
  #obs <- xidep[ytype==2,1] # vector of observations partie survie
  T<-xidep[ytype==2,1] # vector of times partie survie
  Nj <- length(T)
  cens<-which(T==max(T))  # censoring time=30
  init <- which(T==0)
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  b0b = b0[ytype==2]
  b1b = b1[ytype==2]
  h0b = h0[ytype==2]
  alphab = alpha[ytype==2]
  
  haz <- h0b*exp(alphab*(b0b+b1b*T))
  H <- (h0b/(alphab*b1b))*exp((b0b+b1b*T)*alphab)-(h0b/(alphab*b1b))*exp(alphab*b0b)
  
  logpdf <- rep(0,Nj)
  logpdf[cens] <- -H[cens] + H[cens-1]
  logpdf[ind] <- -H[ind] + H[ind-1] + log(haz[ind])
  
  ypred[ytype==2] = logpdf
  return(ypred)
}


# joint TTE  
param<-c(15,0.3,0.01,0.1)
jointTTE<-saemixModel(model=JMmodel,description="JM lin longi one tte",modeltype=c("structural","likelihood"),
                      psi0=matrix(param,ncol=4,byrow=TRUE,dimnames=list(NULL, c("b0","b1","h0","alpha"))),
                      transform.par=c(0,0,1,0), covariance.model=diag(c(1,1,0,0)),
                      fixed.estim = c(1,1,1,1),error.model = "constant",
                      omega.init = diag(c(0.25,0.01,0,0)))



################################################# Running
# Computational function
source(file.path(workDir,"multi_aux2.R"))
source(file.path(workDir,"multi_initializeMainAlgo.R"))
source(file.path(workDir,"multi_estep.R"))
source(file.path(workDir,"multi_mstep.R"))
source(file.path(workDir,"multi_main.R"))
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
