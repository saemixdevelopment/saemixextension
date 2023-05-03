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
data_joint <- read.csv("C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/joint_multilongi_1tte.csv", header=TRUE)
dataJM<-saemixData(name.data=data_joint, name.group=c("id"), name.predictors=c("time"),
                   name.response="obs",name.ytype = "ytype")

#model JM 3 lin + 1 tte
JMmodel<-function(psi,id,xidep) {
  ytype<-xidep$ytype  

  b01 <- psi[id,1]
  b11 <- psi[id,2]
  b02 <- psi[id,3]
  b12 <- psi[id,4]
  b03 <- psi[id,5]
  b13 <- psi[id,6]
  h0 <- psi[id,7]
  alpha1 <- psi[id,8]
  alpha2 <- psi[id,9]
  alpha3 <- psi[id,10]
  
  
  T<-xidep[ytype==4,1] # vector of times partie survie
  Nj <- length(T)
  cens<-which(T==max(T))  # censoring time=30
  init <- which(T==0)
  ind <- setdiff(1:Nj, append(init,cens)) # indices of event
  
  b01b = unique(b01)
  b11b = unique(b11)
  b02b = unique(b02)
  b12b = unique(b12)
  b03b = unique(b03)
  b13b = unique(b13)
  h0b = unique(h0)
  alpha1b = unique(alpha1)
  alpha2b = unique(alpha2)
  alpha3b = unique(alpha3)
  
  f=function(x) seq(0,x,length.out=100)
  tab = mapply(f,T[T!=0])
  tab = t(tab)
  pas = tab[,2]-tab[,1]
  
  haz = h0b*exp(alpha1b*(b01b+b11b*tab)+alpha2b*(b02b+b12b*tab)+alpha3*(b03b+b13b*tab))
  H = apply(haz,1,sum)*pas
  hazt = haz[,100]
  
  logpdf <- rep(0,Nj)
  logpdf[cens] <- -H[cens/2] #+ H1[cens-1] à généraliser par la suite (facile)
  logpdf[ind] <- -H[ind/2] + log(hazt[ind/2]) #+ H1[cens-1]
  
  ypred = rep(NA,length(xidep[,1]))
  
  ypred[ytype==1] = b01[ytype==1]+b11[ytype==1]*xidep[ytype==1,1]
  ypred[ytype==2] = b02[ytype==2]+b12[ytype==2]*xidep[ytype==2,1]
  ypred[ytype==3] = b03[ytype==3]+b13[ytype==3]*xidep[ytype==3,1]

  ypred[ytype==4] = logpdf
  
  return(ypred)
}


# joint TTE  
param<-c(15,0.3,7,-0.1,30,0.8,0.00005,0.1,-0.2,0.15) # p1, g1, p2, g2, alpha1, alpha2, b0, b1 
jointTTE<-saemixModel(model=JMmodel,description="JM 3lin+1tte",modeltype=c("structural","structural","structural","likelihood"),
                      psi0=matrix(param,ncol=10,byrow=TRUE,dimnames=list(NULL, c("b01", "b11", "b02", "b12", "b03", "b13", "h0", "alpha1", "alpha2", "alpha3"))),
                      transform.par=c(0,0,0,0,0,0,1,0,0,0), covariance.model=diag(c(1,1,1,1,1,1,0,0,0,0)),
                      fixed.estim = c(1,1,1,1,1,1,1,1,1,1), omega.init = diag(c(0.25,0.01,1,0.01,16,0.36,0,0,0,0)))



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
