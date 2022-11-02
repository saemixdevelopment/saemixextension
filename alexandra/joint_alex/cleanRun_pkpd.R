#################################################
saemixDir <- "C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension"
workDir <- file.path(saemixDir, "alexandra","ecoJoint")
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
data_pkpd.prop <- read.csv(file.path(saemixDir, "alexandra","jointTTE", "simAlex_pkpd_proportional.csv"), header=TRUE)
dataPKPD<-saemixData(name.data=data_pkpd.prop, name.group=c("id"), name.predictors=c("time","dose"), 
                     name.response="obs",name.ytype = "ytype")
# Simulation parameters
param<-c(2,8,2,100,5)
omega.sim<-c(1.4, 0.15, 0.3, 0.15, 0.2)
sigma.sim <- 0.2

# model
pkpd<-function(psi,id,xidep) {
  ka <- psi[id,1] 
  V <- psi[id,2]
  Cl <- psi[id,3]
  Emax <- psi[id,4]
  EC50 <- psi[id,5]
  
  tim<-xidep[,1] 
  D <- xidep[,2]
  ytype<-xidep$ytype
  
  ypred <- (D/V)*(ka/(ka-Cl/V))*(exp(-(Cl/V)*tim)-exp(-ka*tim))
  ypred[ytype==2] <- Emax[ytype==2]*ypred[ytype==2]/(ypred[ytype==2]+EC50[ytype==2])
  #  ypd <- Emax*ypred/(ypred+EC50)
  #  ypred[ytype==2] <- ypd[ytype==2]
  return(ypred)
}

# Proportional error model
pkpdmodel.prop<-saemixModel(model=pkpd,description="joint pkpd",modeltype=rep("structural",2),
                            psi0=matrix(param,ncol=5,byrow=TRUE,dimnames=list(NULL, c("ka","V","Cl","Emax","EC50"))),
                            transform.par=c(1,1,1,1,1), covariance.model=diag(c(1,1,1,1,1)),fixed.estim = c(1,1,1,1,1),
                            omega.init = diag(rep(0.5,5)),error.model = c("proportional","proportional"),error.init = c(0,1,0,1),
                            name.sigma = c("a.1","b.1","a.2","b.2"))

# Testing and comparing to data - added ytype to second equation
xidep1<-dataPKPD@data[,c(3,4,9)]
id1<-dataPKPD@data$index
psi1<-do.call(rbind,rep(list(param),dataPKPD@N))
fpred1 <- pkpdmodel.prop@model(psi1, id1, xidep1)
par(mfrow=c(1,2))
for(yt in 1:2) {
  plot(fpred1[xidep1$ytype==yt],dataPKPD@data$obs[xidep1$ytype==yt], xlab="Population predictions", ylab="Observations", main=ifelse(yt==1,"PK","PD"), pch=20)
  abline(0,1)
}

# need to adjust by hand for the moment
pkpdmodel.prop@name.sigma <-c(pkpdmodel.prop@name.sigma,"a.2","b.2")

################################################# Running
# Computational function
source(file.path(workDir,"multi_aux.R"))
source(file.path(workDir,"multi_initializeMainAlgo.R"))
source(file.path(workDir,"multi_estep.R"))
source(file.path(workDir,"multi_mstep.R"))
source(file.path(workDir,"multi_main.R"))
saemix.data<-dataPKPD
saemix.model<-pkpdmodel.prop
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
