saemixDir<-"/home/eco/work/saemix/saemixextension"
versionDir<-"/home/eco/work/saemix/versions"
progDir<-file.path(saemixDir,"R")
datDir<-file.path(saemixDir,"data")
figDir<-file.path(saemixDir,"documentation","figs")

# Libraries
library(ggplot2)
library(MASS)
library(rlang)
library(gridExtra)
library(tidyverse)
library(devtools)

source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"main.R"))
source(file.path(progDir,"func_aux.R"))
source(file.path(progDir,"main_initialiseMainAlgo.R"))
source(file.path(progDir,"main_estep.R"))
source(file.path(progDir,"main_mstep.R"))
source(file.path(progDir,"func_FIM.R"))
source(file.path(progDir,"func_npde.R"))
source(file.path(progDir,"func_plots.R"))
source(file.path(progDir,"func_distcond.R"))
source(file.path(progDir,"func_simulations.R"))
source(file.path(progDir,"compute_LL.R"))
source(file.path(progDir,"func_estimParam.R"))
source(file.path(progDir,"backward.R"))
source(file.path(progDir,"forward.R"))
source(file.path(progDir,"stepwise.R"))
source(file.path(progDir,"func_stepwise.R"))
source(file.path(progDir,"func_compare.R"))

# Data
toenail.saemix<-read.table(file.path(datDir, "toenail.saemix.tab"), header=TRUE)

saemix.data<-saemixData(name.data=toenail.saemix,name.group=c("id"),name.predictors=c("time","y"), name.response="y",
                        name.covariates=c("treatment"),name.X=c("time"))

# Model
binary.model<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  logpdf<-rep(0,length(tim))
  P.obs = (y==0)*(1-pevent)+(y==1)*pevent
  logpdf <- log(P.obs)
  return(logpdf)
}

saemix.model<-saemixModel(model=binary.model,description="Binary model",
                          modeltype="likelihood",
                          psi0=matrix(c(-0.5,-.15,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("alpha","beta"))),
                          transform.par=c(0,0), covariate.model=c(0,1),
                          covariance.model=matrix(c(1,0,0,1),ncol=2), omega.init=diag(c(0.5,0.3)))

saemix.model.iiv1<-saemixModel(model=binary.model,description="Binary model",
                          modeltype="likelihood",
                          psi0=matrix(c(-0.5,-.15,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("alpha","beta"))),
                          transform.par=c(0,0), covariate.model=c(0,1),
                          covariance.model=matrix(c(1,0,0,0),ncol=2), omega.init=diag(c(0.5,0.3)))

saemix.model.nocov<-saemixModel(model=binary.model,description="Binary model",
                          modeltype="likelihood",
                          psi0=matrix(c(-0.5,-.15,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("alpha","beta"))),
                          transform.par=c(0,0), covariate.model=c(0,0),
                          covariance.model=matrix(c(1,0,0,1),ncol=2), omega.init=diag(c(0.5,0.3)))

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE)
binary.fit<-saemix(saemix.model,saemix.data,saemix.options)

plot(binary.fit, plot.type="convergence")

yfit.nocov<-saemix(saemix.model.nocov,saemix.data,saemix.options)
yfit.iiv1<-saemix(saemix.model.iiv1,saemix.data,saemix.options)


yfit<- yfit.iiv1

simulBinary<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-1/(1+exp(-logit))
  print(summary(logit))
  ysim<-rbinom(length(tim),size=1, prob=pevent)
  return(ysim)
}

nsim<-100
yfit <- simulateDiscreteSaemix(yfit, simulBinary, nsim=nsim)
simdat <-yfit@sim.data@datasim
simdat$visit<-rep(toenail.saemix$visit,nsim)
simdat$treatment<-rep(toenail.saemix$treatment,nsim)

#########################################################################################
saemixDir <- "/home/eco/work/saemix/saemixextension"

# Libraries
library(saemix)

# Bootstrap code
source(file.path(saemixDir,"bootstrap","saemix_bootstrap.R"))

#########################################################################################
# Fitting toenail data in saemix
data(toenail.saemix)

saemix.data<-saemixData(name.data=toenail.saemix,name.group=c("id"),name.predictors=c("time","y"), name.response="y",
                        name.covariates=c("treatment"),name.X=c("time"))

binary.model<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  logpdf<-rep(0,length(tim))
  P.obs = (y==0)*(1-pevent)+(y==1)*pevent
  logpdf <- log(P.obs)
  return(logpdf)
}

saemix.model<-saemixModel(model=binary.model,description="Binary model",
                          modeltype="likelihood",
                          psi0=matrix(c(-0.5,-.15,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("alpha","beta"))),
                          transform.par=c(0,0), covariate.model=c(0,1),covariance.model=matrix(c(1,0,0,0),ncol=2), omega.init=diag(c(0.5,0.3)))

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE)
binary.fit<-saemix(saemix.model,saemix.data,saemix.options)

plot(binary.fit, plot.type="convergence")
binary.fit <- conddist.saemix(binary.fit, nsamp=100)

#
saemixObject<-binary.fit
nsamp<-100
max.iter<-NULL
displayPlot<-FALSE


