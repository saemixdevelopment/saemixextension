# Folders
saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
datDir<-file.path(saemixDir,"data")

# Libraries
library(MASS)
library(rlang)
library(ggplot2)
library(gridExtra)

# Sourcing saemix functions
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
  source(file.path(progDir,"func_plots.R"))
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

###################################
# Theophylline dataset
theo.saemix<-read.table(file.path(datDir, "theo.saemix.tab"), header=TRUE)
theo.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
                        name.group=c("Id"),name.predictors=c("Dose","Time"),
                        name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
                        units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}
theo.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",
                          psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE,  dimnames=list(NULL, c("ka","V","CL"))),
                          transform.par=c(1,1,1))
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
theo.fit<-saemix(theo.model,theo.data,saemix.options)

# Without simulated annealing
saemix.options$nbiter.sa<-0
theo.fit.noSA<-saemix(theo.model,theo.data,saemix.options)

###################################
# Binary data
toenail.saemix <- read.table(file.path(datDir, "toenail.saemix.tab"), header=TRUE)
toe.data<-saemixData(name.data=toenail.saemix,name.group=c("id"),name.predictors=c("time","y"), name.response="y",
                        name.covariates=c("treatment"))
# Model function
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
# Simulation function
simulBinary<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-1/(1+exp(-logit))
  ysim<-rbinom(length(tim),size=1, prob=pevent)
  return(ysim)
}
toe.model<-saemixModel(model=binary.model, description="Binary model", simulate.function = simulBinary,
                          modeltype="likelihood",
                          psi0=matrix(c(-0.5,-.15,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                          transform.par=c(0,0), covariate.model=c(0,1),covariance.model=matrix(c(1,0,0,0),ncol=2), omega.init=diag(c(0.5,0.3)))
saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE)
binary.fit<-saemix(toe.model,toe.data,saemix.options)
binary.fit2<-saemix(toe.model,toe.data,saemix.options)

# Without simulated annealing
saemix.options$nbiter.sa<-5
binary.fit.noSA<-saemix(toe.model,toe.data,saemix.options)

