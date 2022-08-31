# @Eco
workDir<-"/home/eco/work/saemix/saemixextension/bootstrap"
saemixDir <- "/home/eco/work/saemix/saemixextension"
setwd(workDir)

#  Code saemix
progDir <- file.path(saemixDir,"R")
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

# Bootstrap code
source(file.path(saemixDir, "bootstrap", "saemix_bootstrap.R"))

# Number of bootstrap samples
nboot <- 500

# Data
datDir <- file.path(saemixDir, "data")
toenail.saemix <- read.table(file.path(datDir, "toenail.saemix.tab"), header = T)

saemix.data<-saemixData(name.data=toenail.saemix,name.group=c("id"),name.predictors=c("time","y"), name.response="y",
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

saemix.model<-saemixModel(model=binary.model, description="Binary model", simulate.function = simulBinary,
                          modeltype="likelihood",
                          psi0=matrix(c(-0.5,-.15,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                          transform.par=c(0,0), covariate.model=c(0,1),covariance.model=matrix(c(1,0,0,0),ncol=2), omega.init=diag(c(0.5,0.3)))
saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE)
binary.fit<-saemix(saemix.model,saemix.data,saemix.options)


# Case bootstrap
case.bin <- saemix.bootstrap(binary.fit, method="case", nboot=nboot) 
write.table(case.bin, file.path(workDir, "results", "toenail_caseBootstrap.res"), quote=F, row.names=FALSE)

# Conditional non-parametric bootstrap
if(FALSE) {
  cond.bin <- saemix.bootstrap(binary.fit, method="conditional", nboot=nboot) 
  write.table(cond.bin, file.path(workDir, "results", "toenail_condBootstrap.res"), quote=F, row.names=FALSE)
}
