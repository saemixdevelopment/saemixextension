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
set.seed(42919)

# Number of bootstrap samples
nboot <- 500

# Data
datDir <- file.path(saemixDir, "data")
knee.saemix <- read.table(file.path(datDir, "knee.saemix.tab"), header = T)

saemix.data<-saemixData(name.data=knee.saemix,name.group=c("id"),
                        name.predictors=c("y", "time"), name.X=c("time"),
                        name.covariates = c("Age","Sex","treatment","Age2"),
                        units=list(x="d",y="", covariates=c("yr","-","-","yr2")))

# Model
ordinal.model<-function(psi,id,xidep) {
  y<-xidep[,1]
  time<-xidep[,2]
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  alp4<-psi[id,4]
  beta<-psi[id,5]
  
  logit1<-alp1 + beta*time
  logit2<-logit1+alp2
  logit3<-logit2+alp3
  logit4<-logit3+alp4
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  pge4<-exp(logit4)/(1+exp(logit4))
  pobs = (y==1)*pge1+(y==2)*(pge2 - pge1)+(y==3)*(pge3 - pge2)+(y==4)*(pge4 - pge3)+(y==5)*(1 - pge4)
  logpdf <- log(pobs)
  
  return(logpdf)
}

# Simulation function
simulateOrdinal<-function(psi,id,xidep) {
  y<-xidep[,1]
  time<-xidep[,2]
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  alp4<-psi[id,4]
  beta<-psi[id,5]
  
  logit1<-alp1 + beta*time
  logit2<-logit1+alp2
  logit3<-logit2+alp3
  logit4<-logit3+alp4
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  pge4<-exp(logit4)/(1+exp(logit4))
  x<-runif(length(time))
  ysim<-1+as.integer(x>pge1)+as.integer(x>pge2)+as.integer(x>pge3)+as.integer(x>pge4)
  return(ysim)
}


# Fitting
saemix.options<-list(seed=632545, save=FALSE, fim=FALSE, save.graphs=FALSE, nb.chains=10, nbiter.saemix=c(600,100))
#saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, nb.chains=10, fim=FALSE)

if(FALSE) {
# Saemix model
saemix.model<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                          simulate.function = simulateOrdinal,
                          psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                          transform.par=c(0,1,1,1,1),omega.init=diag(c(100, 1, 1, 1, 1)), covariance.model = diag(c(1,0,0,0,1)))

ord.fit<-saemix(saemix.model,saemix.data,saemix.options)

# Case bootstrap
case.ordinal <- try(saemix.bootstrap(ord.fit, method="case", nboot=nboot))
if(is(case.ordinal,"data.frame"))
  write.table(case.ordinal, file.path(workDir, "results", "knee_caseBootstrap.res"), quote=F, row.names=FALSE)

# Conditional non-parametric bootstrap
cond.ordinal <- try(saemix.bootstrap(ord.fit, method="conditional", nboot=nboot))
if(is(cond.ordinal,"data.frame"))
  write.table(cond.ordinal, file.path(workDir, "results", "knee_condBootstrap.res"), quote=F, row.names=FALSE)
}

    covariate.model <- matrix(data=0, nrow=2, ncol=5)
    covariate.model[1,2]<-covariate.model[1,5]<-covariate.model[2,1]<-1
    ordmodel.cov<-saemixModel(model=ordinal.model,description="Ordinal categorical model",
        modeltype="likelihood",simulate.function=simulateOrdinal, 
        psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5, byrow=TRUE, dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), transform.par=c(0,1,1,1,1),
        omega.init=diag(c(100, 1, 1, 1, 1)), covariate.model=covariate.model, 
        covariance.model = diag(c(1,1,1,1,0)), verbose=FALSE)
    ord.fit.cov<-saemix(ordmodel.cov,saemix.data,saemix.options)

# Case bootstrap
case.ordinal <- try(saemix.bootstrap(ord.fit.cov, method="case", nboot=nboot))
if(is(case.ordinal,"data.frame"))
  write.table(case.ordinal, file.path(workDir, "results", "knee_caseBootstrapCov.res"), quote=F, row.names=FALSE)

# Conditional non-parametric bootstrap
cond.ordinal <- try(saemix.bootstrap(ord.fit.cov, method="conditional", nboot=nboot))
if(is(cond.ordinal,"data.frame"))
  write.table(cond.ordinal, file.path(workDir, "results", "knee_condBootstrapCov.res"), quote=F, row.names=FALSE)
