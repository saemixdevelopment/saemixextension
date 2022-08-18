workDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(workDir,"R")

# Libraries
library(ggplot2)
library(MASS)

# Sourcing saemix functions
{
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
}

# Bootstrap functions
source(file.path(workDir,"bootstrap","saemix_bootstrap.R"))

datDir<-file.path(workDir,"bootstrap","data")
isim<-87
namscen<-"pdhillhigh.rich"
iscenar<-"Hill3"

mod.par<-c(E0=5,Emax=30,ED50=500,gamma=3)
mod.iiv<-c(E0=0.09,Emax=0.49,ED50=0.49,gamma=0)
omega<-diag(mod.iiv)
omega[3,2]<-omega[2,3]<-0.245
sigm<-0.1
modelhill<-function(psi,id,xidep) {
  # input:
  #   psi : matrix of parameters (4 columns, E0, Emax, E50, gamma)
  #   id : vector of indices 
  #   xidep : dependent variables (same nb of rows as length of id)
  # returns:
  #   a vector of predictions of length equal to length of id
  dose<-xidep[,1]
  e0<-psi[id,1]
  emax<-psi[id,2]
  e50<-psi[id,3]
  gamma<-psi[id,4]
  f<-e0+emax*dose**gamma/(e50**gamma+dose**gamma)
  return(f)
}

namfile<-file.path(datDir,paste("data_",namscen,"_sim",isim,".tab",sep=""))
saemix.data<-saemixData(name.data=namfile, header=T,
                        name.group=c("id"),name.predictors=c("time","amt"),
                        name.response=c("conc"),name.covariates=NULL,
                        units=list(x="hr",y="mg/L",covariates=c()))
saemix.model<-saemixModel(model=modelhill,description="Hill model", 
                          psi0=matrix(mod.par,ncol=4, byrow=TRUE,dimnames=list(NULL, c("E0","Emax","ED50","Gamma"))),transform.par=c(1,1,1,1),
                          covariance.model=matrix(c(1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0),ncol=4,byrow=TRUE),omega.init = omega,error.model="proportional", error.init = c(0,sigm))

saemix.options<-list(fix.seed=F,directory="current",displayProgress=FALSE, save.graphs=FALSE,print=FALSE)
saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
print(saemix.fit@results)

# om.estim<-c(diag(saemix.fit@results@omega)[1:3], saemix.fit@results@omega[2,3], saemix.fit@results@respar[2])
# par.estim<-c(saemix.fit@results@fixed.effects,om.estim)

nboot<-100
saemix.bootOpt<-list(fix.seed=F,directory="current",displayProgress=F, save.graphs=F, map=F, ll.is=F, print=FALSE)

start_time <- Sys.time()
boot.case<-saemix.bootstrap(saemix.fit, nboot=nboot, method="case")
case_time <- Sys.time()

print(case_time-start_time)

boot.cNP<-saemix.bootstrap(saemix.fit, nboot=nboot, method="conditional")
cNP_time <- Sys.time()
print(cNP_time-case_time)

boot.NP<-saemix.bootstrap(saemix.fit, nboot=nboot, method="residual")
NP_time <- Sys.time()
print(NP_time-cNP_time)

boot.Par<-saemix.bootstrap(saemix.fit, nboot=nboot, method="parametric")
Par_time <- Sys.time()
print(Par_time-NP_time)

end_time <- Sys.time()
