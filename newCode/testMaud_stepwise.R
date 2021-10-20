saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
datDir<-file.path(saemixDir,"data")
ecoDir<-file.path(saemixDir,"testeco")
belDir<-file.path(saemixDir,"testbelhal")

# Libraries
library(ggplot2)
library(MASS)

# Loading package functions
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

# Stepwise algorithm
source(file.path(progDir,"backward.R"))
source(file.path(progDir,"forward.R"))
source(file.path(progDir,"stepwise.R"))
source(file.path(progDir,"func_compare.R"))
source(file.path(progDir,"func_stepwise.R"))


# Running stepwise algorithm for covariate selection on theophylline data
theo.saemix<-read.table(file.path(datDir, "theo.saemix.tab"), header=T)

saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
                        name.group=c("Id"),name.predictors=c("Dose","Time"),
                        name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
                        units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

# Definition of models to be compared
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

saemix.model1<-saemixModel(model=model1cpt,modeltype="structural", 
                           description="One-compartment model with first-order absorption",
                           psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))), 
                           transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, warnings=FALSE)
saemix.fit1<-saemix(saemix.model1,saemix.data,saemix.options)

covariate.init <- matrix(c(1,0,0,0,1,0),ncol=3,nrow=2)
res.forward <- step.saemix(saemix.fit1, direction = "forward")
res.backward <- step.saemix(saemix.fit1, direction = "backward", covariate.init=covariate.init)
res.stepwise <- step.saemix(saemix.fit1, direction="both", covariate.init=covariate.init)
