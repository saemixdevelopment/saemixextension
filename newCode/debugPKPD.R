# Folders
saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
progDirExt<-file.path(saemixDir,"Rext")
datDir<-file.path(saemixDir,"data")
datDir4<-file.path(saemixDir,"data40")
ecoDir<-file.path(saemixDir,"testeco")
belDir<-file.path(saemixDir,"testbelhal")

# Libraries
library(ggplot2)
library(MASS)

# Unchanged from version 3.0
source(file.path(progDir,"aaa_generics.R"))

############################## Data
# Extension
source(file.path(progDirExt,"SaemixData.R"))

# PK/PD
## data
x<-try(pkpd.saemix<-saemixData(name.data=file.path(datDir4,"warfarinPKPD.tab"), header=T,na=".", name.group=c("id"), name.predictors=c("time","amt"), name.response=c("dv"), name.ytype = "dvid", name.covariates=c("sex", "wt", "age"), units=list(x="hr",y="mg/L"), verbose=TRUE, outcome=c(conc="continuous", PCA="continuous")))

############################## Model
source(file.path(progDirExt,"SaemixModel.R"))
source(file.path(progDirExt,"SaemixOutcome.R"))

## model
model1cptdirect<-function(psi,id,xidep) { 
  tim<-xidep[,1]
  dose<-xidep[,2]
  ytype<-xidep$ytype
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  ic50<-psi[id, 4]
  k<-CL/V
  ypk<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  ypd<-100*(1-ypk/(ypk+ic50))
  ypk[ytype==2]<-ypd[ytype==2]
  return(ypk)
}

## Test the function
xidep<-pkpd.saemix@data[,c(pkpd.saemix@name.predictors, "ytype")]
id1<-pkpd.saemix@data$index
psi0<-c(ka=1, vd=5, cl=0.1, ic50=5)
psi1<-do.call(rbind, rep(list(psi0), pkpd.saemix@N))
ypred<-model1cptdirect(psi1, id1, xidep)

## saemixModel

############################## Algorithm

# Unchanged from version 3.0
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixRes.R"))
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

# Finding initial guesses for parameters

