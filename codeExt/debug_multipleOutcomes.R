# Folders
saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
progDirExt<-file.path(saemixDir,"Rext")
datDir<-file.path(saemixDir,"data")

# Files
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDirExt,"SaemixData.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDirExt,"SaemixModel.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"func_aux.R"))
source(file.path(progDirExt,"SaemixOutcome.R"))

# Other files
if(FALSE) {
  source(file.path(progDir,"main.R"))
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
}

# Setup data and model

## 2 continuous responses
pkpd.saemix<-read.table(file.path(saemixDir,"data40","warfarinPKPD.tab"),header=T,na=".")
pkpd.dat<-saemixData(name.data=pkpd.saemix, header=T,na=".", name.group=c("id"), name.predictors=c("time","amt"), name.response=c("dv"), name.ytype = "dvid", name.covariates=c("sex", "wt", "age"), units=list(x="hr",y="mg/L"), verbose=TRUE, outcome=c(conc="continuous", PCA="continuous"))
summary(pkpd.dat@data)

## Model function
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

## Model predictions for pkpd 
xidep<-pkpd.dat@data[,c(pkpd.dat@name.predictors, "ytype")]
id1<-pkpd.dat@data$index
psi0<-c(ka=1, vd=5, cl=0.1, ic50=5)
psi1<-do.call(rbind, rep(list(psi0), pkpd.dat@N))
ypred<-model1cptdirect(psi1, id1, xidep)

# Error model

## Outcomes
error.model<-c('proportional','power')
out1<-new(Class="SaemixContinuousOutcome", name.outcome=names(pkpd.dat@outcome)[1], error.model="proportional", error.parameters=c(0.2))
out2<-new(Class="SaemixContinuousOutcome", name.outcome=names(pkpd.dat@outcome)[2], error.model="power", error.parameters=c(5,0.3,2))
showall(out2)
nr<-2
outcomes<-c(out1, out2)

## Characteristics
name.error.parameters<-c(out1@error.parameters,out2@error.parameters)
error.npar<-c(out1@error.npar, out2@error.npar)
error.parameters<-c(0.2, 5, 0.3, 2)
ind.res<-1:4
ind.res.fix<-4
index.respar<-list(c(1),c(2:4))

## gpred
gpred1<-gpred2<-ypred
for(i in 1:nr) {
  fpred<-ypred[pkpd.dat@data$ytype==i]
  gpred1[pkpd.dat@data$ytype==i]<-outcomes[[i]]@error.function(fpred, error.parameters[index.respar[[i]]])
  gpred2[pkpd.dat@data$ytype==i]<-outcomes[[i]]@error.function(fpred, outcomes[[i]]@error.parameters)
}
summary(gpred1-gpred2)

# saemixModel with outcomes

# Compute 


# Initialise

