#####################################################################################
# Setting up the library

saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")

# Libraries
library(ggplot2)
library(MASS)

source(file.path(progDir,"aaa_generics.R"))
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

# Folder for results
namsimul<-"countRich50risk"
saveDir<-file.path(saemixDir,"ecoSimul",namsimul)

#####################################################################################
# Pb with eg isim=22
pbsim<-c(22,41,46,60,66,112,114,170,179)
isim<-pbsim[3]

nsim<-200
set.seed(810212100) # Seed for the random number generator
dataseed<-trunc(runif(nsim, 1, 10^8)) # seeds used to generate data
simseed<-trunc(runif(nsim, 1, 10^8)) # seeds used to run saemix (same seed for both sets of CI)

# Settings  
param <- c(39.1, 0.0388, 0.1 )
omega<-c(0.5, 0.5) # SD=50%
paramSimul<-c(param, omega)
parnam<-c("alpha","beta","risk","omega.alpha","omega.beta")

nsuj<-40
xtim<-c(0.0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)

# Model
countData.model<-function(psi,id,xidep) {
  tim <- xidep[,1]
  y <- xidep[,2] 
  alpha <- psi[id,1]
  beta <- psi[id,2]
  lambda <- alpha*exp(-beta*tim)
  
  logpdf <- rep(0,length(tim))
  logpdf <- -lambda + y*( (log(alpha) - beta*tim )) - log(factorial(y))
  return(logpdf) 
}

# Model
#                          omega.init = diag(c(paramSimul[4]*3, paramSimul[5]**3)),
#                          transform.par=c(0,0),covariance.model=matrix(c(1,0,0,1),ncol=2))

saemix.model.true<-saemixModel(model=countData.model,description="Count data model", modeltype="likelihood",
                               psi0=matrix(c(param[1:2],0,param[3]),ncol=2,byrow=TRUE,dimnames=list(NULL,parnam[1:2])),
                               covariate.model=matrix(c(0,1),ncol=2), omega.init = diag(c(0.5,0.5)),
                               transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2))

set.seed(dataseed[isim])
partab<-as.data.frame(matrix(data=0,nrow=nsuj,ncol=2,dimnames=list(NULL,parnam[1:2])))
for(i in 1:2) partab[,i]<-rnorm(nsuj,mean=log(param[i]),sd=omega[i])
partab[(1+nsuj/2):nsuj,2]<-partab[(1+nsuj/2):nsuj,2]+param[3]
for(i in 1:2) partab[,i]<-exp(partab[,i])

psim<-data.frame()
for(itim in xtim) {
  lambda<-partab[,1]*exp(-partab[,2]*itim)
  psim<-rbind(psim,lambda)
}
datsim<-data.frame(id=rep(1:nsuj,each=length(xtim)),time=rep(xtim,nsuj),lambda=unlist(psim))
rownames(datsim)<-NULL
ysim<-rpois(dim(datsim)[1], lambda=datsim$lambda)
#  summary(datsim)
datsim$y<-ysim
datsim$risk<-ifelse(datsim$id>(nsuj/2),1,0)

# Running saemix
saemix.data<-saemixData(name.data=datsim,
                        name.group=c("id"),name.predictors=c("time","y"), name.covariates=c("risk"),name.X=c("time"))

saemix.options<-list(seed=simseed[isim],save=FALSE,save.graphs=FALSE, fim=FALSE, displayProgress=FALSE)
count.fit<-try(saemix(saemix.model.true,saemix.data,saemix.options))

ypred<-predict.saemixmodel(saemix.model.true, predictors=saemix.data@data[,c("time","y")], psi=NA, id=saemix.data@data[,c("id")])
summary(ypred$predictions)
idx<-unique(ypred$predictions$id[is.infinite(ypred$predictions$pred)])

# 
# object<-saemix.model.true
# predictors<-saemix.data@data[,c("time","y")]
# psi<-NA
# id<-saemix.data@data[,c("id")]

saemix.data<-saemixData(name.data=datsim[!(datsim$id %in% idx),],
                        name.group=c("id"),name.predictors=c("time","y"), name.covariates=c("risk"),name.X=c("time"))
ypred<-predict.saemixmodel(saemix.model.true, predictors=saemix.data@data[,c("time","y")], psi=NA, id=saemix.data@data[,c("id")])
summary(ypred$predictions)
count.fit<-try(saemix(saemix.model.true,saemix.data,saemix.options))

