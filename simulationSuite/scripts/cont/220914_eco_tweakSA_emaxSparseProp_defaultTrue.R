# Folders
saemixDir <- "/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
saemixVersion <- "tweakSA"
simulDir <-"/home/eco/work/saemix/saemixextension"

# Libraries
library(ggplot2)
library(MASS)
library(rlang)

# Source functions
source(file.path(progDir,"aaa_generics.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"func_aux.R"))
source(file.path(progDir,"main_initialiseMainAlgo.R"))

# Main files from package
#source(file.path(progDir,"main.R"))
source(file.path(progDir,"main_estep.R"))
#source(file.path(progDir,"main_mstep.R"))
source(file.path(progDir,"func_FIM.R"))
source(file.path(progDir,"func_plots.R"))
source(file.path(progDir,"func_distcond.R"))
source(file.path(progDir,"func_simulations.R"))
source(file.path(progDir,"compute_LL.R"))
source(file.path(progDir,"func_estimParam.R"))

# Testing alternate versions
source(file.path(saemixDir,"newCode","alternateMain.R"))
source(file.path(saemixDir,"newCode","alternate_mstep.R"))

# Details
today <- '220914'
who <- "eco"
settings <- "defaultTrue"
computeIndPar <- FALSE # can be set to TRUE to compute individual parameter estimates
set.seed(1409221901)
zeseed <- trunc(runif(100, 0, 1)*10^7)

# Scenario
scenario <- "emaxSparseProp"
datDir <- file.path(simulDir,"simulationSuite","cont", scenario,"data")
resDir <- file.path(simulDir,"simulationSuite","cont", scenario,"results")
namRes <- file.path(resDir, paste0(today,"_",who,"_",saemixVersion,"_",scenario,"_",settings,".res"))
nsim<-100
namsimdat<-"pdemax"
parpop<-c(5,30,500)
nampar<-c("E0","Emax","ED50")
omega<-diag(c(0.09,0.49,0.49))
omega[3,2]<-omega[2,3]<-0.245
respar<-c(0.1)
pvrai <- c(parpop, diag(omega)[1:2], omega[2,3],omega[3,3], respar)
pfaux<-c(10,60,1000,0.1,0.1,0.01,0.1,sqrt(0.0625))

if(settings=="defaultTrue") {
  psi0 <- pvrai[1:3]
  omega0 <- omega
  res0 <- respar
}
if(settings %in% c("defaultFalse", "longFalse")) {
  psi0 <- pfaux[1:3]
  omega0 <- diag(c(1,1,1))
  res0 <- respar*2
}

modelemax<-function(psi,id,xidep) {
  # input:
  #   psi : matrix of parameters (3 columns, E0, Emax, E50)
  #   id : vector of indices 
  #   xidep : dependent variables (same nb of rows as length of id)
  # returns:
  #   a vector of predictions of length equal to length of id
  dose<-xidep[,1]
  e0<-psi[id,1]
  emax<-psi[id,2]
  e50<-psi[id,3]
  f<-e0+emax*dose/(e50+dose)
  return(f)
}

saemix.model<-saemixModel(model=modelemax, description="PD Emax model",
                          psi0=matrix(psi0, ncol=3, byrow = TRUE, dimnames=list(NULL, nampar)),transform.par=c(1,1,1),
                          covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),
                          omega.init = omega0, error.init=c(0,res0), error.model="proportional")
saemix.options<-list(seed=zeseed[1], save=FALSE, save.graphss=FALSE, displayProgress=FALSE)
if(settings=="longFalse") {
  saemix.options$nchains <- 5
  saemix.options$nbiter.saemix <- c(800, 300)
}

isim <- 1
for(isim in 1:nsim) {
  # data
  namfich<-paste('data_',namsimdat,isim,".tab",sep="")
  saemix.data<-saemixData(file.path(datDir,namfich),header=T,name.group="id",name.predictors="dose",name.response="y",units=list(x="-",y="-"))
  
  yfit <- saemixAlternate(saemix.model, saemix.data, saemix.options)
  if(computeIndPar)
    yfit <- conddist.saemix(yfit)
  # Saving population parameters
  idx.eps<-yfit@model@indx.res
  idx.iiv<-yfit@model@indx.omega
  vec<-summary(yfit@results)
  l1<-unlist(c(isim,vec$fixed.effects[,2],vec$random.effects[,2],vec$logLik[1:2,2]))
  l2<-unlist(c(vec$fixed.effects[,3],vec$random.effects[,3]))
  if(isim==1) {
    headersOrig<-c("Simulation",vec$fixed.effects[,1],as.character(vec$random.effects[,1]),"LL.lin","LL.IS")
    headersSE<-headersOrig[1:(length(headersOrig)-2)]
    headersSE[-c(1)]<-paste("SE.", headersSE[-c(1)],sep="")
    headersOrig <- c(headersOrig, headersSE[-c(1)])
    write(headersOrig,namRes,ncol=length(headersOrig))
  }
  write(c(l1,l2),namRes,ncol=length((headersOrig)),append=T)
}

