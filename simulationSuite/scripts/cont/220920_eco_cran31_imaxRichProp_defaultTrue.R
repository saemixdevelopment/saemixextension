# Folders
saemixDir <- "/home/eco/work/saemix/versions/saemix3.1"
saemixVersion <- "cran31"
simulDir <-"/home/eco/work/saemix/saemixextension"

# libraries
library(devtools)

# Details
today <- '220920'
who <- "eco"
settings <- "defaultTrue"
computeIndPar <- FALSE # can be set to TRUE to compute individual parameter estimates

# Scenario
scenario <- "imaxRichProp"
datDir <- file.path(simulDir,"simulationSuite","cont", scenario,"data")
resDir <- file.path(simulDir,"simulationSuite","cont", scenario,"results")
namRes <- file.path(resDir, paste0(today,"_",who,"_",saemixVersion,"_",scenario,"_",settings,".res"))

nsim<-100
namsimdat<-"pdimax"
parpop<-c(100,0.7,500,2)
parcov<-(-0.5) # ED50=500 in group 0, 300 in group 1
nampar<-c("E0","Imax","ED50","gamma")
omega<-diag(c(0.01,0,0.09,0)) # 10% IIV on E0 and 30% IIV on ID50, no correlation
respar<-c(5,0.2) # additive and proportional error terms

sigma<-respar
if(length(grep("Add", scenario))>0) sigma<-respar[1]
if(length(grep("Prop", scenario))>0) sigma<-respar[2]
pvrai <- c(parpop[1:3], parcov, parpop[4], omega[1,1],omega[3,3], respar)
pfaux <-c(50,0.5,100,0,1,1,1,1,0.5)

if(settings=="defaultTrue") {
  psi0 <- pvrai[c(1:3,5)]
  omega0 <- omega
  res0 <- sigma
}
if(settings %in% c("defaultFalse","longFalse")) {
  psi0 <- pfaux[c(1:3,5)]
  omega0 <- diag(c(1,1,1))
  res0 <- sigma*2
}

sigmoidImax<-function(psi,id,xidep) {
  # input:
  #   psi : matrix of parameters (3 columns, E0, Emax, E50)
  #   id : vector of indices 
  #   xidep : dependent variables (same nb of rows as length of id)
  # returns:
  #   a vector of predictions of length equal to length of id
  dose<-xidep[,1]
  e0<-psi[id,1]
  imax<-psi[id,2]
  ed50<-psi[id,3]
  gamma<-psi[id,4]
  f<-e0*(1-imax*(dose**gamma)/((ed50**gamma)+(dose**gamma)))
  return(f)
}

# Random generator seed
set.seed(860922741)
zeseed <- trunc(runif(nsim, 0, 1)*10^7)

# Estimation
dev_mode() # development mode
install.packages(pkgs=file.path(saemixDir,"saemix_3.1.tar.gz"),repos=NULL)
library(saemix)
library(saemix)

saemix.model.add<-saemixModel(model=sigmoidImax, description="PD Imax model",
                              psi0=matrix(psi0, ncol=4, byrow = TRUE, dimnames=list(NULL, nampar)),transform.par=c(1,1,1,1),
                              covariance.model=diag(c(1,0,1,0)), covariate.model=c(0,0,1,0),
                              omega.init = omega0, error.init=c(sigma,0), error.model="additive")
saemix.model.prop<-saemixModel(model=sigmoidImax, description="PD Imax model",
                               psi0=matrix(psi0, ncol=4, byrow = TRUE, dimnames=list(NULL, nampar)),transform.par=c(1,1,1,1),
                               covariance.model=diag(c(1,0,1,0)),covariate.model=c(0,0,1,0),
                               omega.init = omega0, error.init=c(0,sigma), error.model="proportional")
if(length(grep("Add", scenario))>0) saemix.model<-saemix.model.add
if(length(grep("Prop", scenario))>0) saemix.model<-saemix.model.prop

saemix.options<-list(seed=zeseed[1], save=FALSE, save.graphs=FALSE, displayProgress=FALSE)
if(settings=="longFalse") {
  saemix.options$nchains <- 5
  saemix.options$nbiter.saemix <- c(800, 300)
}

# testing
# saemix.data<-saemixData(xdat,header=T,name.group="id",name.predictors="dose",name.response="eff", name.covariate="group",units=list(x="-",y="-"))
# yfit <- saemix(saemix.model, saemix.data, saemix.options)

isim <- 1
for(isim in 1:nsim) {
  # data
  namfich<-paste('data_',namsimdat,isim,".tab",sep="")
  saemix.data<-saemixData(file.path(datDir,namfich),header=T,name.group="id",name.predictors="dose",name.response="eff", name.covariate='group',units=list(x="-",y="-"))
  
  yfit <- saemix(saemix.model, saemix.data, saemix.options)
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

dev_mode() # development mode
