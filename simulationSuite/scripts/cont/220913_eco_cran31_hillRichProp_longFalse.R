# Folders
saemixDir <- "/home/eco/work/saemix/versions/saemix3.1"
saemixVersion <- "cran31"
simulDir <-"/home/eco/work/saemix/saemixextension"

# libraries
library(devtools)

# Details
today <- '220913'
who <- "eco"
settings <- "longFalse"
computeIndPar <- FALSE # can be set to TRUE to compute individual parameter estimates

# Scenario
scenario <- "hillRichProp"
datDir <- file.path(simulDir,"simulationSuite","cont", scenario,"data")
resDir <- file.path(simulDir,"simulationSuite","cont", scenario,"results")
namRes <- file.path(resDir, paste0(today,"_",who,"_",saemixVersion,"_",scenario,"_",settings,".res"))
nsim<-100
namsimdat<-"pdhillhigh"
parpop<-c(5,30,500,2)
nampar<-c("E0","Emax","ED50","gamma")
omega<-diag(c(0.09,0.49,0.49))
omega[3,2]<-omega[2,3]<-0.245
respar<-c(0.1)

pvrai <- c(parpop, diag(omega)[1:2], omega[2,3],omega[3,3], respar)
#pvrai<-c(5,30,500,0.09,0.490,0.245,0.490,sqrt(0.01))
pfaux<-c(10,60,1000,0.1,0.1,0.01,0.1,sqrt(0.0625))

if(settings=="defaultTrue") {
  psi0 <- pvrai
  omega0 <- diag(c(1,1,1,0.2))
  omega0[1:3,1:3]<-omega
  res0 <- respar
}
if(settings %in% c("defaultFalse", "longFalse")) {
  psi0 <- pfaux
  omega0 <- diag(c(1,1,1,1))
  res0 <- respar*2
}

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

set.seed(350922143)
zeseed <- trunc(runif(100, 0, 1)*10^7)

dev_mode() # development mode
install.packages(pkgs=file.path(saemixDir,"saemix_3.1.tar.gz"),repos=NULL)
library(saemix)

saemix.model<-saemixModel(model=modelhill, description="PD Hill model",
                          psi0=matrix(psi0, ncol=4, byrow = TRUE, dimnames=list(NULL, nampar)),transform.par=c(1,1,1,1),
                          covariance.model=matrix(c(1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0),ncol=4,byrow=TRUE),
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
hoh