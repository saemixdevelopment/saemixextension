########################################################### Setup
saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
progDirExt<-file.path(saemixDir,"Rext")
source(file.path(progDir,"aaa_generics.R"))

library(testthat)

# List of test models
source(file.path(saemixDir,"testecoExt","modelsTested.R"))

########################################################### Generate dataset with covariates and 3 occasions
source(file.path(saemixDir,"newCode","generatePKPDdataCovariateOcc.R"))

# Dataset pkpd.dat
head(pkpd.dat)
# parameters over each variability level as table
head(phimat)
# matching phimat and pkpd.dat
head(var.match)
# Outcome
xout

pkpd.dat$cens<-0
pkpd.dat$mdv<-0
x<-new(Class="SaemixData", name.data="pkpd.dat")
x@outcome<-xout
x@data<-pkpd.dat
x@name.group<-"id"
x@name.predictors<-c("tim","dose")
x@name.X<-"tim"
x@name.response<-c("y")
x@name.ytype<-"ytype"
x@name.covariates<-colnames(pkpd.dat)[7:12]
x@name.cens<-"cens"
x@name.mdv<-"mdv"
x@N<-length(pkpd.dat$id)
nind.obs<-tapply(pkpd.dat$id, pkpd.dat$id, length)
x@nind.obs<-c(nind.obs[match(unique(pkpd.dat$id),names(nind.obs))])
x@ntot.obs<-sum(x@nind.obs)
x@units<-c("hr","mg")
x@var.match <- var.match

saemix.data<-x

########################################################### Generating a model
# Model 2
imodel <- 2
saemix.model <- chooseDebugModel(model=imodel) # PK model 2 (without covariates, no IIV on ka, ka estimated)
psi0<-saemix.model@mu.start

# Model 11
imodel <- 11
saemix.model <- chooseDebugModel(model=imodel) # PK/PD model 11 (3 occasions, no period effect, covariates on both levels )
psi0<-saemix.model@mu.start

# Model 8 - with j0.covariate and i0.omega2 non null
imodel <- 8
saemix.model <- chooseDebugModel(model=imodel) # PK model 8 (covariates, no IIV on ka, ka estimated)
psi0<-saemix.model@mu.start

########################################################### saemix options
# Options and auxiliary files
source(file.path(saemixDir,"Rext","func_aux.R"))
source(file.path(saemixDir,"Rext","saemixControl.R"))

saemix.options<-saemixControl(nb.chains=2)
saemix.options$nbiter.tot<-sum(saemix.options$nbiter.saemix)

########################################################### Initialisation - Variability variables
# match variability levels
namvarlev<-c(saemix.model@var.model[[1]]@variable)
xlabel<-data.frame(level1=saemix.data@data[,saemix.data@name.group]) # labels
xindex<-data.frame(level1=saemix.data@data[,"index"])
if(nvarlevel==1) {
  saemix.model@nphirep<-list(rep(1,saemix.data@N))
  } else {
    saemix.model@nphirep<-vector(mode="list",length=saemix.model@nvarlevel)
    for(ivarlev in 2:saemix.model@nvarlevel) {
      namvarlev<-c(namvarlev, saemix.model@var.model[[ivarlev]]@variable)
      xlabel<- cbind(xlabel,
                     paste(xlabel[,(ivarlev-1)],saemix.data@data[,saemix.model@var.model[[ivarlev]]@variable],sep="#")
    }
    colnames(xlabel)<-paste0("level",1:nvarlevel)
    for(ivarlev in 2:saemix.model@nvarlevel)
      xindex[,ivarlev]<-match(xlabel[,ivarlev], unique(xlabel[,ivarlev]))
    colnames(xindex)<-paste0("level",1:nvarlevel)
    xlabel<-xlabel[!duplicated(xlabel[,nvarlevel]),]
    for(ivarlev in 1:saemix.model@nvarlevel)
          saemix.model@nphirep[[ivarlev]]<-c(tapply(xlabel[,ivarlev], xlabel[,ivarlev],length))
}

# need to add xindex to the model object, can be used to obtain model predictions for different levels of variability included

########################################################### Initialisation - WIP
source(file.path(saemixDir,"newCode","extractDesignMatrix.R"))

# Model structure
for(ivarlev in 1:length(saemix.model@var.model)) {
  newlist <- extractModelStructure(saemix.data, saemix.model, level=ivarlev, verbose=saemix.options$warnings)
  saemix.model@ind.model <- list(newlist$indiv.model)
}
Uargs2 <- newlist$Uargs
indmodel <- newlist$indiv.model

# Initialise model parameters and update omega matrices
newlist <- initialisePhi(saemix.data, saemix.model, nb.chains=saemix.options$nb.chains, verbose=saemix.options$warnings)
saemix.model <- newlist$saemix.model
phiM <- newlist$phiM
Dargs <- newlist$Dargs
DYF <- newlist$DYF

var.eta<-mydiag(indiv.model@omega)
theta0<-c(fixedpsi.ini,var.eta[i1.omega2],pres[indx.res])
l1<-betas.ini
l1[indx.betaI]<-transphi(matrix(l1[indx.betaI],nrow=1),saemix.model["transform.par"])
allpar0<-c(l1,var.eta[i1.omega2],pres[indx.res])

########################################################### 
