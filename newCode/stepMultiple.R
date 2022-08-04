########################################################### Folder and generics
saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
progDirExt<-file.path(saemixDir,"Rext")
source(file.path(progDir,"aaa_generics.R"))

library(testthat)

########################################################### Data
# Creating an SaemixData object with PK/PD data by hand - Responses
source(file.path(progDirExt,"SaemixOutcome.R"))
source(file.path(progDirExt,"SaemixOutcome-methods.R"))
source(file.path(progDirExt,"SaemixData.R"))
source(file.path(progDirExt,"SaemixCovariateModel.R"))
source(file.path(progDirExt,"SaemixCovariate.R"))
pkpd<-read.table(file.path(saemixDir,"data40","warfarinPKPD.tab"), header=T)

x<-new(Class="SaemixData", name.data="pkpd")
x1<-continuousOutcome(model="combined2", start=c(1, 0.2))
x2<-continuousOutcome(start=c(2))
#new(Class="SaemixOutcome",type="continuous")
x@outcome<-list(conc=new(Class="SaemixOutcome",type="continuous"), effect=new(Class="SaemixOutcome",type="continuous"))

zesuj<-unique(pkpd$id)
ind1<-match(pkpd$id, zesuj)
x@data<-data.frame(index=ind1, pkpd[,c(1:4)], ytype=pkpd$dvid, cens=0, mdv=0, pkpd[,c(6:8)])
x@name.group<-"id"
x@name.predictors<-c("time","amt")
x@name.X<-"time"
x@name.response<-c("dv")
x@name.ytype<-"ytype"
x@name.covariates<-colnames(pkpd)[6:8]
x@name.cens<-"cens"
x@name.mdv<-"mdv"
x@N<-length(zesuj)
nind.obs<-tapply(pkpd$id, pkpd$id, length)
x@nind.obs<-c(nind.obs[match(zesuj,names(nind.obs))])
x@ntot.obs<-sum(x@nind.obs)
x@units<-c("hr","mg")

saemix.data<-pkpd.saemix<-x

########################################################### Model
# Model function
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

# Mean value for population parameters
pkpd.psi0<-c(ka=1, vd=5, cl=0.1, ic50=5)
# Table with population parameters for every subject in the dataset)
psi1<-do.call(rbind, rep(list(pkpd.psi0), pkpd.saemix@N))

## Test the function
xidep<-pkpd.saemix@data[,c(pkpd.saemix@name.predictors, "ytype")]
id1<-pkpd.saemix@data$index
ypred<-model1cptdirect(psi1, id1, xidep)

# Model outcomes
pkpd.outcome<-list(conc=saemixOutcome(unit="mg/L", model="combined2", start=c(1, 0.2)), 
                   effect=saemixOutcome(unit="%", start=1))

# Test one error model on all the predictions
f1<-runif(10, 0, 10)
g1<-pkpd.outcome[[1]]@error.function(f1, pkpd.outcome[[1]]@error.parameters)
if(max(abs(g1-sqrt((f1*x1$start[2])**2+x1$start[1]**2)))>10^(-10)) cat("Problem in computing g\n")

########################################################### SaemixModel
# Model class
source(file.path(progDirExt,"SaemixVarLevel.R"))
source(file.path(progDirExt,"SaemixParameter.R"))
source(file.path(progDirExt,"SaemixParameter-methods.R"))
source(file.path(progDirExt,"SaemixModel.R"))
source(file.path(progDirExt,"SaemixModel-methods.R"))

# List of test models
if(FALSE) {
  source(file.path(saemixDir,"testecoExt","modelsTested.R"))
  pkpd.model <- chooseDebugModel(model=7)
}

# Model parameters
lpar <- list(ka=saemixPar(mu.start=pkpd.psi0[1], omega.start=0.3),
             cl=saemixPar(mu.start=pkpd.psi0[2], omega.start=0.5, covariate=c(sex=binCov(beta=0.5))),
             vd=saemixPar(mu.start=pkpd.psi0[3], omega.start=0.5, covariate=c(wt=contCov(beta=1, beta.fix=1)), rho.param=c("cl"), rho=0.5),
             ic50=saemixPar(mu.start=pkpd.psi0[4]))


pkpd.model<-saemixModel(model=model1cptdirect, outcome=pkpd.outcome, parameter=lpar)

saemix.model<-pkpd.model

########################################################### saemix options
# Options and auxiliary files
source(file.path(saemixDir,"Rext","func_aux.R"))
source(file.path(saemixDir,"Rext","saemixControl.R"))

saemix.options<-saemixControl(nb.chains=2)
saemix.options$nbiter.tot<-sum(saemix.options$nbiter.saemix)

########################################################### Old initialisation (to check)
source(file.path(progDirExt,"main_initialiseMainAlgo_old.R"))
mylist<-initialiseMainAlgo.old(saemix.data, saemix.model, saemix.options)
if(FALSE) {
  print(names(mylist))
  print(names(mylist$Uargs))
  print(names(mylist$Dargs))
  print(names(mylist$varList))
}

########################################################### Predictions and error
# ssq, compute.LLy

# Creating input for compute.LLy using several Markov chains
## code from initialise
nchains<-saemix.options$nb.chains
chdat<-new(Class="SaemixRepData",data=pkpd.saemix, nb.chains=nchains)
NM<-chdat["NM"]
IdM<-chdat["dataM"]$IdM
yM<-chdat["dataM"]$yM
XM<-chdat["dataM"][,c(pkpd.saemix["name.predictors"],pkpd.saemix["name.cens"],pkpd.saemix["name.mdv"],pkpd.saemix["name.ytype"]),drop=FALSE]
io<-matrix(data=0,nrow=pkpd.saemix["N"],ncol=max(pkpd.saemix["nind.obs"]))
for(i in 1:pkpd.saemix@N)
  io[i,1:pkpd.saemix["nind.obs"][i]]<-1
ioM<-do.call(rbind,rep(list(io),nchains))
ind.ioM <- which(t(ioM)!=0)
DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])

phiM<-psiM<-do.call(rbind, rep(list(pkpd.psi0), pkpd.saemix@N*nchains)) # h=id() so phiM=psiM

##########################
# Definition of Dargs
## elements extracted from the model and from the duplicated data through chains
## Note: error.parameters could be also in Dargs (model specific)
Dargs<-list(IdM=IdM, XM=XM, yM=yM, NM=NM, ind.ioM=ind.ioM, 
            model=saemix.model@model, transform.par=saemix.model@transform.par, outcome=saemix.model@outcome)

# do we need NM ? yes probably
# check if need nobs, yobs, is.lpdfM

# Checking new Dargs definition w/r to previous elements
test_that("Replicated data structure, ind.ioM moved to Dargs", {
  expect_identical(mylist$Dargs$IdM, Dargs$IdM)
  expect_identical(mylist$Dargs$yM, Dargs$yM)
  expect_identical(mylist$Dargs$XM, Dargs$XM)
  expect_equal(mylist$Dargs$NM, Dargs$NM)
  expect_null(mylist$Dargs$ind.ioM)
  expect_identical(mylist$Uargs$ind.ioM, Dargs$ind.ioM)
})
test_that("Model elements", {
  expect_identical(mylist$Dargs$transform.par, Dargs$transform.par)
  expect_identical(mylist$Dargs$structural.model, Dargs$model)
  for(i in 1:length(Dargs$outcome)) expect_identical(mylist$Dargs$error.model[i], Dargs$outcome[[i]]@error.model)
})

##########################
# debuging ssq, compute.LLy, optim
if(FALSE)
  source(file.path(saemixDir,"newCode","debug_ssqAcceptanceRatio.R"))
xvec2<-compute.LLy(phiM, Dargs, DYF)


##########################
#### Eco 12/04/2022 - stopped here
# Definition
Uargs<-list(nchains=nchains, i0.omega2=i0.omega2, j0.covariate=j0.covariate,
            MCOV0=MCOV0, COV0=COV0)

# Assume that the covariate model applies to all variability levels
## normally, we would determine which covariates apply to each covariate level
var.model<-saemix.model@var.model[[1]]
covariate.model<-saemix.model@covariate.model
level<-1

# index.XXX: index in a vector
# idxmat.XXX: index in a matrix

########################################################### Variability levels
########################## Two variability levels

omega1<-vec2mat(c(1,0.9,0,1,0,0.5))
omega2<-vec2mat(c(0.5,0.2,0,0.5,0,0))
var.list<-list(saemixVarModel(omega=omega1, variable="id"), saemixVarModel(omega=omega2, variable="occ", verbose=T))

repmatrix <- function(mat, nrep) {
  matrep <- NULL
  for(icol in 1:ncol(mat))
    matrep<-cbind(matrep, rep(mat[,icol], times=nrep))
  return(matrep)
#  return(as.data.frame(matrep))
}

par.id <- data.frame(ka=seq(1,2,length.out = 5), cl=seq(5,10, length.out=5), vd=seq(20,50, length.out=5),id=1:5)
nb.occ <- c(2,3,1,2,4)
par.occ<- repmatrix(par.id, nrep=nb.occ)
colnames(par.occ)<-colnames(par.id)
par.occ<-as.data.frame(par.occ)
par.occ <- cbind(par.occ, idocc=paste0(par.id[,4],"-",sequence(nb.occ)))

########################################################### Used in main_mstep.R to optimise beta0 (???)
# Fixed parameters ? (ie covariate effects and parameters without IIV ?)

# Input
## b0: parameters to optimise on in M-step (size args$j0.covariate)
## phiM: matrix of parameters
## args: list, here we use the elements 
### MCOV0
### COV0
### j0.covariate
### i0.omega2
### nchains
## Dargs: list, here we use the elements 
### transform.par, model (structural model), 
### IdM, XM, yM, ind.ioM (a vector of 0/1 indicating which elements of the square matrix DYF have data), 
### list of outcomes (type, error.model, error.function)
## error.parameters: list of the error parameters for the different outcomes in the model (NULL for non-continuous outcome)
## DYF: a matrix of dimension nmax=max(n_i) times (N*nb.chains); the column for subject i has the last nmax-n_i elements set to 0, so that colSums sums on the observations for that subject 
# Output 
## U: vector with either sum(-LL) (continuous models, minus the constant sum(log(1/sqrt(2*pi)))) or sum(-logpdf) (likelihood models)


compute.Uy<-
  #function(b0,phiM,pres,args,Dargs,DYF) { # old
  function(b0, phiM, Dargs, DYF) {
  # Attention, DYF variable locale non modifiee en dehors
  args$MCOV0[args$j0.covariate]<-b0
  phi0<-args$COV0 %*% args$MCOV0
  phiM[,args$i0.omega2]<-do.call(rbind,rep(list(phi0),args$nchains))
  
  psiM<-transphi(phiM,Dargs$transform.par)
  dyf.change<-DYF[args$ind.ioM]
  fpred<-Dargs$model(psiM,Dargs$IdM,Dargs$XM)
  lpred<-fpred
  for(iout in 1:length(Dargs$outcome)) {
    idx1<-which(Dargs$XM$ytype==iout)
    if(Dargs$outcome[[iout]]@type=="continuous") {
      if(Dargs$outcome[[iout]]@error.model=="exponential") fpred[idx1]<-log(cutoff(fpred[idx1]))
      gpred<-Dargs$outcome[[iout]]@error.function(fpred[idx1], Dargs$error.parameters[[iout]])
      dyf.change[idx]<-0.5*((Dargs$yM[idx1]-fpred[idx1])/gpred[idx1])**2+log(gpred[idx1])
    } else {
      dyf.change[idx]<-(-fpred[idx1])
    }
  }
  DYF[args$ind.ioM]<- dyf.change
  U<-sum(DYF)
  return(U)
  }





#   Uargs<-list(nchains=saemix.options$nb.chains,nb.parameters=nb.parameters, nb.betas=nb.betas, nb.etas=nb.etas, 
#               nb.parest=nb.parest,indx.betaC=indx.betaC, indx.betaI=indx.betaI, ind.res=ind.res,
#               indest.omega=indest.omega, i0.omega2=i0.omega2, i1.omega2=i1.omega2,	j.covariate=j.covariate, 
#               ind.fix10=ind.fix10, ind.fix11=ind.fix11, ind.fix1=ind.fix1, ind.fix0=ind.fix0,
#               MCOV0=MCOV0, COV=COV, COV0=COV0, COV1=COV1, LCOV=LCOV, COV2=COV2, dstatCOV=dstatCOV, 
#               Mcovariates=Mcovariates, ind.ioM=ind.ioM)


################################### Testing
# Creating input for compute.LLy.mat
## from initialise
# using several Markov chains
nchains<-2
chdat<-new(Class="SaemixRepData",data=pkpd.saemix, nb.chains=nchains)
NM<-chdat["NM"]
IdM<-chdat["dataM"]$IdM
yM<-chdat["dataM"]$yM
XM<-chdat["dataM"][,c(pkpd.saemix["name.predictors"],pkpd.saemix["name.cens"],pkpd.saemix["name.mdv"],pkpd.saemix["name.ytype"]),drop=FALSE]
io<-matrix(data=0,nrow=pkpd.saemix["N"],ncol=max(pkpd.saemix["nind.obs"]))
for(i in 1:pkpd.saemix@N)
  io[i,1:pkpd.saemix["nind.obs"][i]]<-1
ioM<-do.call(rbind,rep(list(io),nchains))
ind.ioM <- which(t(ioM)!=0)
DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])

phiM<-psiM<-do.call(rbind, rep(list(pkpd.psi0), pkpd.saemix@N*nchains)) # h=id() so phiM=psiM

# Dargs for initial version
Dargs1<-list(transform.par=c(0,0,0,0), IdM=IdM, XM=XM, yM=yM, modeltype="structural", etype.exp=integer(0), structural.model=pkpd.model@model) # normal distribution for parameters
pres1<-c(pkpd.model@outcome[[1]]@error.parameters,pkpd.model@outcome[[2]]@error.parameters,0) # old version: needs 2 parameters for each error model
xvec1<-compute.LLy.old(phiM,args=list(ind.ioM=ind.ioM),Dargs1,DYF,pres1)

# Dargs for new version
Dargs2<-list(transform.par=c(0,0,0,0), IdM=IdM, XM=XM, yM=yM, ind.ioM=ind.ioM, model=pkpd.model@model, outcome=pkpd.model@outcome)
pres2<-list(pkpd.model@outcome[[1]]@error.parameters, pkpd.model@outcome[[2]]@error.parameters)
xvec2<-compute.LLy(phiM, Dargs2, DYF, pres2)
if(max(abs(xvec1-xvec2))>0) cat("Mismatch between old and new versions\n")



########################################################### Used in main_estep to optimise phi1.opti

conditional.distribution_c<-function(phi1,phii,idi,xi,yi,mphi,idx,iomega,trpar,model,pres,err) {
  phii[idx]<-phi1
  psii<-transphi(matrix(phii,nrow=1),trpar)
  if(is.null(dim(psii))) psii<-matrix(psii,nrow=1)
  fi<-model(psii,idi,xi)
  #  if(err=="exponential") # Reverted this bit to the previous version to avoid a compiler error, not sure why it was changed...
  #    fi<-log(cutoff(fi))
  #  if (!(is.null(pres)) && pres[1] == pres) { # package compile throws an error when comparing a vector of length 2 (pres) to a vector of length 1
  # gi <- cutoff(pres[1])
  # } 
  # else{
  #   gi<-error(fi,pres) #    cutoff((pres[1]+pres[2]*abs(fi)))
  # }
  ind.exp<-which(err=="exponential")
  for(ityp in ind.exp) 
    fi[xi$ytype==ityp]<-log(cutoff(fi[xi$ytype==ityp]))
  gi<-error(fi,pres,xi$ytype)      #    cutoff((pres[1]+pres[2]*abs(fi)))
  Uy<-sum(0.5*((yi-fi)/gi)**2+log(gi))
  dphi<-phi1-mphi
  Uphi<-0.5*sum(dphi*(dphi%*%iomega))
  return(Uy+Uphi)
}

conditional.distribution_d<-function(phi1,phii,idi,xi,yi,mphi,idx,iomega,trpar,model) {
  phii[idx]<-phi1
  psii<-transphi(matrix(phii,nrow=1),trpar)
  if(is.null(dim(psii))) psii<-matrix(psii,nrow=1)
  fi<-model(psii,idi,xi)
  Uy <- sum(-fi)
  dphi<- phi1-mphi
  Uphi<- 0.5*sum(dphi*(dphi%*%iomega))
  return(Uy+Uphi)
}


########################################################### One longitudinal + one event

x1<-continuousOutcome(model="proportional", start=c(0.3))
x2<-discreteOutcome(type="event")
x3<-discreteOutcome()
xout<-new(Class="SaemixContinuousOutcome", error.model=x1$error.model, error.npar=x1$error.npar, error.function=x1$error.function, error.parameters=x1$start, error.fix=x1$error.fix)
x2<-new(Class="SaemixDiscreteOutcome", type="event")
x3<-new(Class="SaemixDiscreteOutcome", type="binary")



# Truncated gaussian distribution (verifie par rapport a definition de erf/matlab)
normcdf<-function(x,mu=0,sigma=1)
  cutoff(pnorm(-(x-mu)/sigma,lower.tail=FALSE),1e-30)


