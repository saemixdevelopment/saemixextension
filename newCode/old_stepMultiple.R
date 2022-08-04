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
# Mean value for population parameters
pkpd.psi0<-c(ka=1, vd=5, cl=0.1, ic50=5)

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



# Input
## b0: parameters we want to optimise
## args: list, use j0.covariate, COV0, MCOV0, i0.omega2, nchains, ind.ioM
## phiM: matrix of parameters phi
## pres: vector of residual error parameters
## Dargs: list, use transform.par, modeltype, structural.model(), etype.exp, XM, yM
## DYF: a vector, modified locally but not exported

compute.Uy.old<-function(b0,phiM,pres,args,Dargs,DYF) {
  # Attention, DYF variable locale non modifiee en dehors
  args$MCOV0[args$j0.covariate]<-b0
  phi0<-args$COV0 %*% args$MCOV0
  phiM[,args$i0.omega2]<-do.call(rbind,rep(list(phi0),args$nchains))
  psiM<-transphi(phiM,Dargs$transform.par)
  if (Dargs$modeltype=="structural"){
    fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
    for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
    gpred<-error(fpred,pres,Dargs$XM$ytype)
    DYF[args$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
  } else {
    fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
    for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
    DYF[args$ind.ioM]<- -fpred
  }
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


function(saemix.data, saemix.model, level=1, verbose=FALSE) {
  if(is.numeric(level)) {
    if(level<=length(names(var.model))) level<-names(saemix.model@var.model)[level] else level<-names(saemix.model@var.model)[1]
  }
  var.model<-saemix.model@var.model[[level]]
  covariate.model<-saemix.model@covariate.model
  betaest.model<-rbind(rep(1,saemix.model@npar), covariate.model)
  name.variable<-var.model@variable
  data <- saemix.data@data
  
  #  if(is.character(level)) level<-match(level, names(object@var.model))
  index.omega.novar<-which((1-mydiag(var.model["omega.model"]))>0) # index of parameters without IIV ## i0.omega2
  idxmat.omega<-which(var.model["omega.model"]>0) # index of 1's in omega.model ## indest.omega
  index.omega.var<-which(mydiag(var.model["omega.model"])>0) # index of parameters with IIV ## i1.omega2
  # Matching covariates in the model and the data
  ivar<-data[, name.variable]
  idx<-which(rownames(covariate.model) %in% saemix.data@name.covariates)
  covariate.model<-covariate.model[idx,]
  N<-length(unique(ivar))
  # Covariate model & design matrix
  ## Mcovariate: nrow=N, ncol=1 (id) + nb of covariates in the model
  ### for each subject, 1 in the ID column and the values of the covariates for that subject 
  ## TODO: here or before ? apply transformations before creating Mcovariates [see SaemixcovariateTransform.R]
  ### continuous covariates => replace with value
  ### binary covariates => transform to 0/1
  ### categorical covariates => create (ncat-1) columns
  ### also transform covariate.model (may add rows eg for 2 dummy variables) and adjust everything
  if(length(covariate.model)==0) {
    tab<-data.frame(ivar=ivar) 
    Mcovariates<-data.frame(ivar=rep(1,N))
  } else {
    tab<-data.frame(ivar=ivar,data[, rownames(covariate.model),drop=FALSE])
    temp2<-unique(tab)
    temp<-tab[!duplicated(ivar),,drop=FALSE]
    if(dim(temp)[1]!=dim(temp2)[1]) {
      if(verbose) messsage("Some covariates have time-varying values; only the first is taken into account in the current version of the algorithm.\n")
    }
    Mcovariates<-data.frame(ivar=rep(1,N),temp[,2:dim(temp)[2]])
  }
  colnames(Mcovariates)[1]<-<-name.variable

#  mu.start<-saemix.model@mu.start # initial fixed effects (original parametrization)
  muphi.start<-transpsi(matrix(saemix.model@mu.start,nrow=1),saemix.model["transform.par"]) # starting value for mu on the scale of phi (Gaussian parametrization)
  if(length(covariate.model)>0) {
    beta.start<-covariate.model
    beta.start[covariate.model==1]<- saemix.model@beta.start
    matrixfixedpar.start<-rbind(muphi.start, beta.start)
    fixedpar.start<-matrixfixedpar.start[which(betaest.model>0)]
    fixedpar.start<-matrix(fixedpar.start,ncol=1)
  } else {
    matrixfixedpar.start<-matrix(muphi.start, nrow=1)
    fixedpar.start<-t(matrixfixedpar.start)
  }
  # covariate.estim <-1-covariate.model.fix
  # betas.ini => fixedpar.start
  # ind.covariates => idxmat.fixedpar
  idxmat.fixedpar<-which(saemix.model["betaest.model"]==1)
  nb.fixedpar<-length(idxmat.fixedpar) # total nb of fixed parameters (mu and beta), was nb.betas, should be length(fixedpar.start)

  # the covariates in matrix form to compute the vector of parameters as phi^p = mu^p + sum_q beta_q^p cov_q
  ### order: mu_1 for par1, beta_11, beta_12,... for par1,  mu_2 for par2, beta_21, beta_22,... for par2, ...
  ## LCOV, MCOV: nrow=nb of mu+nb of betas; ncol=nb of parameters in the model (=modpar, eg ka, cl)
  ### LCOV: design matrix to pass from the vector of fixed parameters (mu+beta) to the model parameters 
  ### MCOV: values of fixed parameters in matrix form (same structure)
  ## COV: nrow=N, ncol=nb of mu+beta
  ### when a column corresponds to a muphi, it is filled with 1
  ### when a column corresponds to a beta, the values are those of the corresponding covariate for each subject
  ## COV2: nrow=ncol=nb of mu+beta obtained as t(COV).COV (initialised, not used here)
  ## COV1: nrow=N, ncol=nb of columns of COV associated to fixed parameters in the model (???)
  ## dstatCOV: nrow=N, ncol=nb of parameters in the model
  ### formed as the product of COV.MCOV, removing the columns corresponding to fixed parameters
  ### for each parameter, 
  nb.parameters<-saemix.model@npar
  LCOV<-matrix(data=0,nrow=nb.fixedpar,ncol=nb.parameters)
  colnames(LCOV)<-saemix.model@name.modpar
  j1<-1
  COV<-matrix(nrow=dim(Mcovariates)[1],ncol=0)
  name.fixpar <-c()
#  pfix<-matrix(data=0,nrow=1,ncol=nb.parameters)
  mean.phi<-matrix(data=0,nrow=N,ncol=nb.parameters) # population parameters for subject i = mu + sum_q beta_q cov_i,q
  for(j in 1:nb.parameters) {
    name.fixpar<-c(name.fixpar, paste0("muphi.",saemix.model@name.modpar[j]))
    jcov<-which(betaest.model[,j]==1)
    if(length(jcov)>1) name.fixpar<-c(name.fixpar, paste0("beta.",saemix.model@name.modpar[j],".",rownames(betaest.model)[jcov[jcov>1]]))
    lambdaj<-matrixfixedpar.start[jcov,j]
    aj<-as.matrix(Mcovariates[,jcov])
    COV<-cbind(COV,aj)
    nlj<-length(lambdaj)
    j2<-j1+nlj-1
    LCOV[j1:j2,j]<-matrix(data=1,nrow=nlj,ncol=1)
    j1<-j2+1
    if(length(jcov)<=1) mean.phi[,j]<-aj*lambdaj else mean.phi[,j]<-aj%*%lambdaj
#    pfix[j]<-length(lambdaj)
  }
#   pfix<-colSums(betaest.model) # ? needed ?
  # indx.betaI and indx.betaC (index of mu and of beta) renamed as idxmat.mu and idxmat.beta
  # index.fixedpar.fix which of these parameters are fixed, see if still needed [was ind.fix1]
   if(length(covariate.model)>0) {
     idxmat.mu<-(c(1:nb.parameters)-1)*dim(betaest.model)[1]+1
     idxmat.beta<-which(betaest.model==1)
     idxmat.beta<-setdiff(idxmat.beta, idxmat.mu)
   } else {
     idxmat.mu<-c(1:nb.parameters)
     idxmat.beta<-c()
   }
  rownames(LCOV)<-name.fixpar
  rownames(COV)<-1:dim(COV)[1]
  MCOV<-LCOV
  
  COV2<-t(COV)%*%COV
  # j.covariate => idxmat.mcov.par: index of elements in LCOV present in the model
  idxmat.mcov.par<-which(LCOV==1)
  MCOV[idxmat.mcov.par]<-fixedpar.start
  betas<-fixedpar.start
  
  betaest.model.fix <- rbind(saemix.model@mu.fix, saemix.model@covariate.model.fix)
  index.fixedpar.fix<-which(betaest.model.fix[idxmat.fixedpar]==1)  # fixed fixed effects (fixed elements of mu+beta)
  index.fixedpar.estim<-which(betaest.model.fix[idxmat.fixedpar]==0) # estimated fixed effects (estimated elements of mu+beta)
  COV1 <- COV[,index.fixedpar.fix,drop=FALSE]
  dstatCOV<-COV[,index.fixedpar.estim,drop=FALSE]%*%MCOV[index.fixedpar.estim,,drop=FALSE]

  
  # Index of fixed parameters (mu+beta) on parameters with IIV that are fixed to their starting value
  # ind.fix11 => index.fixedpariiv.fix (unused for the moment)
  betaest.model.fix1 <- betaest.model.fix
  if(length(index.omega.novar)>0) betaest.model.fix1[,index.omega.novar]<-0
  index.fixedpariiv.fix<-which(betaest.model.fix1[idxmat.fixedpar]==1)

  # Index of fixed parameters (mu+beta) on parameters without IIV that are fixed to their starting value
  # ind.fix10 => index.fixedparnoiiv.fix (index in 1:nb(mu+beta))
  # j0.covariate => idxmat.mcov.fixedparnoiiv.fix (index in the matrix MCOV (or LCOV))
  betaest.model.fix0 <- betaest.model.fix
  if(length(index.omega.var)>0) betaest.model.fix0[,index.omega.var]<-0
  index.fixedparnoiiv.fix<-which(betaest.model.fix0[idxmat.fixedpar]==1)
  
  MCOV0<-MCOV[index.fixedparnoiiv.fix,index.omega.novar,drop=FALSE]
  COV0<-COV[,index.fixedparnoiiv.fix, drop=FALSE]
  idxmat.mcov.fixedparnoiiv.fix <- which(LCOV[index.fixedparnoiiv.fix, index.omega.novar]==1)
  # flag.fmin is set to 1 if at least one mu in the model has no IIV and is fixed to its starting value [makes no sense]
  flag.fmin <- as.integer(sum(betaest.model.fix0[1,])>0)
  
  # WIP - 12/05/2022
  
  Uargs<-list(nchains=saemix.options$nb.chains,nb.parameters=nb.parameters, nb.betas=nb.betas, nb.etas=nb.etas, 
              nb.parest=nb.parest,indx.betaC=indx.betaC, indx.betaI=indx.betaI, ind.res=ind.res,indest.omega=indest.omega,
              i0.omega2=i0.omega2, i1.omega2=i1.omega2,j.covariate=j.covariate, j0.covariate=j0.covariate,
              ind.fix10=ind.fix10, ind.fix11=ind.fix11, ind.fix1=ind.fix1, ind.fix0=ind.fix0,
              MCOV0=MCOV0, COV=COV, COV0=COV0, COV1=COV1, LCOV=LCOV, COV2=COV2, dstatCOV=dstatCOV, 
              Mcovariates=Mcovariates, ind.ioM=ind.ioM)
  
  var.model@covariate.model <- covariate.model
  var.model@Mcovariates <- Mcovariates
  return(list(Mcovariates=Mcovariates))
}

## TODO debug and extract elements for Uargs
### indx.res, ind.res => now unused
if(FALSE) {
  
  # ECO TODO: integrate all this section in the object creation ?
  # Initialisation: 
  # create local copies modified of omega.init and covariate.model in saemix.model
  # A la fin: i1.omega2 renomme en indx.omega et ind.res en indx.res
  i0.omega2<-which((1-mydiag(saemix.model["covariance.model"]))>0) # index of parameters without IIV
  indest.omega<-which(saemix.model["covariance.model"]>0)
  #  i1.omega2<-which(mydiag(saemix.model$covariance.model)>0)
  i1.omega2<-saemix.model@indx.omega # index of parameters with IIV
  ind.res<-saemix.model["indx.res"]
  
  # Covariate model & design matrix
  id<-saemix.data["data"][,saemix.data["name.group"]]
  if(length(saemix.data["name.covariates"])==0) tab<-data.frame(id=id) else
    tab<-data.frame(id=id,saemix.data["data"][, saemix.data["name.covariates",drop=FALSE]])
  temp2<-unique(tab)
  temp<-tab[!duplicated(id),,drop=FALSE]
  if(dim(temp)[1]!=dim(temp2)[1]) {
    if(saemix.options$warnings) cat("Some covariates have time-varying values; only the first is taken into account in the current version of the algorithm.\n")
  }
  #temp<-temp[order(temp[,1]),]
  if(length(saemix.data["name.covariates"])>0) {
    Mcovariates<-data.frame(id=rep(1,N),temp[,2:dim(temp)[2]])} else {
      Mcovariates<-data.frame(id=rep(1,N))
    }
  # removing from model unused lines - section removed, now dealt with during the creation of the model
  j.cov<-which(rowSums(saemix.model["betaest.model"])>0)
  betaest.model<-saemix.model["betaest.model"][j.cov,,drop=FALSE]
  Mcovariates<-Mcovariates[,j.cov,drop=FALSE] # eliminate all the unused covariates
  for(icol in 1:dim(Mcovariates)[2])
    if(is.factor(Mcovariates[,icol])) Mcovariates[,icol]<-as.numeric(Mcovariates[,icol])-1
  saemix.model["betaest.model"]<-betaest.model
  temp1<-betaest.model[-c(1),,drop=FALSE]
  #  if(is.null(dim(temp1))) temp1<-matrix(temp1,nrow=1, dimnames=list(rownames(betaest.model)[-c(1)], colnames(betaest.model)))  
  saemix.model["covariate.model"]<-temp1
  
  fixedpsi.ini<-saemix.model["psi0"][1,] # initial fixed effects (original parametrization)
  betaI.ini<-transpsi(matrix(fixedpsi.ini,nrow=1),saemix.model["transform.par"]) #initial fixed effects (Gaussian parametrization)
  fixed.ini<-saemix.model["betaest.model"]*0
  fixed.ini[1,]<-betaI.ini
  nr.psi0<-dim(saemix.model["psi0"])[1]
  nr.cov<-dim(saemix.model["betaest.model"])[1]
  if(nr.psi0>nr.cov) {
    saemix.model["psi0"]<-saemix.model["psi0"][1:nr.cov,,drop=FALSE]
    nr.psi0<-dim(saemix.model["psi0"])[1]
  }
  #t1<-NULL
  if(nr.psi0<nr.cov) {
    #  t1<-t(covariate.model[(nr.psi0+1):nr.cov,])
    psi1<-saemix.model["psi0"][nr.psi0,]
    for(j in (nr.psi0+1):nr.cov)
      saemix.model["psi0"]<-rbind(saemix.model["psi0"],psi1)
    nr.psi0<-dim(saemix.model["psi0"])[1]
  }
  if(nr.psi0>1) fixed.ini[2:nr.psi0,]<-saemix.model["psi0"][2:nr.psi0,]
  
  covariate.estim<-saemix.model["betaest.model"]
  covariate.estim[1,]<-1-saemix.model["mu.fix"]
  
  betas.ini<-fixed.ini[which(saemix.model["betaest.model"]>0)]
  betas.ini<-matrix(betas.ini,ncol=1)
  
  nb.betas<-sum(saemix.model["betaest.model"])
  ind.covariate<-which(saemix.model["betaest.model"]==1)
  #matrix(which(covariate.model==1),nrow=1)
  
  # the covariates
  LCOV<-MCOV<-matrix(data=0,nrow=nb.betas,ncol=nb.parameters)
  j1<-1
  COV<-matrix(nrow=dim(Mcovariates)[1],ncol=0)
  pfix<-matrix(data=0,nrow=1,ncol=nb.parameters)
  mean.phi<-matrix(data=0,nrow=N,ncol=nb.parameters)
  for(j in 1:nb.parameters) {
    jcov<-which(saemix.model["betaest.model"][,j]==1)
    lambdaj<-fixed.ini[jcov,j]
    aj<-as.matrix(Mcovariates[,jcov])
    COV<-cbind(COV,aj)
    nlj<-length(lambdaj)
    j2<-j1+nlj-1
    LCOV[j1:j2,j]<-matrix(data=1,nrow=nlj,ncol=1)
    j1<-j2+1
    if(length(jcov)<=1) mean.phi[,j]<-aj*lambdaj else mean.phi[,j]<-aj%*%lambdaj
    pfix[j]<-length(lambdaj)
  }
  
  if(nb.parameters>1){
    indx.betaI<-cumsum(c(0,pfix[1:(nb.parameters-1)]))+1	
  } else{
    indx.betaI<-1
  }
  
  idx<-1:nb.betas
  indx.betaC<-idx[is.na(match(idx,indx.betaI))]
  saemix.model["indx.fix"]<-indx.betaI
  saemix.model["indx.cov"]<-indx.betaC
  
  COV2<-t(COV)%*%COV
  j.covariate<-which(LCOV==1)
  MCOV[j.covariate]<-betas.ini
  betas<-betas.ini
  
  ind.fix1<-which(covariate.estim[ind.covariate]==1)
  ind.fix0<-which(covariate.estim[ind.covariate]==0)
  COV1<-COV[,ind.fix1]
  #if(length(ind.fix0)==1) dstatCOV<-matrix(COV[,ind.fix0],ncol=1)%*%MCOV[ind.fix0,] else 
  dstatCOV<-COV[,ind.fix0,drop=FALSE]%*%MCOV[ind.fix0,]
  
  covariate.estim1<-covariate.estim
  covariate.estim1[,i0.omega2]<-0
  ind.fix11<-which(covariate.estim1[ind.covariate]==1)
  covariate.estim0<-covariate.estim
  covariate.estim0[,i1.omega2]<-0
  ind.fix10<-which(covariate.estim0[ind.covariate]==1)
  MCOV0<-MCOV[ind.fix10,i0.omega2,drop=FALSE]
  #if(is.null(dim(MCOV0)) & length(MCOV0)>0) MCOV0<-matrix(MCOV0,ncol=1)
  COV0<-COV[,ind.fix10]
  j0.covariate<-which(LCOV[ind.fix10,i0.omega2]==1)
  flag.fmin<-as.integer(sum(covariate.estim0[1,])>0)
  
}
