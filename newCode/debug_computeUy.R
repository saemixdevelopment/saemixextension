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
pkdat<-read.table(file.path(saemixDir,"data","theo.saemix.tab"), header=T)
pkpd<-read.table(file.path(saemixDir,"data40","warfarinPKPD.tab"), header=T)

x<-new(Class="SaemixData", name.data="pkdat")
x1<-continuousOutcome(model="combined2", start=c(1, 0.2))
x@outcome<-list(conc=new(Class="SaemixOutcome",type="continuous"))
colnames(pkdat)<-c("id","dose","time","conc","wt","sex")
pkdat$age<-pkdat$wt+5-2*pkdat$sex

zesuj<-unique(pkdat$id)
ind1<-match(pkdat$id, zesuj)
x@data<-data.frame(index=ind1, pkdat[,c(1,3,2,4)], ytype=1, cens=0, mdv=0, pkdat[,c(5:7)])
x@name.group<-"id"
x@name.predictors<-c("time","dose")
x@name.X<-"time"
x@name.response<-c("conc")
x@name.ytype<-"ytype"
x@name.covariates<-c("sex","wt","age")
x@name.cens<-"cens"
x@name.mdv<-"mdv"
x@N<-length(zesuj)
nind.obs<-tapply(pkdat$id, pkdat$id, length)
x@nind.obs<-c(nind.obs[match(zesuj,names(nind.obs))])
x@ntot.obs<-sum(x@nind.obs)
x@units<-c("hr","mg")
x@nvarlevel<-as.integer(1)
x@nrep.unit<-list(level1=c(nind.obs[match(zesuj,names(nind.obs))]))
x@id.unit<-as.matrix(data.frame(level1=x@data$index))

saemix.data<-x

########################################################### SaemixModel
# Model class
source(file.path(progDirExt,"SaemixVarLevel.R"))
source(file.path(progDirExt,"SaemixParameter.R"))
source(file.path(progDirExt,"SaemixParameter-methods.R"))
source(file.path(progDirExt,"SaemixIndividualModel.R"))
source(file.path(progDirExt,"SaemixModel.R"))
source(file.path(progDirExt,"SaemixModel-methods.R"))

# List of test models
source(file.path(saemixDir,"testecoExt","modelsTested.R"))
# Model 2
imodel <- 2
saemix.model <- chooseDebugModel(model=imodel) # PK model 2 (without covariates, no IIV on ka, ka estimated)
psi0<-saemix.model@mu.start

# Model 7
imodel <- 7
saemix.model <- chooseDebugModel(model=imodel) # PK model 7 (with covariates and fixed effects but no fixed covariate effects, ka fixed, no IIV on ka)
psi0<-saemix.model@mu.start

# Model 5
imodel <- 5
saemix.model <- chooseDebugModel(model=imodel) # PK model 5 (with covariates and fixed effects but no fixed covariate effects)
psi0<-saemix.model@mu.start

# Model 8
imodel <- 8
saemix.model <- chooseDebugModel(model=imodel) # PK model 7 (with covariates and fixed effects, no IIV on ka, fixed and non-fixed covariate effects on ka)
psi0<-saemix.model@mu.start

# Model 9
imodel <- 9
saemix.model <- chooseDebugModel(model=imodel) # PK model 9 (with covariates and fixed effects but no fixed covariate effects, no IIV on ka)
psi0<-saemix.model@mu.start

########################################################### saemix options
# Options and auxiliary files
source(file.path(saemixDir,"R","func_aux.R")) # error model for the old version of compute.Uy

# Options and auxiliary files - overhaul
source(file.path(saemixDir,"Rext","func_aux.R"))
source(file.path(saemixDir,"Rext","saemixControl.R"))

saemix.options<-saemixControl(nb.chains=2)
saemix.options$nbiter.tot<-sum(saemix.options$nbiter.saemix)

########################################################### Old initialisation (to check)
source(file.path(progDirExt,"main_initialiseMainAlgo_old.R"))
mylist <- initialiseMainAlgo.old(saemix.data, saemix.model, saemix.options)
Uargs1<-mylist$Uargs
Dargs1<-mylist$Dargs
varList1<-mylist$varList
opt<-mylist$opt

########################################################### Initialisation - Variability variables
# match variability levels - WIP
# TODO: test with 2 or more levels of variability
namvarlev<-c(saemix.model@var.model[[1]]@variable)
xlabel<-data.frame(level1=saemix.data@data[,saemix.data@name.group]) # labels
xindex<-data.frame(level1=saemix.data@data[,"index"])
if(saemix.model@nvarlevel==1) {
  saemix.model@nphirep<-list(rep(1,saemix.data@N))
  saemix.data@nrep.unit<-list(c(tapply(xindex[,1], xindex[,1], length)))
} else {
  saemix.model@nphirep<-saemix.data@nrep.unit<-vector(mode="list",length=saemix.model@nvarlevel)
  for(ivarlev in 2:saemix.model@nvarlevel) {  # creating labels 
    namvarlev<-c(namvarlev, saemix.model@var.model[[ivarlev]]@variable)
    xlabel<- cbind(xlabel,
                   paste(xlabel[,(ivarlev-1)],saemix.data@data[,saemix.model@var.model[[ivarlev]]@variable],sep="#"))
  }
  colnames(xlabel)<-paste0("level",1:saemix.model@nvarlevel)
  for(ivarlev in 2:saemix.model@nvarlevel) {
    xindex[,ivarlev]<-match(xlabel[,ivarlev], unique(xlabel[,ivarlev]))
  }
  colnames(xindex)<-paste0("level",1:saemix.model@nvarlevel)
  saemix.data@id.unit <- xindex
  xlabel<-xlabel[!duplicated(xlabel[,saemix.model@nvarlevel]),]
  for(ivarlev in 1:saemix.model@nvarlevel) {
    saemix.model@nphirep[[ivarlev]]<-c(tapply(xlabel[,ivarlev], xlabel[,ivarlev],length))
    saemix.data@nrep.unit[[ivarlev]]<-c(tapply(xindex[,ivarlev], xindex[,ivarlev], length))
  }
}
if(saemix.model@nvarlevel==1) {
  expect_equal(length(saemix.model@nphirep[[1]]), saemix.data@N)
}
expect_identical(xindex[,1], saemix.data@data$index)

# added xindex to the model object, can be used to obtain model predictions for different levels of variability included

########################################################### Check omega matrices and compute inverse
# Can't do it before because we need other parts of the model (namely transform) ? or maybe during the extraction
# try to move this to the extraction of the variance model
nvarlevel<-length(saemix.model@var.model)
# mean.phi<-NULL
for(ivarlev in 1:nvarlevel) {
  var.model <- saemix.model@var.model[[ivarlev]]
  omega<-indiv.model@omega
  idzero <- which(mydiag(omega)==0)
  if(length(idzero)>0) {
    for(idz in idzero) {
      omega[idz,idz]<-ifelse(saemix.model@transform.par[idz]==0, saemix.model@mu.fix[idz]*0.5, 1)
    }
    var.model@omega<-omega
  }
  chol.omega<-try(chol(omega[var.model@index.eta,var.model@index.eta]),silent=TRUE)
  if(inherits(chol.omega,"try-error")) {
    #	cat("ind.eta=",ind.eta,"\n")
    #	print(saemix.model["omega.init"])
    #	print(omega[ind.eta,ind.eta])
    chol.omega<-var.model@omega[var.model@index.eta,var.model@index.eta]<-omega[var.model@index.eta,var.model@index.eta]<-mydiag(nrow=length(var.model@index.eta),ncol=length(var.model@index.eta))
    if(verbose) message("Problem inverting covariance matrix, setting initial Omega to diagonal.\n")
  }
  var.model@chol.omega<-chol.omega
  saemix.model@var.model[[ivarlev]] <- var.model
  # mean.phi.var<-indiv.model@COV %*% indiv.model@MCOV
  # mean.phi<-cbind(mean.phi,repmatrix(mean.phi.var,saemix.model@nphirep[[ivarlev]]))
}

# mean.phiM<-do.call(rbind,rep(list(mean.phi),nb.chains))

########################################################### New initialisation - WIP
source(file.path(saemixDir,"newCode","extractDesignMatrix.R"))

saemix.model@ind.model <- vector(mode="list", length=saemix.model@nvarlevel)
for(ivarlev in 1:saemix.model@nvarlevel) {
  newlist <- extractModelStructure(saemix.data, saemix.model, level=ivarlev, verbose=FALSE)
  saemix.model@ind.model[[ivarlev]] <- newlist$indiv.model
}
Uargs2 <- newlist$Uargs # ...

# Compare arguments
if(Uargs1$nb.betas != Uargs2$nb.betas) cat("Difference in nb.betas\n")
if(newlist$flag.fmin != opt$flag.fmin) cat("Difference in flag.fmin, ",opt$flag.fmin," previously versus ",newlist$flag.fmin,"\n")
for(inam in names(Uargs2)[names(Uargs2)!='nb.betas']) {
  if(is.null(Uargs1[[inam]])) cat("Element", inam, "not in Uargs1 \n") else {
    if(!identical(unname(Uargs1[[inam]]), unname(Uargs2[[inam]]))) {
#    if(length(which (which (Uargs1[[inam]] == Uargs2[[inam]]) == FALSE))>0)
       cat("Difference in ",inam,"\n")
      print(Uargs1[[inam]])
      print(Uargs2[[inam]])
    }
  }
}
if(imodel==2) cat("Model 2: difference in Mcovariates but seems identical except row names\n")
if(imodel==7) cat("Model 7: No differences expected\n")
if(imodel==5) cat("Model 5: Expected differences for ind.fix1, ind.fix0, ind.fix11 (2 fixed covariate effects = 4, 7), COV1 (missing wt column) and dstatCOV\n")

########################################################### Initialise phiM
# mean.phi<- indmodel@COV %*% indmodel@MCOV
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

########################################################### Sum_ijk (log(p)) - compute.Uy
# Input
## b0: parameters we want to optimise on
### when more than one level of variability, contains the successive elements in the different variability levels
## phiM: matrix of individual parameters phi for each level (phi_i then phi_ik then ..., nb.modpar columns each)
## indiv.model: list of individual statistical models
## Dargs: list, use model, transform.par, XM, yM, IdM, ind.ioM, outcome, nb.chains
## DYF: a vector, modified locally but not exported
# Output
## Sum_ijk (log(p(y_ijk / psi_i, psi_ik, theta))) 

compute.Uy<-function(b0,phiM,indiv.model,Dargs,DYF) {
  # Caution, DYF is a local variable not modified, only used here for its structure
#  ivarlev<-1
  i1<-1
  for(ivarlev in 1:length(indiv.model)) {
    if(length(indiv.model[[ivarlev]]@idxmat.mcov.fixedpar.optim)>0) {
      indiv.model[[ivarlev]]@MCOV0[indiv.model[[ivarlev]]@idxmat.mcov.fixedpar.optim]<-b0[i1:(i1+length(indiv.model[[ivarlev]]@idxmat.mcov.fixedpar.optim)-1)]
      i1<-i1+length(indiv.model[[ivarlev]]@idxmat.mcov.fixedpar.optim)
      phi0 <- indiv.model[[ivarlev]]@COV0 %*% indiv.model[[ivarlev]]@MCOV0
      phiM[,(indiv.model[[1]]@nb.modpar*(ivarlev-1)+indiv.model[[ivarlev]]@index.omega.novar)]<-do.call(rbind,rep(list(phi0),Dargs$nb.chains))
    }
  }
  if(length(indiv.model)>1) {
    phiMik <- phiM[,1:indiv.model[[1]]@nb.modpar]
    for(ivarlev in 2:length(indiv.model)) {
      i1<-(ivarlev-1)*indiv.model[[1]]@nb.modpar
      phiMik <- phiMik + phiM[,(1+i1):(i1+indiv.model[[1]]@nb.modpar)]
    }
  } else phiMik<-phiM
  psiM<-transphi(phiMik,Dargs$transform.par)
  fpred<-Dargs$model(psiM,Dargs$IdM,Dargs$XM)
  lpred<-fpred
  for(iout in 1:length(Dargs$outcome)) {
    idx1<-which(Dargs$XM$ytype==iout)
    if(Dargs$outcome[[iout]]@type=="continuous") {
      if(Dargs$outcome[[iout]]@error.model=="exponential") fpred[idx1]<-log(cutoff(fpred[idx1]))
      gpred<-Dargs$outcome[[iout]]@error.function(fpred[idx1], Dargs$outcome[[iout]]@error.parameters)
      lpred[idx1]<-0.5*((Dargs$yM[idx1]-fpred[idx1])/gpred[idx1])**2+log(gpred[idx1])
    } else {
      lpred[idx1]<-(-fpred[idx1])
    }
  }
  DYF[Dargs$ind.ioM]<- lpred
  U<-sum(DYF)
  print(U)
  return(U)
}

nb.b0<-length(saemix.model@ind.model[[1]]@idxmat.mcov.fixedpar.optim)
if(nb.b0>0) {
  beta0<-rep(2,nb.b0)
  x1<-compute.Uy(beta0, phiM, saemix.model@ind.model, Dargs, DYF)
  beta0<-saemix.model@ind.model[[1]]@fixedpar[saemix.model@ind.model[[1]]@idxmat.mcov.fixedpar.optim]
  x2<-compute.Uy(beta0, phiM, saemix.model@ind.model, Dargs, DYF)
  cat("x1=",x1," , x2=",x2,"\n")
}
# Trying to optimise

beta0<-rep(2,nb.b0)
beta1<-optim(par=beta0,fn=compute.Uy,phiM=phiM,indiv.model=saemix.model@ind.model,Dargs=Dargs,DYF=DYF,control=list(maxit=100))$par

if(nb.b0>0) {
  beta0<-rep(2,nb.b0)
  x1<-compute.Uy.old(beta0, phiM, pres=varList1$pres, args=Uargs1, Dargs=Dargs1 ,DYF)
  beta0<-saemix.model@ind.model[[1]]@fixedpar[saemix.model@ind.model[[1]]@idxmat.mcov.fixedpar.optim]
  x2<-compute.Uy.old(beta0, phiM, pres=varList1$pres, args=Uargs1, Dargs=Dargs1 ,DYF)
  cat("x1=",x1," , x2=",x2,"\n")
}

source(file.path(progDirExt,"main_initialiseMainAlgo_old.R"))
mylist <- initialiseMainAlgo.old(saemix.data, saemix.model, saemix.options)
Uargs1<-mylist$Uargs
Dargs1<-mylist$Dargs
varList1<-mylist$varList
opt<-mylist$opt

# compute.Uy.old
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
    #    for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
    DYF[args$ind.ioM]<- -fpred
  }
  U<-sum(DYF)
  return(U)
}
########################################################### 



# Parameters: phiM
phiM<-mylist$phiM
psiM<-transphi(phiM, saemix.model@transform.par)
phi<-array(data=0,dim=c(Dargs1$N, Uargs1$nb.parameters, saemix.options$nb.chains))

# need suffStat$statphi1
suffStat<-list(statphi1=0,statphi2=0,statphi3=0,statrese=0)
fpred<-saemix.model@model(psiM, Dargs1$IdM, Dargs1$XM)
# for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
ff<-matrix(fpred,nrow=Dargs1$nobs,ncol=Uargs1$nchains)
for(k in 1:Uargs1$nchains) phi[,,k]<-phiM[((k-1)*Dargs1$N+1):(k*Dargs1$N),]
stat1<-apply(phi[,varList1$ind.eta,,drop=FALSE],c(1,2),sum) # sum on columns ind.eta of phi, across 3rd dimension
suffStat$statphi1 <- stat1/Uargs1$nchains
# suffStat$statphi1<-suffStat$statphi1+opt$stepsize[kiter]*(stat1/Uargs$nchains-suffStat$statphi1)

# need comega, d1.omega
domega<-cutoff(mydiag(varList1$omega[varList1$ind.eta,varList1$ind.eta]),.Machine$double.eps)
omega.eta<-varList1$omega[varList1$ind.eta,varList1$ind.eta,drop=FALSE]
omega.eta<-omega.eta-mydiag(mydiag(varList1$omega[varList1$ind.eta,varList1$ind.eta]))+mydiag(domega)
chol.omega<-try(chol(omega.eta))
d1.omega<-Uargs1$LCOV[,varList1$ind.eta]%*%solve(omega.eta)
d2.omega<-d1.omega%*%t(Uargs1$LCOV[,varList1$ind.eta])
comega<-Uargs1$COV2*d2.omega

# need betas
betas<-mylist$betas

# Old version
temp<-d1.omega[Uargs1$ind.fix11,]*(t(Uargs1$COV1)%*%(suffStat$statphi1-Uargs1$dstatCOV[,varList1$ind.eta]))

# New version: recompute dstatCOV every time
dstatCOV<-Uargs1$COV[,index.fixedpar.estim,drop=FALSE]%*% varList1$MCOV[index.fixedpar.estim,]
temp<-d1.omega[Uargs1$ind.fix11,]*(t(Uargs1$COV1)%*%(suffStat$statphi1-dstatCOV[,varList1$ind.eta]))

print(betas[Uargs1$ind.fix11])
betas[Uargs1$ind.fix11]<-solve(comega[Uargs1$ind.fix11,Uargs1$ind.fix11],rowSums(temp)) 
print(betas[Uargs1$ind.fix11])
# ECO TODO: utiliser optimise dans le cas de la dimension 1
if(length(Uargs1$ind.fix10)>0) {
  beta0<-optim(par=betas[Uargs1$ind.fix10],fn=compute.Uy.old,phiM=phiM,pres=varList1$pres,args=Uargs1,Dargs=Dargs1,DYF=DYF1,control=list(maxit=opt$maxim.maxiter))$par 
} else cat("Nothing to optimise, ind.fix10 empty\n")

if(FALSE) {
  dim(d1.omega[Uargs1$ind.fix11,])
  dim(t(Uargs1$COV1))
  dim(Uargs1$dstatCOV[,varList1$ind.eta])
  dim((suffStat$statphi1))
}
