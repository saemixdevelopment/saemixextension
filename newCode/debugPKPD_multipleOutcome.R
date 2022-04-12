saemixDir<-"/home/eco/work/saemix/saemixextension"
source(file.path(saemixDir,"R","aaa_generics.R"))

########################################################### Data
# Creating an SaemixData object with PK/PD data by hand - Responses
source(file.path(saemixDir,"Rext","SaemixData.R"))
source(file.path(saemixDir,"Rext","SaemixOutcome.R"))
pkpd<-read.table(file.path(saemixDir,"data40","warfarinPKPD.tab"), header=T)
pkpd$lwt<-log(pkpd$wt/60)

x<-new(Class="SaemixData", name.data="pkpd")
x1<-continuousOutcome(model="combined2", start=c(1, 0.2))
x2<-continuousOutcome(start=c(2))
#new(Class="SaemixOutcome",type="continuous")
x@outcome<-list(conc=new(Class="SaemixOutcome",type="continuous"), effect=new(Class="SaemixOutcome",type="continuous"))

zesuj<-unique(pkpd$id)
ind1<-match(pkpd$id, zesuj)
x@data<-data.frame(index=ind1, pkpd[,c(1:4)], ytype=pkpd$dvid, cens=0, mdv=0, pkpd[,c(9,7,8)])
x@name.group<-"id"
x@name.predictors<-c("time","amt")
x@name.X<-"time"
x@name.response<-c("dv")
x@name.ytype<-"ytype"
x@name.covariates<-colnames(pkpd)[c(9,7,8)]
x@name.cens<-"cens"
x@name.mdv<-"mdv"
x@N<-length(zesuj)
nind.obs<-tapply(pkpd$id, pkpd$id, length)
x@nind.obs<-c(nind.obs[match(zesuj,names(nind.obs))])
x@ntot.obs<-sum(x@nind.obs)
x@units<-c("hr","mg")

pkpd.saemix<-x

########################################################### Model
# Beginning of the model class - SaemixStructuralModel
source(file.path(saemixDir,"Rext","SaemixModel.R"))

# Default options
source(file.path(saemixDir,"Rext","saemixControl.R"))

# Model outcomes
xout1<-new(Class="SaemixContinuousOutcome", error.model=x1$error.model, error.npar=x1$error.npar, error.function=x1$error.function, error.parameters=x1$start, error.fix=x1$error.fix)
xout2<-new(Class="SaemixContinuousOutcome", error.model=x2$error.model, error.npar=x2$error.npar, error.function=x2$error.function, error.parameters=x2$start, error.fix=x2$error.fix)
xout<-list(conc=xout1, effect=xout2)

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
#pkpd.model<-new(Class="SaemixStructuralModel", model=model1cptdirect, outcome=list(conc=continuousOutcome(model="combined2", start=c(1, 0.2)), effect=continuousOutcome(start=c(2))))

pkpd.model<-new(Class="SaemixModel", model=model1cptdirect, outcome=pkpd.outcome, parameter=pkpd.psi0)

########################################################### Computational functions needed in e-step
source(file.path(saemixDir,"Rext","func_aux.R"))

################################### Setting up the models
## 
saemix.data<-pkpd.saemix
saemix.model<-pkpd.model

## Starting values for parameters
pkpd.psi0 <- c(ka=1, vd=5, cl=0.1, ic50=5)
nb.parameters<-length(pkpd.psi0)

# Distribution - PK log-normal, PD normal
transform.par<-rep(1,length(pkpd.psi0))
transform.par[4]<-0

# Covariance model: all parameters have IIV, correlation between CL and V
omega.init<-covariance.model<-diag(4)
covariance.model[2,3]<-covariance.model[3,2]<-1
covariance.model.fix<-omega.init*0

## Covariate model - here assume lwt on CL and V and sex on ka, beta_V,lwt=1 fixed, the other two estimated
covariate.model<-matrix(c(0,1,1,0,1,0,0,0,0,0,0,0), byrow=T, ncol=4)
covariate.model.fix<-matrix(c(0,1,0,0,0,0,0,0,0,0,0,0), byrow=T, ncol=4)
covariate.model.start<-c(1,.75,0.5)

################################### Parameter model
betaest.model<-matrix(c(rep(1,nb.parameters), c(t(covariate.model))),ncol=nb.parameters,byrow=TRUE)
colnames(betaest.model)<-names(pkpd.psi0)
rownames(betaest.model)<-rep("",dim(betaest.model)[1])
rownames(betaest.model)[1]<-"mu.psi"

################################### Starting parameters
fixedpsi.ini<-pkpd.psi0
betaI.ini<-transpsi(matrix(fixedpsi.ini,nrow=1),transform.par) #initial fixed effects (Gaussian parametrization)
fixed.ini<-betaest.model*0
fixed.ini[1,]<-betaI.ini

beta2<-t(betaest.model)
beta2[t(betaest.model)==1]<-c(betaI.ini,covariate.model.start)
betaest.model.start<-t(beta2)

################################### Covariance model
# i0.omega2 => index.omega.novar
# indest.omega => index.omega.mat
# i1.omega2 => index.omega.estim
index.omega.novar<-which((1-mydiag(covariance.model))>0) # index of parameters without IIV
index.omega.estim<-which(mydiag(covariance.model)>0)  # index of parameters with IIV
index.omega.mat<-which(covariance.model>0)
index.omega.fix<-which(covariance.model.fix>0)

# i0.omega2<-which((1-mydiag(saemix.model["covariance.model"]))>0) # index of parameters without IIV
# indest.omega<-which(saemix.model["covariance.model"]>0)
# i1.omega2<-saemix.model@indx.omega # index of parameters with IIV
# ind.res<-saemix.model["indx.res"]

################################### Design matrix
icol<-pkpd.saemix@name.covariates[which(rowSums(covariate.model)>0)]
Mcovariates<-cbind(id=1:pkpd.saemix@N, pkpd.saemix@data[!duplicated(pkpd.saemix@data[,pkpd.saemix@name.group]),icol])

#### Eco 12/04/2022 - stopped here

################################### Initialise phi

# Initialisation of phiM
if(length(i0.omega2)>0) {
  xmat<-covariate.estim[,i0.omega2]
  if(is.null(dim(xmat))) xmat<-matrix(xmat,ncol=length(i0.omega2))
  i0.temp<-which(colSums(xmat)==0)
  ind0.eta<-i0.omega2[i0.temp] # ind0.eta: index of parameters without IIV
} else ind0.eta<-c()
if(length(ind0.eta)>0) { # ind.eta: index of parameters with IIV
  idx<-1:nb.parameters
  ind.eta<-idx[-c(ind0.eta)]
} else ind.eta<-1:nb.parameters
nb.etas<-length(ind.eta)

itest.phi<-1:NM
ltest.phi<-length(itest.phi)
phiM<-matrix(data=0,nrow=NM,ncol=nb.parameters)
etaM<-matrix(data=0,nrow=NM,ncol=nb.etas)

mean.phiM<-do.call(rbind,rep(list(mean.phi),saemix.options$nb.chains))
kt<-0
omega<-saemix.model["omega.init"]
chol.omega<-try(chol(omega[ind.eta,ind.eta]),silent=TRUE)
if(inherits(chol.omega,"try-error")) {
  #	cat("ind.eta=",ind.eta,"\n")
  #	print(saemix.model["omega.init"])
  #	print(omega[ind.eta,ind.eta])
  chol.omega<-saemix.model["omega.init"][ind.eta,ind.eta]<-omega[ind.eta, ind.eta]<-mydiag(nrow=length(ind.eta),ncol=length(ind.eta))
  if(saemix.options$warnings) cat("Problem inverting covariance matrix, setting initial Omega to diagonal.\n")
}

# Find a valid set of parameters wrt to the structural.model.
# Any parameter set that does not generate NaN, inf or imaginary numbers
# will satisfy this criteria.
phiMc<-mean.phiM
while (ltest.phi>0) {
  kt<-kt+1
  if (kt==100) 
    stop("stats:fit.saemix:FailedInitialParameterGuess\nFailed to find a valid initial parameter guess\n")
  end   
  etaMc<-0.5*matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
  phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
  etaM[itest.phi,]<-etaMc[itest.phi,]
  phiM[itest.phi,]<-phiMc[itest.phi,]
  psiM<-transphi(phiM,saemix.model["transform.par"])
  fpred<-structural.model(psiM, IdM, XM)
  inan<-(is.na(fpred)+is.infinite(fpred)+(Im(fpred)!=0))
  itest.phi<-unique(IdM[inan])
  ltest.phi<-length(itest.phi)
}



################################### Creating args, Dargs, varList from initialise

# using 2 Markov chains
nchains<-2
chdat<-new(Class="SaemixRepData",data=saemix.data, nb.chains=nchains)
NM<-chdat["NM"]
IdM<-chdat["dataM"]$IdM
yM<-chdat["dataM"]$yM
XM<-chdat["dataM"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
io<-matrix(data=0,nrow=saemix.data["N"],ncol=max(saemix.data["nind.obs"]))
for(i in 1:saemix.data@N)
  io[i,1:saemix.data["nind.obs"][i]]<-1
ioM<-do.call(rbind,rep(list(io),nchains))
ind.ioM <- which(t(ioM)!=0)
DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])
phiM<-psiM<-do.call(rbind, rep(list(pkpd.psi0), saemix.data@N*nchains)) # h=id() so phiM=psiM

Dargs<-list(transform.par=c(0,0,0,0), IdM=IdM, XM=XM, yM=yM, ind.ioM=ind.ioM, model=pkpd.model@model, outcome=pkpd.model@outcome)
pres<-list(pkpd.model@outcome[[1]]@error.parameters, pkpd.model@outcome[[2]]@error.parameters)


###################################
# Data - passed on to functions, unchanged
Dargs<-list(IdM=IdM, XM=XM, yM=yM, NM=NM, N=N, ind.ioM=ind.ioM, 
            yobs=saemix.data["data"][,saemix.data["name.response"]], nobs=saemix.data["ntot.obs"],
            transform.par=rep(0, length(pkpd.psi0)),
#            transform.par=saemix.model["transform.par"], 
            model=saemix.model["structural.model"], outcome=saemix.model["outcome"],
            is.lpdfM=is.lpdfM)

flag.fmin<-0

# List of indices and variables (fixed) - passed on to functions, unchanged
nb.parest<-sum(covariate.estim)+ sum(saemix.model["covariance.model"][upper.tri(saemix.model["covariance.model"], diag=TRUE)])+1+as.integer(saemix.model["error.model"]=="combined")

Uargs<-list(nchains=saemix.options$nb.chains,nb.parameters=nb.parameters, nb.betas=nb.betas, nb.etas=nb.etas, 
            nb.parest=nb.parest,indx.betaC=indx.betaC, indx.betaI=indx.betaI, ind.res=ind.res,indest.omega=indest.omega,
            i0.omega2=i0.omega2, i1.omega2=i1.omega2,j.covariate=j.covariate, j0.covariate=j0.covariate,
            ind.fix10=ind.fix10, ind.fix11=ind.fix11, ind.fix1=ind.fix1, ind.fix0=ind.fix0,
            MCOV0=MCOV0, COV=COV, COV0=COV0, COV1=COV1, LCOV=LCOV, COV2=COV2, dstatCOV=dstatCOV, 
            Mcovariates=Mcovariates, ind.ioM=ind.ioM)
# Variability-related elements
omega.eta<-omega[ind.eta,ind.eta] # IIV matrix for estimated parameters
pres<-vector(length=length(saemix.model["outcome"]), mode="list")
for(iout in 1:length(pres))
  if(saemix.model@outcome@type=="continuous") pres[[iout]]<-saemix.model@outcome[[iout]]@error.parameters
varList<-list(pres=pres,ind0.eta<-ind0.eta
varList$ind.eta<-ind.eta
varList$omega<-omega
varList$MCOV=MCOV
varList$domega2<-do.call(cbind,rep(list((sqrt(mydiag(omega.eta)))*saemix.options$rw.ini),nb.etas))
varList$diag.omega<-mydiag(omega)

# List of options and settings (fixed) - passed on to functions, unchanged
saemix.options<-saemixControl()
saemix.options$nbiter.tot<-sum(saemix.options$nbiter.saemix)
stepsize<-rep(1,saemix.options$nbiter.tot)
stepsize[(saemix.options$nbiter.saemix[1]+1):saemix.options$nbiter.tot]<-1/
  (1:saemix.options$nbiter.saemix[2])
stepsize[1:saemix.options$nbiter.burn]<-0
opt<-list(stepsize.rw=saemix.options$stepsize.rw,stepsize=stepsize,
          proba.mcmc=saemix.options$proba.mcmc,nbiter.mcmc=saemix.options$nbiter.mcmc,
          nbiter.sa=saemix.options$nbiter.sa,nbiter.map=saemix.options$nbiter.map,alpha1.sa=saemix.options$alpha.sa,
          alpha0.sa=10^(-3/saemix.options$nbiter.sa),nbiter.saemix=saemix.options$nbiter.saemix,
          maxim.maxiter=saemix.options$maxim.maxiter,flag.fmin=flag.fmin)

################################### estep

kiter<-1

estep(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
  # E-step - simulate unknown parameters
  # Input: kiter, Uargs, mean.phi (unchanged)
  # Output: varList, DYF, phiM (changed)
  
