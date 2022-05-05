# Data and model definition from newCode/stepMultiple.R
saemix.model<-pkpd.model
saemix.data<-pkpd.saemix

# Initialise
# List of options and settings (fixed) - passed on to functions, unchanged
saemix.options<-saemixControl()
saemix.options$nbiter.tot<-sum(saemix.options$nbiter.saemix)
stepsize<-rep(1,saemix.options$nbiter.tot)
stepsize[(saemix.options$nbiter.saemix[1]+1):saemix.options$nbiter.tot]<-1/
  (1:saemix.options$nbiter.saemix[2])
stepsize[1:saemix.options$nbiter.burn]<-0

# TODO +++
## check indices
## check flag.fmin

### flag.fmin: 1 if at least one parameter has IIV (???)
#### how to generalise to more than one level, should flag.fmin become a vector ?
ivarlevel<-1
flag.fmin<-as.integer(sum(diag(saemix.model@var.model[[ivarlevel]]@omega.model *(1-saemix.model@var.model[[ivarlevel]]@omega.model.fix)))>0)

opt<-list(stepsize.rw=saemix.options$stepsize.rw,stepsize=stepsize,
          proba.mcmc=saemix.options$proba.mcmc,nbiter.mcmc=saemix.options$nbiter.mcmc,
          nbiter.sa=saemix.options$nbiter.sa,nbiter.map=saemix.options$nbiter.map,alpha1.sa=saemix.options$alpha.sa,
          alpha0.sa=10^(-3/saemix.options$nbiter.sa),nbiter.saemix=saemix.options$nbiter.saemix,
          maxim.maxiter=saemix.options$maxim.maxiter,flag.fmin=flag.fmin)
                  
# E-step
index.eta.fix<-index.eta<-domega2<-diag.omega<-omega<-list()
for(i in length(saemix.model@var.model)) {
  index.eta[[i]]<-saemix.model@var.model[[i]]@index.eta
  index.eta.fix[[i]]<-which(diag(saemix.model@var.model[[i]]@omega.model.fix)==1)
  nb.etas<-length(index.eta[[i]])
  omega[[i]]<-saemix.model@var.model[[i]]@omega
  diag.omega[[i]]<-mydiag(omegalist[[i]])
  domega2[[i]]<-do.call(cbind,rep(list((sqrt(mydiag(omega.eta)))*saemix.options$rw.ini),nb.etas))
}

#   varList<-list(pres=pres,ind0.eta=ind0.eta,ind.eta=ind.eta,omega=omega, MCOV=MCOV,
#                 domega2=do.call(cbind,rep(list((sqrt(mydiag(omega.eta)))*saemix.options$rw.ini),nb.etas)),diag.omega=mydiag(omega))

varList<-list(index.eta=index.eta, index.eta.fix=index.eta.fix, omega=omega, diag.omega=diag.omega, domega2=domega2)

# Dargs
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
Dargs<-list(transform.par=saemix.model@transform.par, IdM=IdM, XM=XM, yM=yM, ind.ioM=ind.ioM, NM=NM,
            model=saemix.model@model, outcome=saemix.model@outcome)

# Uargs (not sure we still need it, TBC)
## ind.ioM => dans Dargs
Uargs<-list(nchains=nchains)

# DYF
DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])

# error.parameters
error.parameters<-list()
for(iout in 1:saemix.model@nb.outcome)
  error.parameters[[iout]]<-saemix.model@outcome[[iout]]@error.parameters

# Parameters at iteration kiter
nsuj<-saemix.data@N
phiM<-mean.phi<-do.call(rbind,rep(list(saemix.model@mu.start), nsuj))
iiv.sd<-diag(varList$omega[[ivarlevel]])
for(icol in 1:dim(phiM)[2]) {
  phiM[,icol]<-phiM[,icol]*exp(rnorm(nsuj, mean=0, sd=sqrt(iiv.sd[icol])))
}
phiM<-do.call(rbind,rep(list(phiM), Uargs$nchains))

# Setup E-step
nb.etas<-length(varList$index.eta[[ivarlevel]])
omega.eta<-varList$omega[[ivarlevel]][varList$index.eta[[ivarlevel]], varList$index.eta[[ivarlevel]], drop=FALSE]
domega<-cutoff(mydiag(omega.eta),.Machine$double.eps)
omega.eta<-omega.eta-mydiag(mydiag(omega.eta))+mydiag(domega)
chol.omega<-try(chol(omega.eta))
somega<-solve(omega.eta)

VK<-rep(c(1:nb.etas),2)
mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))
# ? why this line ? (meaning of ind0.eta)
# phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]
# not equivalent (previously no index.eta.fix, so something else)
# phiM[,varList$index.eta.fix[[ivarlevel]]]<-mean.phiM[,varList$index.eta.fix[[ivarlevel]]]

U.y<-compute.LLy(phiM, Dargs, DYF, error.parameters)

etaM<-phiM[,varList$index.eta[[ivarlevel]]]-mean.phiM[,varList$index.eta[[ivarlevel]],drop=FALSE]
phiMc<-phiM
for(u in 1:opt$nbiter.mcmc[1]) { # 1er noyau
  etaMc<-matrix(rnorm(Dargs$NM*nb.etas),ncol=nb.etas)%*%chol.omega
  phiMc[,varList$index.eta[[ivarlevel]]]<-mean.phiM[,varList$index.eta[[ivarlevel]]]+etaMc
  Uc.y<-compute.LLy(phiMc, Dargs, DYF, error.parameters)
  deltau<-Uc.y-U.y
  ind<-which(deltau<(-1)*log(runif(Dargs$NM)))
  etaM[ind,]<-etaMc[ind,]
  U.y[ind]<-Uc.y[ind]
}
U.eta<-0.5*rowSums(etaM*(etaM%*%somega))


#   Dargs<-list(IdM=IdM, XM=XM, yM=yM, NM=NM, N=N, nobs=saemix.data["ntot.obs"],
#               yobs=saemix.data["data"][,saemix.data["name.response"]],transform.par=saemix.model["transform.par"],
#               error.model=saemix.model["error.model"],structural.model=structural.model, is.lpdf=is.lpdf)
#   Uargs<-list(nchains=saemix.options$nb.chains,nb.parameters=nb.parameters, nb.betas=nb.betas, nb.etas=nb.etas, 
#               nb.parest=nb.parest,indx.betaC=indx.betaC, indx.betaI=indx.betaI, ind.res=ind.res,
#               indest.omega=indest.omega, i0.omega2=i0.omega2, i1.omega2=i1.omega2,	j.covariate=j.covariate, 
#               ind.fix10=ind.fix10, ind.fix11=ind.fix11, ind.fix1=ind.fix1, ind.fix0=ind.fix0,
#               MCOV0=MCOV0, COV=COV, COV0=COV0, COV1=COV1, LCOV=LCOV, COV2=COV2, dstatCOV=dstatCOV, 
#               Mcovariates=Mcovariates, ind.ioM=ind.ioM)



compute.LLy<-function(phiM, Dargs, DYF, error.parameters) {
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$model(psiM,Dargs$IdM,Dargs$XM)
  lpred<-fpred
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  for(iout in 1:length(Dargs$outcome)) {
    idx1<-which(Dargs$XM$ytype==iout)
    # print(summary(fpred[idx1]))
    if(Dargs$outcome[[iout]]@type=="continuous") {
      if(Dargs$outcome[[iout]]@error.model=="exponential") fpred[idx1]<-log(cutoff(fpred[idx1]))
      gpred<-Dargs$outcome[[iout]]@error.function(fpred[idx1], error.parameters[[iout]])
      lpred[idx1]<-0.5*((Dargs$yM[idx1]-fpred[idx1])/gpred)**2+log(gpred)
      # print(summary(gpred))
    } else lpred[idx1]<- (-fpred[idx1])
  }
  DYF[Dargs$ind.ioM]<-lpred
  U<-colSums(DYF)
  return(U)
}


