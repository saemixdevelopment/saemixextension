# Gradient of the likelihood
computeJacFIMinv.lin<-function(saemixObject) {
  # Estimate the Fisher Information Matrix 
  # Input: saemix object from a fit
  # Output: FIM (currently) => TODO: updated saemix object
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  xind<-saemix.data["data"][,saemix.data["name.predictors"],drop=FALSE]
  yobs<-saemix.data["data"][,saemix.data["name.response"]]
  
  #  covariance.model<-0*saemix.model["covariance.model"]
  covariance.model<-saemix.model["covariance.model"]
  omega<-saemix.res["omega"]
  omega.null<-0*omega
  #  diag(covariance.model)<-mydiag(saemix.model["covariance.model"])
  #  omega<-0*saemix.res["omega"] # Why use only diag(omega) ???
  #  diag(omega)<-mydiag(saemix.res["omega"])
  hat.phi<-saemix.res["cond.mean.phi"]
  nphi<-dim(hat.phi)[2]
  nomega<-sum(covariance.model[lower.tri(covariance.model,diag=TRUE)])
  if (saemixObject["model"]["modeltype"]=="structural"){
    nres<-length(saemix.res["indx.res"])
  } else{
    nres <- 0
  }
  nytype<-length(unique(saemix.data["data"]["ytype"]))
  dphi<-cutoff(abs(colMeans(hat.phi))*1e-4,1e-10)
  coefphi<-c(0,-1,1)
  
  F<-array(data=0,dim=c(saemix.data["ntot.obs"],nphi,length(coefphi)))
  gs<-matrix(0,saemix.data["ntot.obs"],4)
  etype.exp<-which(saemix.model["error.model"]=='exponential')
  
  for (l in 1:length(coefphi)) {
    for (j in 1:nphi) {
      phi<-hat.phi
      phi[,j]<-phi[,j]+coefphi[l]*dphi[j]
      psi<-transphi(phi,saemix.model["transform.par"])
      f <- saemix.model["model"](psi, saemix.data["data"][,"index"],xind)
      for(ityp in etype.exp) f[saemix.data["data"][,saemix.data["name.ytype"]]==ityp]<-log(cutoff(f[saemix.data["data"][,saemix.data["name.ytype"]]==ityp]))    
      F[,j,l]<-f
    }
  }
  #  DF<-(F[,,3]-F[,,2])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi))/2 
  DF<-(F[,,3]-F[,,1])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi)) #gradient of f (changed from F[,,2] to F[,,1])
  z<-matrix(0,saemix.data["ntot.obs"],1)
  
  ind.covariates<-which(saemix.model["betaest.model"]>0)
  f0<-F[,1,1]
  
  # Covariance model
  omega <- saemix.res["omega"]
  idx.omega<-which(omega[lower.tri(omega, diag=T)]>0)
  
  # Design matrices
  Mcov <- vector(mode='list', length=saemix.data["N"])
  j2<-0
  for (i in 1:saemix.data["N"]) {
    ni<-saemix.data["nind.obs"][i]
    j1<-j2+1
    j2<-j2+ni
    z[j1:j2]<-yobs[j1:j2] - f0[j1:j2] + DF[j1:j2,,drop=FALSE]%*%hat.phi[i,]
    Ai<-kronecker(diag(nphi),as.matrix(saemix.model["Mcovariates"][i,]))
    Mcov[[i]]<-Ai[,ind.covariates,drop=FALSE]
  }
  ndat.exp<-0
  for(ityp in etype.exp) ndat.exp<-ndat.exp+sum(yobs[saemix.data["data"][,saemix.data["name.ytype"]]==ityp])
  
  theta<-c(saemix.res@fixed.effects, omega[lower.tri(omega,diag=T)][idx.omega])
  if(saemix.model@modeltype=="structural") theta<-c(theta, saemix.res["respar"][saemix.res["indx.res"]])
  
  #  saemix.LLlin(theta, modeltype=saemix.model@modeltype, etype=saemix.data["data"][["ytype"]], f0=f0, DF0=DF, z=z, Mcov=Mcov, nfix=length(saemix.res@betas), nomega=dim(omega)[1], nind.obs=saemix.data["nind.obs"], ndat.exp=ndat.exp, idx.fix=saemix.res["indx.fix"],idx.omega=idx.omega, idx.res=saemix.res["indx.res"], tr.fix=saemix.model["transform.par"]) 
  
  # Gradient (needs numDeriv)
  jacobLL<-try(numDeriv::jacobian(saemix.LLlin, theta, modeltype=saemix.model@modeltype, etype=saemix.data["data"][["ytype"]], f0=f0, DF0=DF, z=z, Mcov=Mcov, nfix=length(saemix.res@betas), nomega=dim(omega)[1], nind.obs=saemix.data["nind.obs"], ndat.exp=ndat.exp, idx.fix=saemix.res["indx.fix"], idx.omega=idx.omega, idx.res=saemix.res["indx.res"], tr.fix=saemix.model["transform.par"], method.args=list(d=1e-6)))
  
  if(inherits(jacobLL,"try-error")) {
    if(saemixObject@options$warnings) cat("Error computing the gradient of the log-likelihood.\n")
    invisible(saemixObject)
  }
  jacMF <- t(jacobLL) %*% jacobLL
  invMF <- try(solve(jacMF))
#  print(invMF)
  
  # ECO ici modifie car role de covariate.estim pas clair
  # covariate.estim=si un parametre (et ses covariables associees) sont estimees ou non
  #  covariate.estim<-matrix(rep(saemix.model["fixed.estim"], dim(saemix.model["betaest.model"])[1]),byrow=TRUE, ncol=length(saemix.model["fixed.estim"]))*saemix.model["betaest.model"] # 29/05/20
  if(FALSE) {
    if(inherits(invMF,"try-error")) { # try solving per block
      invMF<-jacMF*0
      Fth<--jacMF[c(1:npar),c(1:npar)]
      Cth<-try(solve(Fth))
      FO<--jacMF[-c(1:npar),-c(1:npar)]
      CO<-try(solve(FO))
      if(!inherits(Cth,"try-error")) invMF[1:npar, 1:npar]<-Cth
      if(!inherits(CO,"try-error")) invMF[-c(1:npar),-c(1:npar)]<-CO
      if(inherits(Cth,"try-error") |inherits(CO,"try-error")) {
        successFIM <- FALSE
        if(saemixObject@options$warnings) cat("Error computing the Fisher Information Matrix: singular system.\n")
      } else {
        successFIM <- TRUE
        if(saemixObject@options$warnings) cat("Full FIM cannot be inverted, separating blocs for fixed effects and variance parameters.\n")
      }
    } else {
      successFIM <- TRUE
      Cth<-invMF[c(1:npar),c(1:npar)] # bloc fixed effect
      CO<-invMF[-c(1:npar),-c(1:npar)] # bloc random effect
      Ccor<-invMF[-c(1:npar),c(1:npar)] # correlation bloc
    }
  }
  # provisionally return inverse of MF
  return((jacobLL))
}

initialiseComputationLL <- function(saemixObject) {
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  xind<-saemix.data["data"][,saemix.data["name.predictors"],drop=FALSE]
  yobs<-saemix.data["data"][,saemix.data["name.response"]]
  
  #  covariance.model<-0*saemix.model["covariance.model"]
  covariance.model<-saemix.model["covariance.model"]
  omega<-saemix.res["omega"]
  omega.null<-0*omega
  #  diag(covariance.model)<-mydiag(saemix.model["covariance.model"])
  #  omega<-0*saemix.res["omega"] # Why use only diag(omega) ???
  #  diag(omega)<-mydiag(saemix.res["omega"])
  hat.phi<-saemix.res["cond.mean.phi"]
  nphi<-dim(hat.phi)[2]
  nomega<-sum(covariance.model[lower.tri(covariance.model,diag=TRUE)])
  if (saemixObject["model"]["modeltype"]=="structural"){
    nres<-length(saemix.res["indx.res"])
  } else{
    nres <- 0
  }
  nytype<-length(unique(saemix.data["data"]["ytype"]))
  dphi<-cutoff(abs(colMeans(hat.phi))*1e-4,1e-10)
  coefphi<-c(0,-1,1)
  
  F<-array(data=0,dim=c(saemix.data["ntot.obs"],nphi,length(coefphi)))
  gs<-matrix(0,saemix.data["ntot.obs"],4)
  etype.exp<-which(saemix.model["error.model"]=='exponential')
  
  for (l in 1:length(coefphi)) {
    for (j in 1:nphi) {
      phi<-hat.phi
      phi[,j]<-phi[,j]+coefphi[l]*dphi[j]
      psi<-transphi(phi,saemix.model["transform.par"])
      f <- saemix.model["model"](psi, saemix.data["data"][,"index"],xind)
      for(ityp in etype.exp) f[saemix.data["data"][,saemix.data["name.ytype"]]==ityp]<-log(cutoff(f[saemix.data["data"][,saemix.data["name.ytype"]]==ityp]))    
      F[,j,l]<-f
    }
  }
  #  DF<-(F[,,3]-F[,,2])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi))/2 
  DF<-(F[,,3]-F[,,1])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi)) #gradient of f (changed from F[,,2] to F[,,1])
  z<-matrix(0,saemix.data["ntot.obs"],1)
  
  ind.covariates<-which(saemix.model["betaest.model"]>0)
  f0<-F[,1,1]
  
  # Covariance model
  omega <- saemix.res["omega"]
  idx.omega<-which(omega[lower.tri(omega, diag=T)]>0)
  
  # Design matrices
  Mcov <- vector(mode='list', length=saemix.data["N"])
  j2<-0
  for (i in 1:saemix.data["N"]) {
    ni<-saemix.data["nind.obs"][i]
    j1<-j2+1
    j2<-j2+ni
    z[j1:j2]<-yobs[j1:j2] - f0[j1:j2] + DF[j1:j2,,drop=FALSE]%*%hat.phi[i,]
    Ai<-kronecker(diag(nphi),as.matrix(saemix.model["Mcovariates"][i,]))
    Mcov[[i]]<-Ai[,ind.covariates,drop=FALSE]
  }
  ndat.exp<-0
  for(ityp in etype.exp) ndat.exp<-ndat.exp+sum(yobs[saemix.data["data"][,saemix.data["name.ytype"]]==ityp])
  
  theta<-c(saemix.res@fixed.effects, omega[lower.tri(omega,diag=T)][idx.omega])
  if(saemix.model@modeltype=="structural") theta<-c(theta, saemix.res["respar"][saemix.res["indx.res"]])
  
  return(list(theta=theta, modeltype=saemix.model@modeltype, etype=saemix.data["data"][["ytype"]], f0=f0, DF0=DF, z=z, Mcov=Mcov, nfix=length(saemix.res@betas), nomega=dim(omega)[1], nind.obs=saemix.data["nind.obs"], ndat.exp=ndat.exp, idx.fix=saemix.res["indx.fix"], idx.omega=idx.omega, idx.res=saemix.res["indx.res"], tr.fix=saemix.model["transform.par"]))
}
