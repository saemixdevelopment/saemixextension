# WIP...
# library(numDeriv)
# library(pracma)

# Computation of the FIM
computeFIM.lin<-function(saemixObject) {
  # Estimate the Fisher Information Matrix 
  # Input: saemix object from a fit
  # Output: updated saemix object
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
  
  linLL <- saemix.LLlin(theta, modeltype=saemix.model@modeltype, etype=saemix.data["data"][["ytype"]], f0=f0, DF0=DF, z=z, Mcov=Mcov, nfix=length(saemix.res@betas), nomega=dim(omega)[1], nind.obs=saemix.data["nind.obs"], ndat.exp=ndat.exp, idx.fix=saemix.res["indx.fix"],idx.omega=idx.omega, idx.res=saemix.res["indx.res"], tr.fix=saemix.model["transform.par"]) 
  
  # Hessian (needs pracma)
  # hessLL<-try(pracma::hessian(saemix.LLlin, theta, modeltype=saemix.model@modeltype, etype=saemix.data["data"][["ytype"]], f0=f0, DF0=DF, z=z, Mcov=Mcov, nfix=length(saemix.res@betas), nomega=dim(omega)[1], nind.obs=saemix.data["nind.obs"], ndat.exp=ndat.exp, idx.fix=saemix.res["indx.fix"], idx.omega=idx.omega, idx.res=saemix.res["indx.res"], tr.fix=saemix.model["transform.par"]))
  # Hessian (needs numDeriv)
  hessLL<-try(numDeriv::hessian(saemix.LLlin, theta, modeltype=saemix.model@modeltype, etype=saemix.data["data"][["ytype"]], f0=f0, DF0=DF, z=z, Mcov=Mcov, nfix=length(saemix.res@betas), nomega=dim(omega)[1], nind.obs=saemix.data["nind.obs"], ndat.exp=ndat.exp, idx.fix=saemix.res["indx.fix"], idx.omega=idx.omega, idx.res=saemix.res["indx.res"], tr.fix=saemix.model["transform.par"]))
  
  
  if(inherits(hessLL,"try-error")) { # try solving per block
    if(saemixObject@options$warnings) cat("Error computing the hessian of the log-likelihood.\n")
    invisible(saemixObject)
  }
  
  # ECO ici modifie car role de covariate.estim pas clair
  # covariate.estim=si un parametre (et ses covariables associees) sont estimees ou non
  #  covariate.estim<-matrix(rep(saemix.model["fixed.estim"], dim(saemix.model["betaest.model"])[1]),byrow=TRUE, ncol=length(saemix.model["fixed.estim"]))*saemix.model["betaest.model"] # 29/05/20
  
  # Parameter names
  covariate.estim<-saemix.model["betaest.model"]
  covariate.estim[1,]<-saemix.model["fixed.estim"]
  
  j<-which(saemix.model["betaest.model"]>0)
  ind.fixed.est<-(covariate.estim[j]>0)
  npar<-sum(ind.fixed.est)
  # Tracking indices for covariances
  myidx.omega<-c()
  myidx.cor<-c()
  name.rand1<-name.rand2<-c()
  myidx.track<-NULL
  ipar<-npar
  for(iom in 1:dim(covariance.model)[1]) {
    for(jom in iom:dim(covariance.model)[1]) {
      if(covariance.model[iom,jom]==1) {
        ipar<-ipar+1
        myidx.track<-rbind(myidx.track,c(ipar,iom,jom))
        if(iom==jom) {
          myidx.omega<-c(myidx.omega,ipar)
          name.rand1<-c(name.rand1,paste("Var",saemixObject@model@name.modpar[iom],sep="."))
          name.rand2<-c(name.rand2,paste("SD",saemixObject@model@name.modpar[iom],sep="."))
        }
        else {
          myidx.cor<-c(myidx.cor,ipar)
          name.rand1<-c(name.rand1,paste("Cov",saemixObject@model@name.modpar[iom],saemixObject@model@name.modpar[jom],sep="."))
          name.rand2<-c(name.rand2,paste("Corr",saemixObject@model@name.modpar[iom],saemixObject@model@name.modpar[jom],sep="."))
        }
      }
    }
  }
  if(length(myidx.cor)>0) {
    track.var<-myidx.track[myidx.track[,1] %in% myidx.omega,]
    for(i in myidx.cor) {
      ij<-which(myidx.track[,1]==i)
      myidx.track[ij,2]<-track.var[track.var[,2]==myidx.track[ij,2],1]
      myidx.track[ij,3]<-track.var[track.var[,2]==myidx.track[ij,3],1]
    }
  }
  if(saemixObject@model@modeltype=="structural")
    namallpar<-c(saemixObject@results@name.fixed,name.rand1, saemixObject@results@name.sigma[saemixObject@results@indx.res], name.rand2) else
      namallpar<-c(saemixObject@results@name.fixed,name.rand1, name.rand2)
  
  invMF <- try(solve(-hessLL))
  if(inherits(invMF,"try-error")) { # try solving per block
    invMF<-hessLL*0
    Fth<--hessLL[c(1:npar),c(1:npar)]
    Cth<-try(solve(Fth))
    FO<--hessLL[-c(1:npar),-c(1:npar)]
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
  
  sTHest<-sqrt(mydiag(Cth))
  #sTH<-matrix(0,1,length(saemix.res["betas"]))
  sTH<-rep(0,length(saemix.res["betas"]))
  sTH[ind.fixed.est]<-sTHest
  se.fixed<-sTH
  
  sO<-sqrt(mydiag(CO))
  se.omega<-matrix(0,nphi,1)
  se.sdcor<-se.cov<-matrix(0,nphi,nphi)
  se.omega[saemix.model["indx.omega"]]<-sO[myidx.omega-npar]
  se.res<-matrix(0,2*nytype,1)
  if(saemixObject@model@modeltype=="structural") se.res[saemix.res["indx.res"]]<-sO[(nomega+1):length(sO)]    
  
  # Table with SE, CV and confidence intervals
  estpar<-c(saemixObject@results@fixed.effects)
  estSE<-c(se.fixed)
  est1<-est2<-se1<-se2<-c()
  if(length(myidx.cor)>0) {
    ipar<-npar
    for(iom in 1:nphi) {
      for(jom in iom:nphi) {
        if(covariance.model[iom,jom]==1) {
          ipar<-ipar+1
          se.cov[iom,jom]<-se.cov[jom,iom]<-sO[(ipar-npar)]
          est1<-c(est1,omega[iom,jom])
          se1<-c(se1,sO[ipar-npar])
          if(iom==jom) {
            se.sdcor[iom,iom]<-sqrt(CO[iom,iom])/2/sqrt(omega[iom,iom])
            est2<-c(est2,sqrt(omega[iom,iom]))
            se2<-c(se2,se.sdcor[iom,iom])
          } else { # compute correlation and SE on correlation using the delta-method
            ebet<-c(omega[iom,jom],omega[iom,iom],omega[jom,jom])
            varbet<-CO[(myidx.track[myidx.track[,1]==ipar,]-npar),(myidx.track[myidx.track[,1]==ipar,]-npar)]
            rho<-ebet[1]/sqrt(ebet[2]*ebet[3])
            debet<-c(1/sqrt(ebet[2]*ebet[3]), -ebet[1]/(ebet[2]**(3/2))/sqrt(ebet[3])/2, -ebet[1]/(ebet[3]**(3/2))/sqrt(ebet[2])/2)
            se.sdcor[iom,jom]<-se.sdcor[jom,iom]<-t(debet) %*% varbet %*% debet
            est2<-c(est2,rho)
            se2<-c(se2,se.sdcor[iom,jom])
          }
        }
      }
    }
    if(saemixObject@model@modeltype=="structural") estpar<-c(estpar,est1,saemixObject@results@respar[saemixObject@results@indx.res],est2) else
      estpar<-c(estpar,est1,est2)
    if(saemixObject@model@modeltype=="structural") estSE<-c(estSE,se1,se.res[saemixObject@results@indx.res],se2) else estSE<-c(estSE,se1,se2)
  } else {
    diag(se.cov)<-se.omega
    if(saemixObject@model@modeltype=="structural")  
      estpar<-c(estpar,diag(omega)[saemixObject@results@indx.omega],saemixObject@results@respar[saemixObject@results@indx.res], sqrt(diag(omega)[saemixObject@results@indx.omega])) else estpar<-c(estpar,diag(omega)[saemixObject@results@indx.omega], sqrt(diag(omega)[saemixObject@results@indx.omega])) 
      if(saemixObject@model@modeltype=="structural")
        estSE<-c(estSE,se.omega[saemixObject@results@indx.omega],se.res[saemixObject@results@indx.res],se.omega[saemixObject@results@indx.omega]/2/sqrt(diag(omega)[saemixObject@results@indx.omega])) else estSE<-c(estSE,se.omega[saemixObject@results@indx.omega],se.omega[saemixObject@results@indx.omega]/2/sqrt(diag(omega)[saemixObject@results@indx.omega]))
  }
  conf.int<-data.frame(name=namallpar, estimate=estpar, se=estSE)
  conf.int$cv<-100*conf.int$se/conf.int$estimate
  conf.int$lower<-conf.int$estimate - 1.96*conf.int$se
  conf.int$upper<-conf.int$estimate + 1.96*conf.int$se
  saemix.res["se.fixed"]<-se.fixed
  saemix.res["se.omega"]<-c(se.omega)
  saemix.res["se.cov"]<-se.cov
  if(saemixObject@model@modeltype=="structural") saemix.res["se.respar"]<-c(se.res)
  saemix.res["conf.int"]<-conf.int
  saemix.res["ll.lin"]<-c(linLL)
  saemix.res["fim"]<-(-hessLL)
  saemix.res["aic.lin"]<-(-2)*saemix.res["ll.lin"]+ 2*saemix.res["npar.est"]
  saemix.res["bic.lin"]<-(-2)*saemix.res["ll.lin"]+ log(saemix.data["N"])*saemix.res["npar.est"]
  saemix.res["bic.covariate.lin"]<-(-2)*saemix.res["ll.lin"]+ log(saemix.data["N"])*saemix.res["nbeta.random"]+log(sum(saemix.data["nind.obs"]))*saemix.res["nbeta.fixed"]
  
  ##################################
  #delete(hw)
  saemixObject["results"]<-saemix.res
  return(saemixObject)
}


# inverse FIM as hessian of the log-likelihood - on the phi (eta) scale
computeFIMinv.lin<-function(saemixObject) {
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
  
  # Hessian (needs pracma)
  # hessLL<-try(pracma::hessian(saemix.LLlin, theta, modeltype=saemix.model@modeltype, etype=saemix.data["data"][["ytype"]], f0=f0, DF0=DF, z=z, Mcov=Mcov, nfix=length(saemix.res@betas), nomega=dim(omega)[1], nind.obs=saemix.data["nind.obs"], ndat.exp=ndat.exp, idx.fix=saemix.res["indx.fix"], idx.omega=idx.omega, idx.res=saemix.res["indx.res"], tr.fix=saemix.model["transform.par"]))
  # Hessian (needs numDeriv)
  hessLL<-try(numDeriv::hessian(saemix.LLlin, theta, modeltype=saemix.model@modeltype, etype=saemix.data["data"][["ytype"]], f0=f0, DF0=DF, z=z, Mcov=Mcov, nfix=length(saemix.res@betas), nomega=dim(omega)[1], nind.obs=saemix.data["nind.obs"], ndat.exp=ndat.exp, idx.fix=saemix.res["indx.fix"], idx.omega=idx.omega, idx.res=saemix.res["indx.res"], tr.fix=saemix.model["transform.par"]))
  
  if(inherits(hessLL,"try-error")) { # try solving per block
    if(saemixObject@options$warnings) cat("Error computing the hessian of the log-likelihood.\n")
    return()
  }
  
  # ECO ici modifie car role de covariate.estim pas clair
  # covariate.estim=si un parametre (et ses covariables associees) sont estimees ou non
  #  covariate.estim<-matrix(rep(saemix.model["fixed.estim"], dim(saemix.model["betaest.model"])[1]),byrow=TRUE, ncol=length(saemix.model["fixed.estim"]))*saemix.model["betaest.model"] # 29/05/20
  
  npar<-length(saemix.res@betas)
  invMF <- try(solve(-hessLL))
  if(inherits(invMF,"try-error")) { # try solving per block
    invMF<-hessLL*0
    Fth<--hessLL[c(1:npar),c(1:npar)]
    Cth<-try(solve(Fth))
    FO<--hessLL[-c(1:npar),-c(1:npar)]
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
  # provisionally return inverse of MF
  return((invMF))
}

# Linearised log-likelihood
computeLL.lin<-function(saemixObject) {
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
  
  LLlin <- saemix.LLlin(theta, modeltype=saemix.model@modeltype, etype=saemix.data["data"][["ytype"]], f0=f0, DF0=DF, z=z, Mcov=Mcov, 
                        nfix=length(saemix.res@betas), nomega=dim(omega)[1], nind.obs=saemix.data["nind.obs"], ndat.exp=ndat.exp, 
                        idx.fix=saemix.res["indx.fix"],idx.omega=idx.omega, idx.res=saemix.res["indx.res"], tr.fix=saemix.model["transform.par"]) 
  
  return(LLlin)
}

# Linearised log-likelihood for a vector of parameters theta - on psi scale
saemix.LLlin <- function(theta, modeltype="structural", etype, f0, DF0, z, Mcov, nfix, nomega, nind.obs, ndat.exp=0, idx.fix, idx.omega, idx.res, tr.fix) {
  # Input
  ## theta: population parameters as a vector, in the following order
  ### fixed effects (mu and betas)
  ### omega
  ### residual error parameters (saemix.res@respar,saemix.data["data"][["ytype"]]) if modeltype is structural
  ## f0: f(ti, phi_i) where phi_i=conditional mean estimates of the individual parameters
  ## DF0: df(ti, phi_i)
  ## z: vector of values y_i-f(ti,phi_i)+df(ti,phi_i).phi_i
  ## Mcov: Mcovariate matrices for each subject
  ## nfix: number of mu+beta on phi scale (to multiply with design matrix for each subject)
  ## nomega: number of random effects
  ## nind.obs: vector with the number of observations for each subject
  ## ndat.exp: number of observations for which an exponential error model is used
  ## idx.fix: indices of which parameters are mu
  ## idx.omega: indices of lower triangular elements to replace in matrix omega 
  ## idx.res: indices of residual error parameters 
  ## tr.fix: h transformation (from psi to phi)
  # Output
  ## linearised log-likelihood
  ntheta<-length(theta)
  betas <- theta[1:nfix]
  betas[idx.fix] <- transpsi(t(betas[idx.fix]), tr.fix)
  nsuj <- length(nind.obs)
  omega <- matrix(0, nrow=nomega, ncol=nomega)
  # g0<-cutoff(saemix.res["respar"][1]+saemix.res["respar"][2]*abs(f0))
  if (modeltype=="structural"){
    respar<-c(0,0)
    respar[idx.res]<-theta[(ntheta-(length(idx.res))+1):ntheta]
    g0<-error(f0, respar, etype) 
    omega[lower.tri(omega, diag=T)][idx.omega] <- theta[(nfix+1):(ntheta-(length(idx.res)))]
  }
  else 
    omega[lower.tri(omega, diag=T)][idx.omega] <- theta[(nfix+1):ntheta]
  omega<-omega+t(omega)-diag(diag(omega)) # reconstructing symmetrical variance matrix
  
  invVi<-Gi<-list() # Individual variance matrices
  j2<-0
  for (i in 1:nsuj) {
    ni<-nind.obs[i]
    j1<-j2+1
    j2<-j2+ni
    #    z[j1:j2]<-yobs[j1:j2] - f0[j1:j2] + DF0[j1:j2,,drop=FALSE]%*%hat.phi[i,] # computed before the call to the function and passed in z
    if (modeltype=="structural"){
      Vi<- DF0[j1:j2,,drop=FALSE] %*% omega %*% t(DF0[j1:j2,,drop=FALSE]) + mydiag((g0[j1:j2])^2, nrow=ni)
    } else{
      Vi<- DF0[j1:j2,,drop=FALSE] %*% t(DF0[j1:j2,,drop=FALSE])+ mydiag(1, nrow=ni) # ?? should be 0 ?? anyway shouldn't be computed for non-Gaussian models; and where is omega ?
    }
    #    invVi[[i]]<-solve(Vi[[i]])
    # Invert avoiding numerical problems
    Gi[[i]]<-round(Vi*1e12)/1e12
    VD<-try(eigen(Gi[[i]]))
    if(inherits(VD,"try-error") || det(Gi[[i]])==0) {
      cat("Unable to compute the FIM by linearisation.\n") # si matrice de variance non inversible
      stop()
    }
    D<-Re(VD$values)
    V<-Re(VD$vectors)
    invVi[[i]] <- V%*%mydiag(1/D,nrow=length(D))%*%t(V)
  } 
  LLlin<- -0.5*sum(nind.obs)*log(2*pi)
  j2<-0
  for (i in 1:nsuj) {
    ni<-nind.obs[i]
    j1<-j2+1
    j2<-j2+ni
    DFi<-DF0[j1:j2,,drop=FALSE]
    zi<-z[j1:j2]
    Ai <- Mcov[[i]]
    DFAi<-DFi%*%Ai
    Dzi<-zi-DFAi %*% betas
    
    LLlin <- LLlin - 0.5*log(det(Gi[[i]])) - 0.5*t(Dzi)%*% invVi[[i]] %*%Dzi 
  }
  LLlin<-LLlin-ndat.exp
  return(LLlin)
}

# Linearised individual log-likelihoods for a vector of parameters theta - on psi scale
saemix.LLlin.ind <- function(theta, modeltype="structural", etype, f0, DF0, z, Mcov, nfix, nomega, nind.obs, idx.fix, idx.omega, idx.res, tr.fix) {
  # Input
  ## theta: population parameters as a vector, in the following order
  ### fixed effects (mu and betas)
  ### omega
  ### residual error parameters (saemix.res@respar,saemix.data["data"][["ytype"]]) if modeltype is structural
  ## etype: error model type
  ## f0: f(ti, phi_i) where phi_i=conditional mean estimates of the individual parameters
  ## DF0: df(ti, phi_i)
  ## z: vector of values y_i-f(ti,phi_i)+df(ti,phi_i).phi_i
  ## Mcov: Mcovariate matrices for each subject
  ## nfix: number of mu+beta on phi scale (to multiply with design matrix for each subject)
  ## nomega: number of random effects
  ## nind.obs: vector with the number of observations for each subject
  ## idx.fix: indices of which parameters are mu
  ## idx.omega: indices of lower triangular elements to replace in matrix omega 
  ## idx.res: indices of residual error parameters 
  ## tr.fix: h transformation (from psi to phi)
  # Output
  ## linearised log-likelihood for each subject (vector), without the constant ndat.exp
  ntheta<-length(theta)
  betas <- theta[1:nfix]
  betas[idx.fix] <- transpsi(t(betas[idx.fix]), tr.fix)
  nsuj <- length(nind.obs)
  omega <- matrix(0, nrow=nomega, ncol=nomega)
  # g0<-cutoff(saemix.res["respar"][1]+saemix.res["respar"][2]*abs(f0))
  if (modeltype=="structural"){
    respar<-c(0,0)
    respar[idx.res]<-theta[(ntheta-(length(idx.res))+1):ntheta]
    g0<-error(f0, respar, etype) 
    omega[lower.tri(omega, diag=T)][idx.omega] <- theta[(nfix+1):(ntheta-(length(idx.res)))]
  }
  else 
    omega[lower.tri(omega, diag=T)][idx.omega] <- theta[(nfix+1):ntheta]
  omega<-omega+t(omega)-diag(diag(omega)) # reconstructing symmetrical variance matrix
  
  invVi<-Gi<-list() # Individual variance matrices
  j2<-0
  LLlin <- c()
  for (i in 1:nsuj) {
    ni<-nind.obs[i]
    j1<-j2+1
    j2<-j2+ni
    #    z[j1:j2]<-yobs[j1:j2] - f0[j1:j2] + DF0[j1:j2,,drop=FALSE]%*%hat.phi[i,] # computed before the call to the function and passed in z
    if (modeltype=="structural"){
      Vi<- DF0[j1:j2,,drop=FALSE] %*% omega %*% t(DF0[j1:j2,,drop=FALSE]) + mydiag((g0[j1:j2])^2, nrow=ni)
    } else{
      Vi<- DF0[j1:j2,,drop=FALSE] %*% t(DF0[j1:j2,,drop=FALSE])+ mydiag(1, nrow=ni) # ?? should be 0 ?? anyway shouldn't be computed for non-Gaussian models; and where is omega ?
    }
    #    invVi[[i]]<-solve(Vi[[i]])
    # Invert avoiding numerical problems
    Gi[[i]]<-round(Vi*1e12)/1e12
    VD<-try(eigen(Gi[[i]]))
    if(inherits(VD,"try-error") || det(Gi[[i]])==0) {
      cat("Unable to compute the FIM by linearisation.\n") # si matrice de variance non inversible
      stop()
    }
    D<-Re(VD$values)
    V<-Re(VD$vectors)
    invVi[[i]] <- V%*%mydiag(1/D,nrow=length(D))%*%t(V)
  } 
  LLlin<- -0.5*sum(nind.obs)*log(2*pi)
  j2<-0
  for (i in 1:nsuj) {
    ni<-nind.obs[i]
    j1<-j2+1
    j2<-j2+ni
    DFi<-DF0[j1:j2,,drop=FALSE]
    zi<-z[j1:j2]
    Ai <- Mcov[[i]]
    DFAi<-DFi%*%Ai
    Dzi<-zi-DFAi %*% betas
    
    LLlin <- c(LLlin, - 0.5*log(det(Gi[[i]])) - 0.5*t(Dzi)%*% invVi[[i]] %*%Dzi)
  }
# log-likelihood: sum(LLlin)-ndat.exp [where ndat.exp=nb of observations exponentially transformed]
  return(LLlin)
}

