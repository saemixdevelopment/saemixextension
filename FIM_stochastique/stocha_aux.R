deriv_mu_omega_i <- function(i, COV, LCOV, phiM2, mu_iiv, elim,
                             covariance.model, Omega_inv, Omega_Kronecker, Omega,
                             n_mu_iiv){
  # For one subject, return derivative of individual log_likelihood in:
  #     fixed effects with iiv -> deriv_mu_i
  #     Omega's components -> deriv_omega_i
  if(n_mu_iiv>1){
    C_i.t = diag(COV[i,])%*%LCOV
  }else{
    C_i.t = diag(COV[i])%*%LCOV
  }
  
  C_i = t(C_i.t)
  phi_Cmu = phiM2[i,] - C_i %*% mu_iiv
  
  deriv_mu_i =  C_i.t %*% Omega_inv %*% phi_Cmu 
  
  deriv_omega_i = 0.5 * elim %*% Omega_Kronecker %*% as.vector( phi_Cmu%*%t(phi_Cmu) - Omega )
  
  return(list( deriv_mu_i = as.vector(deriv_mu_i),
               deriv_omega_i = as.vector(deriv_omega_i)
  ))
}

hess_mu_omega_i <- function(i, COV, LCOV, phiM2, mu_iiv, elim, elim.t,
                             covariance.model, Omega_inv, Omega_Kronecker, Omega,
                             n_mu_iiv){
  # For one subject, return second derivative of individual log_likelihood in:
  #     deriv_mu_mu_i - deriv_omega_mu_i - deriv_omega_omega_i

  if(n_mu_iiv>1){
    C_i.t = diag(COV[i,])%*%LCOV
  }else{
    C_i.t = diag(COV[i])%*%LCOV
  }
  
  C_i = t(C_i.t)
  phi_Cmu = phiM2[i,] - C_i %*% mu_iiv
  phi_Cmu.t = t(phi_Cmu)
  S_i = phi_Cmu %*% phi_Cmu.t
  
  deriv_mu_mu_i = - C_i.t %*% Omega_inv %*% C_i
  
  deriv_omega_mu_i  = - kronecker(phi_Cmu.t %*% Omega_inv, C_i.t %*% Omega_inv) %*% elim.t
  
  deriv_omega_omega_i = - elim %*% kronecker(Omega_inv, Omega_inv %*% (S_i - 0.5*Omega) %*% Omega_inv) %*% elim.t
  
  return(list( deriv_mu_mu_i = deriv_mu_mu_i,
               deriv_omega_mu_i = deriv_omega_mu_i, 
               deriv_omega_omega_i = deriv_omega_omega_i
  ))
}



# compute.LLy mais avec moyenne sur les chaines
compute.LLy.2<-function(phiM,args,Dargs,DYF,pres) {
  DYF = DYF[,1:Dargs$N]
  args$ind.ioM =  args$ind.ioM[1:length(Dargs$yobs)]
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM[1:length(Dargs$yobs)],Dargs$XM[1:length(Dargs$yobs),])
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  if (Dargs$modeltype=="structural"){
    gpred<-error(fpred,pres,Dargs$XM$ytype)
    DYF[args$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
  } else {
    DYF[args$ind.ioM]<- -fpred
    # cat("compute.LLy.2 \n")
  }
  U<-colSums(DYF)
  return(U)
}

transfFim.sa <-function(saemixObject, FIM.stocha) {
  # Transform the Fisher Information Matrix
  # and compute the s.e. of the estimated parameters  

  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]

  covariance.model<-saemix.model["covariance.model"]
  omega<-saemix.res["omega"]
  omega.null<-0*omega

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
 
  ind.covariates<-which(saemix.model["betaest.model"]>0)

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
  
  MF = FIM.stocha
  
  if (sum(ind.fixed.est)>0) {
    Mparam<-matrix(0,dim(saemix.model["betaest.model"])[1], dim(saemix.model["betaest.model"])[2])
    Mparam[1,]<-saemix.model["transform.par"]
    Mtp<-Mparam[saemix.model["betaest.model"]>0]    
    Mtp<-Mtp[ind.fixed.est]
    dbetas <- dtransphi(saemix.res["betas"][ind.fixed.est],Mtp)
    
    J_diag = rep(1, times = ncol(MF))
    J_diag[which(ind.fixed.est)] = 1/dbetas
    J = diag(J_diag)
    MF = t(J)%*%MF%*%J
  } 
  
  fim <- varCovMat(MF) 
  saemix.res["fim.sa"]<-fim
  
  inv_MF<-try(as.inverse(fim)) #try(solve(fim))
  if(inherits(inv_MF,"try-error")) {
    if(saemixObject@options$warnings) cat("Error computing the Fisher Information Matrix: singular system.\n")
    inv_MF<-NA*fim
  }
  
  sTHest<-sqrt(mydiag(inv_MF[c(1:npar),c(1:npar)]))
  sTH<-rep(0,length(saemix.res["betas"]))
  sTH[which(ind.fixed.est)]<-sTHest
  se.fixed<-sTH

  sO<-sqrt(mydiag(inv_MF[-c(1:npar),-c(1:npar)]))
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
            se.sdcor[iom,iom]<-se.cov[iom,iom]/2/sqrt(omega[iom,iom])
            est2<-c(est2,sqrt(omega[iom,iom]))
            se2<-c(se2,se.sdcor[iom,iom])
          } else { # compute correlation and SE on correlation using the delta-method
            ebet<-c(omega[iom,jom],omega[iom,iom],omega[jom,jom])
            varbet<-inv_MF[myidx.track[myidx.track[,1]==ipar,],myidx.track[myidx.track[,1]==ipar,]]
            rho<-ebet[1]/sqrt(ebet[2]*ebet[3])
            debet<-c(1/sqrt(ebet[2]*ebet[3]), -ebet[1]/(ebet[2]**(3/2))/sqrt(ebet[3])/2, -ebet[1]/(ebet[3]**(3/2))/sqrt(ebet[2])/2)
            se.sdcor[iom,jom]<-se.sdcor[jom,iom]<-sqrt( t(debet) %*% varbet %*% debet )
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
  saemix.res["se.fixed.sa"]<-se.fixed
  saemix.res["se.omega.sa"]<-c(se.omega)
  saemix.res["se.cov.sa"]<-se.cov
  if(saemixObject@model@modeltype=="structural") saemix.res["se.respar.sa"]<-c(se.res)
  saemix.res["conf.int.sa"]<-conf.int
 
  saemixObject["results"]<-saemix.res
  return(saemixObject)
  #  return(list(ll.lin,fim,DFi, Dzi, invVi))
}

transfFim.sa_hess <-function(saemixObject, FIM.stocha) {
  # Transform the Fisher Information Matrix
  # and compute the s.e. of the estimated parameters  
  
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  
  covariance.model<-saemix.model["covariance.model"]
  omega<-saemix.res["omega"]
  omega.null<-0*omega
  
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
  
  ind.covariates<-which(saemix.model["betaest.model"]>0)
  
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
  
  MF = FIM.stocha
  
  if (sum(ind.fixed.est)>0) {
    Mparam<-matrix(0,dim(saemix.model["betaest.model"])[1], dim(saemix.model["betaest.model"])[2])
    Mparam[1,]<-saemix.model["transform.par"]
    Mtp<-Mparam[saemix.model["betaest.model"]>0]    
    Mtp<-Mtp[ind.fixed.est]
    dbetas <- dtransphi(saemix.res["betas"][ind.fixed.est],Mtp)
    
    J_diag = rep(1, times = ncol(MF))
    J_diag[which(ind.fixed.est)] = 1/dbetas
    J = diag(J_diag)
    MF = t(J)%*%MF%*%J
  } 
  
  fim <- varCovMat(MF) 
  saemix.res["fim.sa_hess"]<-fim
  
  inv_MF<-try(as.inverse(fim)) #try(solve(fim))
  if(inherits(inv_MF,"try-error")) {
    if(saemixObject@options$warnings) cat("Error computing the Fisher Information Matrix: singular system.\n")
    inv_MF<-NA*fim
  }
  
  sTHest<-sqrt(mydiag(inv_MF[c(1:npar),c(1:npar)]))
  sTH<-rep(0,length(saemix.res["betas"]))
  sTH[which(ind.fixed.est)]<-sTHest
  se.fixed<-sTH
  
  sO<-sqrt(mydiag(inv_MF[-c(1:npar),-c(1:npar)]))
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
            se.sdcor[iom,iom]<-se.cov[iom,iom]/2/sqrt(omega[iom,iom])
            est2<-c(est2,sqrt(omega[iom,iom]))
            se2<-c(se2,se.sdcor[iom,iom])
          } else { # compute correlation and SE on correlation using the delta-method
            ebet<-c(omega[iom,jom],omega[iom,iom],omega[jom,jom])
            varbet<-inv_MF[myidx.track[myidx.track[,1]==ipar,],myidx.track[myidx.track[,1]==ipar,]]
            rho<-ebet[1]/sqrt(ebet[2]*ebet[3])
            debet<-c(1/sqrt(ebet[2]*ebet[3]), -ebet[1]/(ebet[2]**(3/2))/sqrt(ebet[3])/2, -ebet[1]/(ebet[3]**(3/2))/sqrt(ebet[2])/2)
            se.sdcor[iom,jom]<-se.sdcor[jom,iom]<-sqrt( t(debet) %*% varbet %*% debet )
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
  saemix.res["se.fixed.sa_hess"]<-se.fixed
  saemix.res["se.omega.sa_hess"]<-c(se.omega)
  saemix.res["se.cov.sa_hess"]<-se.cov
  if(saemixObject@model@modeltype=="structural") saemix.res["se.respar.sa_hess"]<-c(se.res)
  saemix.res["conf.int.sa_hess"]<-conf.int
  
  saemixObject["results"]<-saemix.res
  return(saemixObject)
  #  return(list(ll.lin,fim,DFi, Dzi, invVi))
}