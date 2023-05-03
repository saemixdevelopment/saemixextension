################## Stochastic approximation - compute sufficient statistics (M-step) #####################
mstep.multi<-function(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat, deltai) {
  # M-step - stochastic approximation
  # Input: kiter, Uargs, structural.model, DYF, phiM (unchanged)
  # Output: varList, phi, betas, suffStat (changed)
  #					mean.phi (created)
  
  # Update variances - TODO - check if here or elsewhere
  nb.etas<-length(varList$ind.eta)
  domega<-cutoff(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]),.Machine$double.eps)
  omega.eta<-varList$omega[varList$ind.eta,varList$ind.eta,drop=FALSE]
  omega.eta<-omega.eta-mydiag(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]))+mydiag(domega)
  #  print(varList$omega.eta)
  chol.omega<-try(chol(omega.eta))
  d1.omega<-Uargs$LCOV[,varList$ind.eta]%*%solve(omega.eta)
  d2.omega<-d1.omega%*%t(Uargs$LCOV[,varList$ind.eta])
  comega<-Uargs$COV2*d2.omega
  
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)
  for(itype in Dargs$etype.exp) fpred[Dargs$XM$ytype==itype]<-log(cutoff(fpred[Dargs$XM$ytype==itype]))
  #	if(Dargs$error.model=="exponential")
  #		fpred<-log(cutoff(fpred))
  ff<-matrix(fpred,nrow=Dargs$nobs,ncol=Uargs$nchains)
  for(k in 1:Uargs$nchains) phi[,,k]<-phiM[((k-1)*Dargs$N+1):(k*Dargs$N),]
  # overall speed similar
  #    phi<-aperm(array(phiM,c(N,nchains,3)),c(1,3,2))
  stat1<-apply(phi[,varList$ind.eta,,drop=FALSE],c(1,2),sum) # sum on columns ind.eta of phi, across 3rd dimension
  stat2<-matrix(data=0,nrow=nb.etas,ncol=nb.etas)
  stat3<-apply(phi**2,c(1,2),sum) #  sum on phi**2, across 3rd dimension
  if(length(suffStat)>3) statr<-rep(0,length(suffStat)-3)
  for(k in 1:Uargs$nchains) {
    tab_ytype = Dargs$XM$ytype[1:length(Dargs$yobs)]   ### pour gérer le nombre de chaînes et que ytype soit de la meme longueur que yobs et fk 
    phik<-phi[,varList$ind.eta,k]
    stat2<-stat2+t(phik)%*%phik
    fk<-ff[,k]
    if(length(suffStat)>3) { # at least one response is not given by its LL
      for(itype in 1:length(Dargs$modeltype)) {
        if(Dargs$modeltype[itype]=="structural") {
          if(!is.na(match(Dargs$error.model[itype],c("constant","exponential"))))
            resk<-sum((Dargs$yobs[tab_ytype==itype]-fk[tab_ytype==itype])**2) else {
              resk<-sum((Dargs$yobs[tab_ytype==itype]-fk[tab_ytype==itype])**2/cutoff(fk[tab_ytype==itype]**2,.Machine$double.eps))
            }
          statr[itype]<-statr[itype]+resk
        }
      }
    }
  }
  # Update sufficient statistics
  suffStat$statphi1<-suffStat$statphi1+opt$stepsize[kiter]*(stat1/Uargs$nchains-suffStat$statphi1)
  suffStat$statphi2<-suffStat$statphi2+opt$stepsize[kiter]*(stat2/Uargs$nchains-suffStat$statphi2)
  suffStat$statphi3<-suffStat$statphi3+opt$stepsize[kiter]*(stat3/Uargs$nchains-suffStat$statphi3)
  if(length(suffStat)>3) {
    for(i in 4:length(suffStat))
      suffStat[[i]] <- suffStat[[i]]+opt$stepsize[kiter]*(statr[i-3]/Uargs$nchains-suffStat[[i]])
  }
  
  ### delta_ik pour le calcul des SE 
  mu0 = betas[1,1]
  mu1 = betas[2,1]
  omega0 = varList$omega[1,1]
  omega1 = varList$omega[2,2]
  sigma = varList$pres[1]
  
  ni = sapply(unique(Dargs$IdM),function(i) length(which(Dargs$IdM==i))) 
  err = (Dargs$yobs-fpred)^2
  err2 = sapply(unique(Dargs$IdM), function(i) sum(err[Dargs$IdM==i]/sigma**3))
  
  deltamu0 = suffStat$statphi1[,1]/omega0-(mu0/omega0)
  deltamu1 = suffStat$statphi1[,2]/omega1-(mu1/omega1)
  delta_omega0 = -1/(2*omega0)+mu0**2/(2*omega0**2)+suffStat$statphi3[,1]/(2*omega0**2)-(mu0*suffStat$statphi1[,1])/omega0**2
  delta_omega1 = -1/(2*omega1)+mu1**2/(2*omega1**2)+suffStat$statphi3[,2]/(2*omega1**2)-(mu1*suffStat$statphi1[,2])/omega1**2
  delta_sigma = -ni/sigma + err2
  deltai_new = matrix(c(deltamu0,deltamu1,delta_omega0,delta_omega1,delta_sigma),ncol=length(unique(Dargs$IdM)),byrow = T)
  deltaik = (1-opt$stepsize[kiter])*deltai + opt$stepsize[kiter]*deltai_new
  
  ############# Maximisation
  ##### fixed effects
  
  if (opt$flag.fmin && kiter>=opt$nbiter.sa) {
    temp<-d1.omega[Uargs$ind.fix11,]*(t(Uargs$COV1)%*%(suffStat$statphi1-Uargs$dstatCOV[,varList$ind.eta]))
    betas[Uargs$ind.fix11]<-solve(comega[Uargs$ind.fix11,Uargs$ind.fix11],rowSums(temp))
    # ECO TODO: utiliser optimise dans le cas de la dimension 1
    #		if(length(Uargs$ind.fix10)>1)
    suppressWarnings(beta0<-optim(par=betas[Uargs$ind.fix10],fn=compute.Uy.multi,phiM=phiM,pres=varList$pres,args=Uargs,Dargs=Dargs,DYF=DYF,control=list(maxit=opt$maxim.maxiter))$par) # else
    #		beta0<-optimize(f=compute.Uy, interval=c(0.01,100)*betas[Uargs$ind.fix10],phiM=phiM,pres=varList$pres,args=Uargs,Dargs=Dargs,DYF=DYF)
    #		if(kiter==opt$nbiter.sa) {
    #		  cat("ind.fix10=",Uargs$ind.fix10,"ind.fix11=",Uargs$ind.fix11,"ind.fix1=",Uargs$ind.fix1,"ind.fix0=",Uargs$ind.fix0,"\n")
    #  		cat(betas,"\n")
    #		}
    betas[Uargs$ind.fix10]<-betas[Uargs$ind.fix10]+opt$stepsize[kiter]*(beta0-betas[Uargs$ind.fix10])
  } else {
    temp<-d1.omega[Uargs$ind.fix1,]*(t(Uargs$COV1)%*%(suffStat$statphi1-Uargs$dstatCOV[,varList$ind.eta]))
    betas[Uargs$ind.fix1]<-solve(comega[Uargs$ind.fix1,Uargs$ind.fix1],rowSums(temp))
  }
  
  varList$MCOV[Uargs$j.covariate]<-betas
  mean.phi<-Uargs$COV %*% varList$MCOV
  e1.phi<-mean.phi[,varList$ind.eta,drop=FALSE]
  
  # Covariance of the random effects
  omega.full<-matrix(data=0,nrow=Uargs$nb.parameters,ncol=Uargs$nb.parameters)
  omega.full[varList$ind.eta,varList$ind.eta]<-suffStat$statphi2/Dargs$N + t(e1.phi)%*%e1.phi/Dargs$N - t(suffStat$statphi1)%*%e1.phi/Dargs$N - t(e1.phi)%*%suffStat$statphi1/Dargs$N
  varList$omega[Uargs$indest.omega]<-omega.full[Uargs$indest.omega]
  
  # Simulated annealing (applied to the diagonal elements of omega)
  if (kiter<=opt$nbiter.sa) {
    diag.omega.full<-mydiag(omega.full)
    vec1<-diag.omega.full[Uargs$i1.omega2]
    vec2<-varList$diag.omega[Uargs$i1.omega2]*opt$alpha1.sa
    idx<-as.integer(vec1<vec2)
    varList$diag.omega[Uargs$i1.omega2]<-idx*vec2+(1-idx)*vec1
    varList$diag.omega[Uargs$i0.omega2]<-varList$diag.omega[Uargs$i0.omega2]*opt$alpha0.sa
  } else {
    varList$diag.omega<-mydiag(varList$omega)
  }
  varList$omega<-varList$omega-mydiag(mydiag(varList$omega))+mydiag(varList$diag.omega)
  
  # Residual error
  ytype = Dargs[["XM"]][["ytype"]][1:length(Dargs$yobs)]
  if(length(grep("structural",saemix.model["modeltype"]))>0) {
    i1<-0
    for(itype in 1:length(saemix.model["modeltype"])) {
      if(Dargs$modeltype[itype]=="structural") {
        i1<-i1+1
        if (Dargs$error.model[itype] %in% c("constant","exponential")) {
          sig2<-suffStat[[i1+3]]/length(ytype[ytype==itype])
          #sig2<-suffStat[[i1+3]]/Dargs$nobs
          if (kiter<=opt$nbiter.sa) {
            varList$pres[1+(i1-1)*2]<-max(varList$pres[1+(i1-1)*2]*opt$alpha1.sa,sqrt(sig2))
          } else {
            varList$pres[1+(i1-1)*2]<-sqrt(sig2)
          }
        }
        if (Dargs$error.model[itype]=="proportional") {
          sig2<-suffStat[[i1+3]]/length(ytype[ytype==itype])
          #sig2<-suffStat[[i1+3]]/Dargs$nobs
          if (kiter<=opt$nbiter.sa) {
            varList$pres[2+(i1-1)*2]<-max(varList$pres[2+(i1-1)*2]*opt$alpha1.sa,sqrt(sig2))
          } else {
            varList$pres[2+(i1-1)*2]<-sqrt(sig2)
          }
        }
        if(Dargs$error.model[itype]=="combined") { ## UNTESTED; need to compute the ssq for the response itype !!!
          #	if (Dargs$error.model=="combined") {
          # ECO TODO: check and secure (when fpred<0 => NaN)
          # JR: using lower=0 in the call to optim does not work, as L-BFGS-B
          # does not cope with non-finite function values that we are obviously
          # getting. Therefore we just take the absolute values after optimizing
          suppressWarnings(ABres<-abs(optim(par=varList$pres[c(1:2)+(i1-1)*2],fn=ssq,y=Dargs$yM[Dargs$XM$ytype==itype],f=fpred[Dargs$XM$ytype==itype],etype=Dargs$XM$ytype[Dargs$XM$ytype==itype])$par))
          #suppressWarnings(ABres<-abs(optim(par=varList$pres[c(1:2)+(i1-1)*2],fn=ssq,y=Dargs$yM,f=fpred,etype=Dargs$XM$ytype)$par))
          if (kiter<=opt$nbiter.sa) {
            for(i in 1:length(varList$pres)) varList$pres[i]<-max(varList$pres[i]*opt$alpha1.sa,ABres[i])
          }  else {
            if (kiter<=opt$nbiter.saemix[1]) {
              for(i in 1:length(varList$pres)) varList$pres[i]<-ABres[i]
            } else {
              for(i in 1:length(varList$pres)) varList$pres[i]<-varList$pres[i]+opt$stepsize[kiter]*(ABres[i]-varList$pres[i])
            }
          }
        }
      }
    }
  }
  # Modified to add SA to constant and exponential residual error models (Edouard Ollier 10/11/2016)
  
  return(list(varList=varList,mean.phi=mean.phi,phi=phi,betas=betas,suffStat=suffStat,deltai=deltaik))
}
