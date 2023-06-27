mstep<-function(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat) {
  
}

mstep<-function(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat) {
  # M-step - stochastic approximation
  # Input: kiter, Uargs, structural.model, DYF, phiM (unchanged)
  # Output: varList, phi, betas, suffStat (changed)
  #					mean.phi (created)
  
  # Update variances - TODO - check if here or elsewhere
  nb.etas<-length(varList$ind.eta)
  # replace diagonal elements lower than .Machine$double.eps with .Machine$double.eps => omega.eta
  domega<-cutoff(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]),.Machine$double.eps)
  omega.eta<-varList$omega[varList$ind.eta,varList$ind.eta,drop=FALSE]
  omega.eta<-omega.eta-mydiag(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]))+mydiag(domega)
  #  print(varList$omega.eta)
  # chol.omega<-try(chol(omega.eta)) # not used...
  d1.omega<-Uargs$LCOV[,varList$ind.eta]%*%solve(omega.eta)
  d2.omega<-d1.omega%*%t(Uargs$LCOV[,varList$ind.eta])
  comega<-Uargs$COV2*d2.omega
  
  # Parameters
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  #	if(Dargs$error.model=="exponential")
  #		fpred<-log(cutoff(fpred))
  ff<-matrix(fpred,nrow=Dargs$nobs,ncol=Uargs$nchains)
  for(k in 1:Uargs$nchains) phi[,,k]<-phiM[((k-1)*Dargs$N+1):(k*Dargs$N),]
  # overall speed similar
  #    phi<-aperm(array(phiM,c(N,nchains,3)),c(1,3,2))
  stat1<-apply(phi[,varList$ind.eta,,drop=FALSE],c(1,2),sum) # sum on columns ind.eta of phi, across 3rd dimension
  stat2<-matrix(data=0,nrow=nb.etas,ncol=nb.etas)
  stat3<-apply(phi**2,c(1,2),sum) #  sum on phi**2, across 3rd dimension
  statr<-0
  for(k in 1:Uargs$nchains) {
    phik<-phi[,varList$ind.eta,k]
    stat2<-stat2+t(phik)%*%phik
    fk<-ff[,k]
    if(length(Dargs$error.model)==1) {
      if(!is.na(match(Dargs$error.model,c("constant","exponential"))))
        resk<-sum((Dargs$yobs-fk)**2) else {
          if(Dargs$error.model=="proportional") {
            #		        idx.okpred<-which(fk>.Machine$double.eps)
            #		        vec<-(Dargs$yobs-fk)**2/cutoff(fk**2,.Machine$double.eps)
            #		        resk<-sum(vec[idx.okpred])
            resk<-sum((Dargs$yobs-fk)**2/cutoff(fk**2,.Machine$double.eps))
          } else resk<-0
        }
    } else resk<-0
    statr<-statr+resk
  }
  # Update sufficient statistics
  suffStat$statphi1<-suffStat$statphi1+opt$stepsize[kiter]*(stat1/Uargs$nchains-suffStat$statphi1)
  suffStat$statphi2<-suffStat$statphi2+opt$stepsize[kiter]*(stat2/Uargs$nchains-suffStat$statphi2)
  suffStat$statphi3<-suffStat$statphi3+opt$stepsize[kiter]*(stat3/Uargs$nchains-suffStat$statphi3)
  suffStat$statrese<-suffStat$statrese+opt$stepsize[kiter]*(statr/Uargs$nchains-suffStat$statrese)
  
  
  ############# Maximisation
  ##### fixed effects

  # Eco: changed COV1 to COV[,index.fixedpar.fix], both elements will be in SaemixIndividualModel  
  # in main.R, COV1 redefined as COV[,ind.fix11] after being COV[,ind.fix1]
  # redefinition occurs at 
  # if(opt$flag.fmin && kiter==saemix.options$nbiter.sa) => here we can adjust without needing COV1
  
  # Big problem: in original algorithm, dstatCOV is never updated, and that doesn't make sense (unless we always refer to the original value of dstatCOV ??)
  # originally defined as COV[,index.fixedpar.estim,drop=FALSE]%*%MCOV[index.fixedpar.estim,,drop=FALSE]
  # MCOV updated in M-step but not dstatCOV... [note: also the case in the original algorithm]
  # Proposal to compute dstatCOV anew in M-step, compare  results ?
  ## NO ! in fact dstatCOV is the submatrix corresponding to parameters not estimated, so they never change !
  
  # dstatCOV <- Uargs1$dstatCOV
  # dstatCOV <- individualmodel@dstatCOV
#  dstatCOV <- Uargs$COV[,Uargs$index.fixedpar.fix,drop=FALSE] %*% Uargs$MCOV[Uargs$index.fixedpar.fix,,drop=FALSE]
  if (opt$flag.fmin && kiter>=opt$nbiter.sa) {
#    temp<-d1.omega[Uargs$ind.fix11,]*(t(Uargs$COV[,Uargs$index.fixedpar.fix])%*%(suffStat$statphi1-Uargs$dstatCOV[,varList$ind.eta]))
    temp<-d1.omega[Uargs$ind.fix11,] * (t(Uargs$COV[,Uargs$index.fixedpar.fix]) %*% (suffStat$statphi1-dstatCOV[,varList$ind.eta]))
    beta0<-optim(par=betas[Uargs$ind.fix10],fn=compute.Uy,phiM=phiM,pres=varList$pres,args=Uargs,Dargs=Dargs,DYF=DYF,control=list(maxit=opt$maxim.maxiter))$par
    betas[Uargs$ind.fix10]<-betas[Uargs$ind.fix10]+opt$stepsize[kiter]*(beta0-betas[Uargs$ind.fix10])
    } else {
    temp<-d1.omega[Uargs$ind.fix1,] * (t(Uargs$COV[,Uargs$ind.fix1]) %*% (suffStat$statphi1-dstatCOV[,varList$ind.eta]))
    betas[Uargs$ind.fix1]<-solve(comega[Uargs$ind.fix1,Uargs$ind.fix1],rowSums(temp)) 
  }

  
}