################## Stochastic approximation - compute sufficient statistics (M-step) #####################
mstep.fim_stoch <-function(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat, deltai, covariance.model, task_fim.sa) {
  # M-step - stochastic approximation
  # M-step - stochastic approximation
  # Input: kiter, Uargs, structural.model, DYF, phiM, covariance.model
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
  
  
  
  if(task_fim.sa){
    #### numerical gradient 
    coef = c(0,1)
    d = 0.000001
    
    nchains = Uargs$nchains
    mphiM = apply(phi,c(1,2),mean)  # mean of phi over all chains  
    
    if(length(Uargs$ind.res)>0){
      # Residual error parameters: numerical gradient
     
      # ly = array(data = NA, dim=c(length(unique(Dargs$IdM)),length(coef),length(varList$pres[varList$pres!=0 & is.na(varList$pres)==F])))
      ly = array(data = NA, dim=c(Dargs$N,length(coef),length(varList$pres[varList$pres!=0 & is.na(varList$pres)==F])))
      for(l in 1:length(coef)){
        for(j in 1:length(varList$pres[varList$pres!=0 & is.na(varList$pres)==F])){
          w = which(varList$pres!=0 & is.na(varList$pres)==F)[j]
          pres = varList$pres
          pres[w]<-varList$pres[w]+coef[l]*d
          ly[,l,j] = compute.LLy.2(mphiM,Uargs,Dargs,DYF,pres)
        }
      }
      
      delta_sigma  = -(ly[,2,]-ly[,1,])/d
    }else{
      delta_sigma = NULL
    }

    
    # Fixed effects with iiv and Omega: formula
    phiM2 = matrix(mphiM[,Uargs$i1.omega2],ncol=length(Uargs$i1.omega2))
    
    mu_iiv = betas[Uargs$ind.fix11]
    n_mu_iiv = length(Uargs$ind.fix11)
    
    if(opt$flag.fmin && kiter>=opt$nbiter.sa){ #  omega.eta already contains only parameters with IIV
      Omega = omega.eta
    }else{
      Omega = omega.eta[Uargs$i1.omega2, Uargs$i1.omega2]
    }
    ncol_Omega = ncol(Omega)
    if(is.null(ncol_Omega)){
      if(!is.null(Omega)){
        ncol_Omega = 1
      }else{
        ncol_Omega = 0
      }
    }
    n_omega = 0.5*ncol_Omega*(ncol_Omega+1)
    Omega_inv = solve(Omega)
    Omega_Kronecker = kronecker(Omega_inv, Omega_inv)
    
    if(ncol_Omega <= 1){
      elim = 1
    }else{
      # Elimination matrix
      elim = elimination.matrix(n = ncol(Omega))
      elim.t = t(elim)
    }
    
    
    # Omega_estim.idx_mat  : Index of component in Omega matrix that are estimated
    Omega_estim.idx_mat = (TRUE==covariance.model[Uargs$i1.omega2, Uargs$i1.omega2])&(lower.tri(covariance.model[Uargs$i1.omega2, Uargs$i1.omega2], diag = TRUE))
    # Omega_estim.idx_vec  : Index of component in the vector of lowerTri Omega matrix that are estimated
    Omega_estim.idx_vec = Omega_estim.idx_mat[lower.tri(Omega_estim.idx_mat,diag=TRUE)]
    n_estim_omega = sum(Omega_estim.idx_vec) # number of omega components to be estimated
    
    res_deriv_mu_omega = do.call(rbind, lapply(X = 1:Dargs$N, FUN = deriv_mu_omega_i, 
                                               COV = Uargs$COV[,Uargs$ind.fix11],
                                               LCOV = Uargs$LCOV[Uargs$ind.fix11,Uargs$i1.omega2],
                                               phiM2, mu_iiv, elim, covariance.model, Omega_inv, Omega_Kronecker, Omega,
                                               n_mu_iiv
    )
    )
    
    # N x n_mu_iiv
    delta_mu = matrix(unlist(res_deriv_mu_omega[,1]), ncol = n_mu_iiv, byrow = TRUE)
    
    delta_omega = matrix(unlist(res_deriv_mu_omega[,2]), ncol = n_omega, byrow = TRUE)
    # keep only columns corresponding to omegas that are estimated -> N x n_estim_omega
    delta_omega = delta_omega[,Omega_estim.idx_vec]
    
    
    # Fixed effects without iiv: numerical gradient
    mcovTot = t(t(Uargs$LCOV)%*%diag(c(betas)))
    
    LCOV2 = Uargs$LCOV
    LCOV2[Uargs$ind.fix0,] = 0
    LCOV2[Uargs$ind.fix11, Uargs$i1.omega2] = 0 #set to 0 param with iiv
    ind = which(LCOV2==1) # index in LCOV2 of param without iiv

    
    if (length(Uargs$ind.fix10)>0){
      lT = array(data = NA, dim=c(length(unique(Dargs$IdM[1:length(Dargs$yobs)])),length(coef),length(c(Uargs$ind.fix10))))
      for(l in 1:length(coef)){
        w = 1
        for(j_ in 1:length(ind)){
          j = ind[j_] # index in mcovTot
          
          indcol = which(LCOV2[sort(c(Uargs$ind.fix10,Uargs$ind.fix0))[j_], ]==1) # column in phi 
          phiM3 = mphiM
          phiM3[,Uargs$i0.omega2] = matrix(data=rep(betas[Uargs$indx.betaI[Uargs$i0.omega2],1],length(unique(Dargs$IdM[1:length(Dargs$yobs)]))),ncol = length(Uargs$i0.omega2),byrow = T)

          mcov2 = mcovTot
          mcov2[j] = mcov2[j]+ d*coef[l]
          fix = (Uargs$COV%*%mcov2) #
          phiM3[,indcol] = fix[,indcol] 
          lT[,l,w] = compute.LLy.2(phiM3,Uargs,Dargs,DYF,pres)
          w = w+1
        }
      }
      
      delta_fix  = -(lT[,2,]-lT[,1,])/(d)
      deltai_new = matrix(c(delta_mu,delta_fix,delta_omega,delta_sigma),ncol=Dargs$N,byrow = T)
      
      # order
      deltai_new[Uargs$ind.fix11, ] = delta_mu
      deltai_new[sort(c(Uargs$ind.fix10,Uargs$ind.fix0)), ] = delta_fix
      
    }else{
      deltai_new = matrix(c(delta_mu,delta_omega,delta_sigma),ncol=Dargs$N,byrow = T)
    }
    
    
    deltaik = (1-opt$stepsize[kiter])*deltai + opt$stepsize[kiter]*deltai_new
  }else{
    deltaik = 0
  }
  
  
  ############# Maximisation
  ##### fixed effects
  
  if (opt$flag.fmin && kiter>=opt$nbiter.sa) {
    temp<-d1.omega[Uargs$ind.fix11,]*(t(Uargs$COV1)%*%(suffStat$statphi1-Uargs$dstatCOV[,varList$ind.eta]))
    betas[Uargs$ind.fix11]<-solve(comega[Uargs$ind.fix11,Uargs$ind.fix11],rowSums(temp))
    # ECO TODO: utiliser optimise dans le cas de la dimension 1
    #		if(length(Uargs$ind.fix10)>1)
    suppressWarnings(beta0<-optim(par=betas[Uargs$ind.fix10],fn=compute.Uy,phiM=phiM,pres=varList$pres,args=Uargs,Dargs=Dargs,DYF=DYF,control=list(maxit=opt$maxim.maxiter))$par) # else
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
  # Modified to add SA to constant and exponential residual error models (Edouard Ollier 10/11/2016)
  if(Dargs$modeltype=="structural") {
    if(length(Uargs$ind.res)==1) { # necessarily only one error model
      if (Dargs$error.model[1] %in% c("constant","exponential")) {
        sig2<-suffStat$statrese/Dargs$nobs
        if (kiter<=opt$nbiter.sa) {
          varList$pres[1]<-max(varList$pres[1]*opt$alpha1.sa,sqrt(sig2))
        } else {
          varList$pres[1]<-sqrt(sig2)
        }
      }
      if (Dargs$error.model[1]=="proportional") {
        sig2<-suffStat$statrese/Dargs$nobs
        if (kiter<=opt$nbiter.sa) {
          varList$pres[2]<-max(varList$pres[2]*opt$alpha1.sa,sqrt(sig2))
        } else {
          varList$pres[2]<-sqrt(sig2)
        }
      }
      
    } else {
      #	if (Dargs$error.model=="combined") {
      # ECO TODO: check and secure (when fpred<0 => NaN)
      # JR: using lower=0 in the call to optim does not work, as L-BFGS-B
      # does not cope with non-finite function values that we are obviously
      # getting. Therefore we just take the absolute values after optimizing
      suppressWarnings(ABres<-abs(optim(par=varList$pres,fn=ssq,y=Dargs$yM,f=fpred,etype=Dargs$XM$ytype)$par))
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
  
  return(list(varList=varList,mean.phi=mean.phi,phi=phi,betas=betas,suffStat=suffStat,deltai=deltaik))
}
