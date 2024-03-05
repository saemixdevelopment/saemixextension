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
    tab_ytype = Dargs$XM$ytype[1:length(Dargs$yobs)]   ### pour g?rer le nombre de cha?nes et que ytype soit de la meme longueur que yobs et fk 
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
  
  
  #### numerical gradient 
  coef = c(0,1)
  d = 0.000001
  
  if (kiter < opt$nbiter.saemix[1]){
    deltai_new = 0
  }
  else{
  
  nchains = Uargs$nchains
  mphiM = apply(phi,c(1,2),mean)  # mean of phi over all chains  
  ly = array(data = NA, dim=c(Dargs$N,length(coef),length(varList$pres[varList$pres!=0 & is.na(varList$pres)==F])))
  for(l in 1:length(coef)){
    for(j in 1:length(varList$pres[varList$pres!=0 & is.na(varList$pres)==F])){
      w = which(varList$pres!=0 & is.na(varList$pres)==F)[j]
      pres = varList$pres
      pres[w]<-varList$pres[w]+coef[l]*d
      ly[,l,j] = compute.LLy.multi2(mphiM,Uargs,Dargs,DYF,pres)
    }
  }
  
  delta_sigma  = -(ly[,2,]-ly[,1,])/d
  
  
  # mu and omega. Warnning: we suppose that omega is a diagonal matrix !! 
  
  phiM2 = matrix(mphiM[,Uargs$i1.omega2],ncol=length(Uargs$i1.omega2))
  
  lmu = array(data=NA,dim=c(Dargs$N,length(Uargs$ind.fix11),length(coef))) 
  mcov = t(t(Uargs$LCOV[Uargs$ind.fix11,Uargs$i1.omega2])%*%diag(c(betas[Uargs$ind.fix11]))) 
  ind = which(Uargs$LCOV[Uargs$ind.fix11,Uargs$i1.omega2]==1)
  compute.LLtheta = function(mu,omega,phiM2) -log(sqrt(omega))-(phiM2-(mu))**2 /(2*omega) # ici le mu représente mu+betaCOV
  omega = diag(omega.eta[1:length(Uargs$i1.omega2),1:length(Uargs$i1.omega2)])
  omegab = matrix(data=omega,nrow=Dargs$N,byrow = T,ncol = length(omega)) # mat des omega, même dimension que mu 
  for (l in 1:length(coef)){
    w = 1
    for (j in ind){
      indcol = ceiling(j/nrow(mcov))  # colonne dans phi, mu et omega correspondante
      mcov2 = mcov
      mcov2[j] = mcov[j]+ d*coef[l]
      mu = (Uargs$COV[,Uargs$ind.fix11]%*%mcov2)
      lmu[,w,l] = compute.LLtheta(mu[,indcol],omegab[,indcol],phiM2[,indcol])
      w = w+1
    }
  }
  delta_mu = (lmu[,,2]-lmu[,,1])/d
  # for omega now 
  
  mu = (Uargs$COV[,Uargs$ind.fix11]%*%mcov)
  delta_omega = (compute.LLtheta(mu,omegab+d,phiM2)-compute.LLtheta(mu,omegab,phiM2))/d
  
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
            indcol = which(LCOV2[sort(c(Uargs$ind.fix10))[j_], ]==1) # column in phi 
	    if (Dargs$transform.par[indcol]!=3){
              phiM3 = mphiM
              #phiM3[,Uargs$i0.omega2] = matrix(data=rep(betas[Uargs$indx.betaI[Uargs$i0.omega2],1],length(unique(Dargs$IdM[1:length(Dargs$yobs)]))),ncol = length(Uargs$i0.omega2),byrow = T)

              mcov2 = mcovTot
              mcov2[j] = mcov2[j]+ d*coef[l]
              fix = (Uargs$COV%*%mcov2) #
              phiM3[,indcol] = fix[,indcol] 
              lT[,l,w] = compute.LLy.multi2(phiM3,Uargs,Dargs,DYF,pres)
            }
            else{
              phiM3 = phiM
              #phiM3[,Uargs$ind.fix10] = matrix(data=rep(betas[Uargs$ind.fix10,1],length(unique(Dargs$IdM))),ncol = length(Uargs$ind.fix10),byrow = T)
              lT[,l,w] = compute.LLy.multi_selog(phiM3,Uargs,Dargs,DYF,pres,indcol,coef,l,d)
            }
          w = w+1
        }
      }
    
    delta_fix  = -(lT[,2,]-lT[,1,])/(d)
    deltai_new = matrix(c(delta_mu,delta_fix,delta_omega,delta_sigma),ncol=Dargs$N,byrow = T)
  }
  else{
    deltai_new = matrix(c(delta_mu,delta_omega,delta_sigma),ncol=Dargs$N,byrow = T)
  }
  }  
  deltaik = (1-opt$stepsize[kiter])*deltai + opt$stepsize[kiter]*deltai_new
  #print(deltaik)
  
  ############# Maximisation
  ##### fixed effects
  
  if (opt$flag.fmin && kiter>=opt$nbiter.sa) {
    temp<-d1.omega[Uargs$ind.fix11,]*(t(Uargs$COV1)%*%(suffStat$statphi1-Uargs$dstatCOV[,varList$ind.eta]))
    betas[Uargs$ind.fix11]<-solve(comega[Uargs$ind.fix11,Uargs$ind.fix11],rowSums(temp))
    # ECO TODO: utiliser optimise dans le cas de la dimension 1
    #		if(length(Uargs$ind.fix10)>1)
    ### coxpenalized 
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


### delta_ik pour le calcul des SE 
#mu0 = betas[1,1]
#mu1 = betas[2,1]
#h0 = exp(betas[3,1])
#alpha = betas[4,1]
#omega0 = varList$omega[1,1]
#omega1 = varList$omega[2,2]
#sigma = varList$pres[1]

#ni = sapply(unique(Dargs$IdM[1:length(Dargs$yobs)]),function(i) length(which(Dargs$IdM[1:length(Dargs$yobs)]==i))) - length(Dargs$modeltype[Dargs$modeltype=="likelihood"])# -1 pour donn?e survie 
#mpsiM<-transphi(mphiM,Dargs$transform.par)
#fpred<-structural.model(psiM, Dargs$IdM[1:length(Dargs$yobs)], Dargs$XM[1:length(Dargs$yobs),])
#err = (Dargs$yobs-fpred)^2
#err2 = sapply(unique(Dargs$IdM[1:length(Dargs$yobs)]), function(i) sum(err[Dargs$IdM[1:length(Dargs$yobs)]==i & Dargs$XM$ytype[1:length(Dargs$yobs)]==1]/sigma**3))

#di = Dargs$yobs[Dargs$XM$ytype==2]
#Ti = Dargs$XM$time[Dargs$XM$ytype==2]

#f_struc = function(t) suffStat$statphi1[,1]+suffStat$statphi1[,2]*t

#deltamu0 = suffStat$statphi1[,1]/omega0-(mu0/omega0)
#deltamu1 = suffStat$statphi1[,2]/omega1-(mu1/omega1)
#delta_omega0 = -1/(2*omega0)+mu0**2/(2*omega0**2)+suffStat$statphi3[,1]/(2*omega0**2)-(mu0*suffStat$statphi1[,1])/omega0**2
#delta_omega1 = -1/(2*omega1)+mu1**2/(2*omega1**2)+suffStat$statphi3[,2]/(2*omega1**2)-(mu1*suffStat$statphi1[,2])/omega1**2
#delta_sigma = -ni/sigma + err2
#delta_h0 = di/h0 - (exp(alpha*f_struc(Ti))-exp(alpha*suffStat$statphi1[,1]))/(alpha*suffStat$statphi1[,2])
#delta_alpha = di*f_struc(Ti)-(h0/(alpha**2*suffStat$statphi1[,2]))*(exp(alpha*f_struc(Ti))*(alpha*f_struc(Ti)-1)+exp(alpha*suffStat$statphi1[,1])*(1-alpha*suffStat$statphi1[,1]))
#deltai_new = matrix(c(deltamu0,deltamu1,delta_h0,delta_alpha,delta_omega0,delta_omega1,delta_sigma),ncol=length(unique(Dargs$IdM)),byrow = T)
#deltaik = (1-opt$stepsize[kiter])*deltai + opt$stepsize[kiter]*deltai_new

