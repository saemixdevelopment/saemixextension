  for (kiter in 1:5) {
    xmcmc<-estep.multi(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
    varList<-xmcmc$varList
    DYF<-xmcmc$DYF
    phiM<-xmcmc$phiM
    # no M-step as stepsize==0
    allpar[(kiter+1),]<-allpar[kiter,]
    
    # Actually useless, theta itself never used ?? => used here to print results
    if(length(grep("structural",Dargs$modeltype))>0)
      theta<-c(fixed.psi,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res]) else
        theta<-c(fixed.psi,var.eta[Uargs$i1.omega2])
    #print(theta)
  }
  
