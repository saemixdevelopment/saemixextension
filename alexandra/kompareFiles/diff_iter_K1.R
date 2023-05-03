 ################################################# SA iterations after burn-in (kiter=6 to 150)
  for (kiter in 6:149) {
    xmcmc<-estep.multi(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
    varList<-xmcmc$varList
    DYF<-xmcmc$DYF
    phiM<-xmcmc$phiM
    
    if(opt$stepsize[kiter]>0) {
      ############# Stochastic Approximation
      xstoch<-mstep.multi(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat, deltai)
      varList<-xstoch$varList
      mean.phi<-xstoch$mean.phi
      phi<-xstoch$phi
      betas<-xstoch$betas
      suffStat<-xstoch$suffStat
      deltai = xstoch$deltai
      
      beta.I<-betas[Uargs$indx.betaI]
      fixed.psi<-transphi(matrix(beta.I,nrow=1),saemix.model["transform.par"])
      betaC<-betas[Uargs$indx.betaC]
      var.eta<-mydiag(varList$omega)
      l1<-betas.ini
      l1[Uargs$indx.betaI]<-fixed.psi
      l1[Uargs$indx.betaC]<-betaC
      
      if(length(grep("structural",Dargs$modeltype))>0)
        allpar[(kiter+1),]<-c(l1,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res]) else
          allpar[(kiter+1),]<-c(l1,var.eta[Uargs$i1.omega2])
      
    } else  #end of loop on if(opt$stepsize[kiter]>0)
      allpar[(kiter+1),]<-allpar[kiter,]
    print(allpar[(kiter+1),])
  }
