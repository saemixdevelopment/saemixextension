    # Burn-in - resetting sufficient statistics
    if(opt$flag.fmin && kiter==saemix.options$nbiter.sa) {
      Uargs$COV1<-Uargs$COV[,Uargs$ind.fix11]
      ind.prov<-!(varList$ind.eta %in% Uargs$i0.omega2)
      varList$domega2<-varList$domega2[ind.prov,ind.prov,drop=FALSE] # keep in domega2 only indices of parameters with IIV
      varList$ind0.eta<-Uargs$i0.omega2
      varList$ind.eta<-1:(Uargs$nb.parameters)  	
      if(length(varList$ind0.eta)>0) varList$ind.eta<-varList$ind.eta[!(varList$ind.eta %in% varList$ind0.eta)] # update ind.eta, now only parameters with IIV
      Uargs$nb.etas<-length(varList$ind.eta)
      suffStat$statphi1<-0
      suffStat$statphi2<-0
      suffStat$statphi3<-0
    }
    
    # E-step
    xmcmc<-estep.multi(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
    varList<-xmcmc$varList
    DYF<-xmcmc$DYF
    phiM<-xmcmc$phiM
    #  psiM<-transphi(phiM,saemix.model["transform.par"])
    
    # M-step
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
    } else { #end of loop on if(opt$stepsize[kiter]>0)
      allpar[(kiter+1),]<-allpar[kiter,]
    }
    # End of loop on kiter
  }
  
  ## inserer calcul Delattre & Kuhn 
  d = array(data=NA,dim=c(Uargs$nb.parest[1],Uargs$nb.parest[1],Dargs$N))
  for (i in 1:Dargs$N){
    ddi = deltai[,i] %*% t(deltai[,i])
    d[,,i] = ddi
  }
  
  fim = 0
  for (i in 1:Dargs$N){
    fim = fim+d[,,i]
  }
  inv_fim=try(solve(fim))
  
