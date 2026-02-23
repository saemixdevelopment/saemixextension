
########################### Computational functions to recode

if(FALSE) {
  
  compute.Uy<-function(b0,phiM,pres,Margs,Dargs,DYF) {
    # Attention, DYF variable locale non modifiee en dehors
    args$MCOV0[args$j0.covariate]<-b0
    phi0<-args$COV0 %*% args$MCOV0
    phiM[,args$i0.omega2]<-do.call(rbind,rep(list(phi0),args$nchains))
    psiM<-transphi(phiM,Dargs$transform.par)
    if (Margs$model.type[[i]]=="continuous"){
      fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
      for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
      gpred<-error(fpred,pres,Dargs$XM$ytype)
      DYF[args$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
    } else {
      fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
      for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
      DYF[args$ind.ioM]<- -fpred
    }
    U<-sum(DYF)
    return(U)
  }
  
  conditional.distribution_c<-function(phi1,phii,idi,xi,yi,mphi,idx,iomega,trpar,model,pres,err) {
    phii[idx]<-phi1
    psii<-transphi(matrix(phii,nrow=1),trpar)
    if(is.null(dim(psii))) psii<-matrix(psii,nrow=1)
    fi<-model(psii,idi,xi)
    #  if(err=="exponential") # Reverted this bit to the previous version to avoid a compiler error, not sure why it was changed...
    #    fi<-log(cutoff(fi))
    #  if (!(is.null(pres)) && pres[1] == pres) { # package compile throws an error when comparing a vector of length 2 (pres) to a vector of length 1
    # gi <- cutoff(pres[1])
    # } 
    # else{
    #   gi<-error(fi,pres) #    cutoff((pres[1]+pres[2]*abs(fi)))
    # }
    ind.exp<-which(err=="exponential")
    for(ityp in ind.exp) 
      fi[xi$ytype==ityp]<-log(cutoff(fi[xi$ytype==ityp]))
    gi<-error(fi,pres,xi$ytype)      #    cutoff((pres[1]+pres[2]*abs(fi)))
    Uy<-sum(0.5*((yi-fi)/gi)**2+log(gi))
    dphi<-phi1-mphi
    Uphi<-0.5*sum(dphi*(dphi%*%iomega))
    return(Uy+Uphi)
  }
  
  conditional.distribution_d<-function(phi1,phii,idi,xi,yi,mphi,idx,iomega,trpar,model) {
    phii[idx]<-phi1
    psii<-transphi(matrix(phii,nrow=1),trpar)
    if(is.null(dim(psii))) psii<-matrix(psii,nrow=1)
    fi<-model(psii,idi,xi)
    Uy <- sum(-fi)
    dphi<- phi1-mphi
    Uphi<- 0.5*sum(dphi*(dphi%*%iomega))
    return(Uy+Uphi)
  }
  
  
}