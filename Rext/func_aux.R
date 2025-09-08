########################### Computational functions
# error<-function(f,ab,etype) { # etype: error model
#   g<-f
#   for(ityp in sort(unique(etype))) {
#     g[etype==ityp]<-error.typ(f[etype==ityp],ab[((ityp-1)*2+1):(ityp*2)])
#   }
#   return(g)
# }

# lerror.model: a list of SaemixErrorModel (for continuous outcome) and "none" (for non-continuous outcome)
error<-function(f,lerror.model,etype) { # etype: identifiant outcome
  g<-f
  for(ityp in sort(unique(etype))) {
    if(is(lerror.model[[ityp]],"SaemixErrorModel"))
      g[etype==ityp]<- lerror.model[[ityp]]@model(f[etype==ityp], lerror.model[[ityp]]@param)
  }
  return(g)
}

# error.model: error model for one outcome
ssq<-function(ab,y,f,error.model) { # Sum of squares; need to put ab first as these parameters are optimised by optim
  g<-error.model@model(f,ab)
  e<-sum(((y-f)**2/g**2)+2*log(g))
  return(e)
}


# not needed anymore
# error.typ<-function(f,ab) {
#   #  g<-cutoff(ab[1]+ab[2]*abs(f))
#   g<-cutoff(sqrt(ab[1]^2+ab[2]^2*f^2))  # Johannes 02/21
#   return(g)
# }

# error.model now a list of SaemixErrorModel (for continuous outcome) and "none" (for non-continuous outcome)
## pres now contained in Dargs$error.model
compute.LLy <- function(phiM, args, Dargs, DYF) {
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  if (Dargs$modeltype=="structural"){
    gpred<-error(fpred,Dargs$error.model,Dargs$XM$ytype)
    DYF[args$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
  } else {
    DYF[args$ind.ioM]<- -fpred
  }
  U<-colSums(DYF)
  return(U)
}

########################### Computational functions to recode

if(FALSE) {
  
  compute.Uy<-function(b0,phiM,pres,args,Dargs,DYF) {
    # Attention, DYF variable locale non modifiee en dehors
    args$MCOV0[args$j0.covariate]<-b0
    phi0<-args$COV0 %*% args$MCOV0
    phiM[,args$i0.omega2]<-do.call(rbind,rep(list(phi0),args$nchains))
    psiM<-transphi(phiM,Dargs$transform.par)
    if (Dargs$modeltype=="structural"){
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