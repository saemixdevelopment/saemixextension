
compute.LLy<-function(phiM,args,Dargs,DYF,pres) {
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  if (Dargs$modeltype=="structural"){
    gpred<-error(fpred,pres,Dargs$XM$ytype)
    DYF[args$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
  } else {
    DYF[args$ind.ioM]<- -fpred
  }
  U<-colSums(DYF)
  return(U)
}


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

# Generate design matrices at each varlevel
## dans data ou dans une fonction qui associe data et model ?

