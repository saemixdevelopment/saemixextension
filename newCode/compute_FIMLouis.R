# WIP... 
# Compute LLi (complete likelihood)


# LLy: compute.LLy
compute.LLy<-function(phiM,args,Dargs,DYF,pres) {
  # input
  ## phiM: individual parameters
  ## Dargs: structural model, data (XM and IdM)
  ## args: which elements of DYF are to be filled
  ## pres: residual error model (when y Gaussian)
  ## DYF: matrix storing the likelihood for the individual observations
  # returns
  ## U: vector with the contribution of each individual to the likelihood
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

# LLphi: portion of the complete log-likelihood for the parameter
compute.LLphi<-function(phiM, theta, omega) {
  # input
  ## phiM: individual parameters (a matrix with nsuj row and npar parameters)
  ## theta: fixed effects (mu and beta)
  ## omega: variance-covariance matrix
  # returns
  ## llphi: vector with the contribution of each subject's individual parameters to LLi
  
  
  
  return(llphi)
}

gradient.LLphi<-function(phiM,args,Dargs) {
  # returns a vector with the contribution of each subject's individual parameter to d_theta(LLi)
  
  return(dLLphi)
}
hessian.LLphi<-function(phiM,args,Dargs) {
  # returns a vector with the contribution of each subject's individual parameter to d2_theta(LLi)
  
  return(d2LLphi)
}
