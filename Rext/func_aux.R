########################### Computational functions
# lerror.model: a list of SaemixErrorModel (for continuous outcome) and "discrete" (for non-continuous outcome)
# error.model: error model for one outcome
error<-function(f,lerror.model,etype) { # etype: identifiant outcome
  g<-f
  for(ityp in sort(unique(etype))) {
    if(is(lerror.model[[ityp]],"SaemixErrorModel"))
      g[etype==ityp]<- lerror.model[[ityp]]@model(f[etype==ityp], lerror.model[[ityp]]@param)
  }
  return(g)
}
# ToDo: debug error, adjust to multiple outcomes

# not needed anymore
# error.typ<-function(f,ab) {
#   #  g<-cutoff(ab[1]+ab[2]*abs(f))
#   g<-cutoff(sqrt(ab[1]^2+ab[2]^2*f^2))  # Johannes 02/21
#   return(g)
# }

# Sum of squares; need to put ab first as these parameters are optimised by optim
## pres now contained in Dargs$error.model
# error.model: error model for one outcome, contains a model() function computing gpred given f and parameters ab
ssq<-function(ab,y,f,error.model) {
  g<-error.model@model(f,ab)
  e<-sum(((y-f)**2/g**2)+2*log(g))
  return(e)
}

# transformPar (defined in SaemixIndivModel) replaces:
## transphi (with param=phi and transform=model@transform)
## transpsi (with param=psi and transform=model@invtransform)
## dtransphi (with param=phi and transform=model@dtransform) [check format and maybe create different function]

############## Likelihood for observations y
# transphi replaced by transformPar

compute.LLy <- function(phiM, Margs, Dargs, DYF) {
  psiM<-transformPar(phiM,Margs$transform)
  fpred<-Margs$structural.model(psiM,Dargs$IdM,Dargs$XM)
  for(ityp in 1:length(Margs$model.type)) {
    if (Margs$model.type[[ityp]]=="continuous"){
      if(Margs$error.model[[ityp]]@name=="exponential") fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
      gpred<-error(fpred,Margs$error.model,Dargs$XM$ytype)
      DYF[Dargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
    } else {
      DYF[Dargs$ind.ioM]<- -fpred
    }  
  }
  U<-colSums(DYF)
  return(U)
}

## functions recoded separating lists by type and usage
## Margs: model related list
## Dargs: data related list

