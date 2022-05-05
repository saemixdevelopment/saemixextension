saemixDir<-"/home/eco/work/saemix/saemixextension"
source(file.path(saemixDir,"R","aaa_generics.R"))

########################################################### Data
# Creating an SaemixData object with PK/PD data by hand - Responses
source(file.path(saemixDir,"Rext","SaemixData.R"))
source(file.path(saemixDir,"Rext","SaemixOutcome.R"))
pkpd<-read.table(file.path(saemixDir,"data40","warfarinPKPD.tab"), header=T)

x<-new(Class="SaemixData", name.data="pkpd")
x1<-continuousOutcome(model="combined2", start=c(1, 0.2))
x2<-continuousOutcome(start=c(2))
#new(Class="SaemixOutcome",type="continuous")
x@outcome<-list(conc=new(Class="SaemixOutcome",type="continuous"), effect=new(Class="SaemixOutcome",type="continuous"))

zesuj<-unique(pkpd$id)
ind1<-match(pkpd$id, zesuj)
x@data<-data.frame(index=ind1, pkpd[,c(1:4)], ytype=pkpd$dvid, cens=0, mdv=0, pkpd[,c(6:8)])
x@name.group<-"id"
x@name.predictors<-c("time","amt")
x@name.X<-"time"
x@name.response<-c("dv")
x@name.ytype<-"ytype"
x@name.covariates<-colnames(pkpd)[6:8]
x@name.cens<-"cens"
x@name.mdv<-"mdv"
x@N<-length(zesuj)
nind.obs<-tapply(pkpd$id, pkpd$id, length)
x@nind.obs<-c(nind.obs[match(zesuj,names(nind.obs))])
x@ntot.obs<-sum(x@nind.obs)
x@units<-c("hr","mg")

pkpd.saemix<-x

########################################################### Model
# Model outcomes
xout1<-new(Class="SaemixContinuousOutcome", error.model=x1$error.model, error.npar=x1$error.npar, error.function=x1$error.function, error.parameters=x1$start, error.fix=x1$error.fix)
xout2<-new(Class="SaemixContinuousOutcome", error.model=x2$error.model, error.npar=x2$error.npar, error.function=x2$error.function, error.parameters=x2$start, error.fix=x2$error.fix)
pkpd.outcome<-list(conc=xout1, effect=xout2)

# Model function
model1cptdirect<-function(psi,id,xidep) { 
  tim<-xidep[,1]
  dose<-xidep[,2]
  ytype<-xidep$ytype
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  ic50<-psi[id, 4]
  k<-CL/V
  ypk<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  ypd<-100*(1-ypk/(ypk+ic50))
  ypk[ytype==2]<-ypd[ytype==2]
  return(ypk)
}

## Test the function
xidep<-pkpd.saemix@data[,c(pkpd.saemix@name.predictors, "ytype")]
id1<-pkpd.saemix@data$index
psi0<-c(ka=1, vd=5, cl=0.1, ic50=5)
psi1<-do.call(rbind, rep(list(psi0), pkpd.saemix@N))
ypred<-model1cptdirect(psi1, id1, xidep)

# Test one error model on all the predictions
f1<-runif(10, 0, 10)
g1<-xout1@error.function(f1, xout1@error.parameters)
if(max(abs(g1-sqrt((f1*x1$start[2])**2+x1$start[1]**2)))>10^(-10)) cat("Problem in computing g\n")

# Beginning of the model class - SaemixStructuralModel
source(file.path(saemixDir,"Rext","SaemixVarLevel.R"))
source(file.path(saemixDir,"Rext","SaemixModel.R"))
source(file.path(saemixDir,"Rext","SaemixModel-methods.R"))
# pkpd.model<-new(Class="SaemixStructuralModel", model=model1cptdirect, outcome=list(conc=continuousOutcome(model="combined2", start=c(1, 0.2)), effect=continuousOutcome(start=c(2))))

# Model parameters for the following (population parameters for every subject in the dataset)
pkpd.psi0<-c(ka=1, vd=5, cl=0.1, ic50=5)
psi1<-do.call(rbind, rep(list(pkpd.psi0), pkpd.saemix@N))

# Variability model
mat1.start<-mat1<-diag(4)
mat1[2,3]<-mat1[3,2]<-1
mat1.fix<-diag(c(0,0,0,1))
mat1.start[mat1==1]<-c(0.3, 0.5, 0.35, 0.35, 0.5, 1)
var1<-saemixVarModel(omega=mat1.start,omega.model=mat1, omega.model.fix=mat1.fix)

########################################################### SaemixModel

pkpd.model<-saemixModel(model=model1cptdirect, outcome=pkpd.outcome, parameter=pkpd.psi0, var.model=var1)

########################################################### Predictions and error
# ssq, compute.LLy
source(file.path(saemixDir,"Rext","func_aux.R"))
source(file.path(saemixDir,"Rext","saemixControl.R"))

# debuging ssq, compute.LLy, optim
if(FALSE)
  source(file.path(saemixDir,"newCode","debug_ssqAcceptanceRatio.R"))

#### Eco 12/04/2022 - stopped here

########################################################### Used in main_mstep.R to optimise beta0 (???)
# Fixed parameters ? (ie covariate effects and parameters without IIV ?)

# Input
## b0: parameters we want to optimise
## args: list, use j0.covariate, COV0, MCOV0, i0.omega2, nchains, ind.ioM
## phiM: matrix of parameters phi
## pres: vector of residual error parameters
## Dargs: list, use transform.par, modeltype, structural.model(), etype.exp, XM, yM
## DYF: a vector, modified locally but not exported

compute.Uy.old<-function(b0,phiM,pres,args,Dargs,DYF) {
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

# Input
## phiM: matrix of parameters
## Dargs: list, here we use the elements transform.par, model (structural model), IdM, XM, yM, ind.ioM (a vector of 0/1 indicating which elements of the square matrix DYF have data), list of outcomes (type, error.model, error.function)
## error.parameters: list of the error parameters for the different outcomes in the model (NULL for non-continuous outcome)
## DYF: a matrix of dimension nmax=max(n_i) times (N*nb.chains); the column for subject i has the last nmax-n_i elements set to 0, so that colSums sums on the observations for that subject 
# Output 
## U: vector with either sum(-LL) (continuous models, minus the constant sum(log(1/sqrt(2*pi)))) or sum(-logpdf) (likelihood models)

compute.Uy.old<-function(b0,phiM,pres,args,Dargs,DYF) {
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


compute.LLy<-function(phiM, Dargs, DYF, error.parameters) {
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$model(psiM,Dargs$IdM,Dargs$XM)
  lpred<-fpred
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  for(iout in 1:length(Dargs$outcome)) {
    idx1<-which(Dargs$XM$ytype==iout)
    # print(summary(fpred[idx1]))
    if(Dargs$outcome[[iout]]@type=="continuous") {
      if(Dargs$outcome[[iout]]@error.model=="exponential") fpred[idx1]<-log(cutoff(fpred[idx1]))
      gpred<-Dargs$outcome[[iout]]@error.function(fpred[idx1], error.parameters[[iout]])
      lpred[idx1]<-0.5*((Dargs$yM[idx1]-fpred[idx1])/gpred)**2+log(gpred)
      # print(summary(gpred))
    } else lpred[idx1]<- (-fpred[idx1])
  }
  DYF[Dargs$ind.ioM]<-lpred
  U<-colSums(DYF)
  return(U)
}


########################################################### Used in main_estep to optimise phi1.opti

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


########################################################### One longitudinal + one event

x1<-continuousOutcome(model="proportional", start=c(0.3))
x2<-discreteOutcome(type="event")
x3<-discreteOutcome()
xout<-new(Class="SaemixContinuousOutcome", error.model=x1$error.model, error.npar=x1$error.npar, error.function=x1$error.function, error.parameters=x1$start, error.fix=x1$error.fix)
x2<-new(Class="SaemixDiscreteOutcome", type="event")
x3<-new(Class="SaemixDiscreteOutcome", type="binary")



# Truncated gaussian distribution (verifie par rapport a definition de erf/matlab)
normcdf<-function(x,mu=0,sigma=1)
  cutoff(pnorm(-(x-mu)/sigma,lower.tail=FALSE),1e-30)


