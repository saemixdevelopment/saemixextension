########################################################### Folder and generics
saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
progDirExt<-file.path(saemixDir,"Rext")
source(file.path(progDir,"aaa_generics.R"))

########################################################### Data
# Creating an SaemixData object with PK/PD data by hand - Responses
source(file.path(progDirExt,"SaemixOutcome.R"))
source(file.path(progDirExt,"SaemixOutcome-methods.R"))
source(file.path(progDirExt,"SaemixData.R"))
source(file.path(progDirExt,"SaemixCovariateModel.R"))
source(file.path(progDirExt,"SaemixCovariate.R"))
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

# Mean value for population parameters
pkpd.psi0<-c(ka=1, vd=5, cl=0.1, ic50=5)
# Table with population parameters for every subject in the dataset)
psi1<-do.call(rbind, rep(list(pkpd.psi0), pkpd.saemix@N))

## Test the function
xidep<-pkpd.saemix@data[,c(pkpd.saemix@name.predictors, "ytype")]
id1<-pkpd.saemix@data$index
ypred<-model1cptdirect(psi1, id1, xidep)

# Model outcomes
pkpd.outcome<-list(conc=saemixOutcome(unit="mg/L", model="combined2", start=c(1, 0.2)), 
                   effect=saemixOutcome(unit="%", start=1))

gpred<-ypred
for(iout in 1:length(pkpd.outcome)) {
  if(pkpd.outcome[[iout]]@type=="continuous")
    gpred[pkpd.saemix@data$ytype==iout] <- pkpd.outcome[[iout]]@error.function(ypred[pkpd.saemix@data$ytype==iout], pkpd.outcome[[iout]]@error.parameters) else gpred[pkpd.saemix@data$ytype==iout]<-NA
}

########################################################### Optimising fixed effects
# Computing SD of predictions for structural models
source(file.path(saemixDir,"Rext","func_aux.R"))
source(file.path(saemixDir,"Rext","saemixControl.R"))

# Applying ssq to the two outcomes
iout<-1
ssq(pkpd.outcome[[iout]]@error.parameters, pkpd.saemix@data[pkpd.saemix@data$ytype==iout, pkpd.saemix@name.response], ypred[pkpd.saemix@data$ytype==iout],pkpd.outcome[[iout]]@error.function)

iout<-2
ssq(pkpd.outcome[[iout]]@error.parameters, pkpd.saemix@data[pkpd.saemix@data$ytype==iout, pkpd.saemix@name.response], ypred[pkpd.saemix@data$ytype==iout],pkpd.outcome[[iout]]@error.function)

# Optimising parameters
iout<-1
y1<-pkpd.saemix@data[pkpd.saemix@data$ytype==iout, pkpd.saemix@name.response]
fpred1<-ypred[pkpd.saemix@data$ytype==iout]
ab<-pkpd.outcome[[iout]]@error.parameters
ABres<-optim(par=ab,fn=ssq,y=y1,f=fpred1, error.function=pkpd.outcome[[iout]]@error.function)$par

## Checking
par(mfrow=c(1,2))
x<-seq(0.1,1, length.out=20)
ssq1<-c()
par1<-ABres
for(ix in x) {
  ssq1<-c(ssq1, ssq(c(ix,par1[2]), y1, fpred1, pkpd.outcome[[iout]]@error.function))
}
plot(x, ssq1, type="b", main="Varying a with b fixed to estimate")
points(ABres[1], ssq(c(ABres[1],par1[2]), y1, fpred1, pkpd.outcome[[iout]]@error.function), pch=20, col="red")

ssq2<-c()
par1<-ABres
for(ix in x) {
  ssq2<-c(ssq2, ssq(c(par1[1],ix), y1, fpred1, pkpd.outcome[[iout]]@error.function))
}

plot(x, ssq2, type="b", main="Varying b with a fixed to estimate")
points(ABres[2], ssq(c(par1[1],ABres[2]), y1, fpred1, pkpd.outcome[[iout]]@error.function), pch=20, col="red")

## Quick check 2 (we could probably plot a 3-D curve but not really worth it)
minSSQ<-ssq(ABres, y1, fpred1, pkpd.outcome[[iout]]@error.function)
a1<-rnorm(100, mean=ABres[1], sd=ABres[1]*.3)
b1<-rnorm(100, mean=ABres[2], sd=ABres[2]*.3)
ssq1<-c()
for(ia in a1) {
  for(ib in b1)
    ssq1<-c(ssq1, ssq(c(ia, ib), y1, fpred1, pkpd.outcome[[iout]]@error.function))
}
if(min(ssq1)<minSSQ) cat("Optimisation didn't work\n")

# Sufficient statistic for constant error model
ABres.add<-optim(par=c(1),fn=ssq,y=y1,f=fpred1, error.function=constantErrorModel, lower=0, upper=10, method="Brent")$par

ssq1<-ssq(1, y1, fpred1, constantErrorModel)

cat("Optimising with optim and method='Brent', sigma=", ABres.add,"\n")
cat("Computing the sufficient statistic, sigma=",sqrt(ssq1/length(y1)),"\n")

# Sufficient statistic for proportional error model - actually, no sufficient statistic because of the log(1/sqrt(2*pi)/(sigma*f) term...)
ABres.prop<-optim(par=c(0.3),fn=ssq,y=y1,f=fpred1, error.function=proportionalErrorModel, lower=0, upper=10, method="Brent")$par

ssq1<-sumSquare(1, y1, fpred1, proportionalErrorModel)
sqrt(ssq1/length(y1))

cat("Optimising with optim and method='Brent', sigma=", ABres.prop,"\n")
cat("Computing the sufficient statistic, sigma=",sqrt(ssq1/length(y1)),"\n")

# Optimisation with only 1 residual error parameter
user.error2<-function(f,ab) {
  g<-cutoff(sqrt((0.5+ab[1]*f)^2))
  return(g)
}

pstart<-0.3
ABres.one<-optim(par=pstart,fn=ssq,y=y1,f=fpred1, error.function=user.error2, lower=0, upper=pstart*10, method="Brent")$par
if(ABres.one>pstart*9.99) cat("Change upper boundary\n")

optim(par=pstart,fn=ssq,y=y1,f=fpred1, error.function=user.error2, lower=0, upper=10^8, method="Brent")$par

####
# So, need to replace the optimisation step
## for each continuous response
## use sufficient statistics for constant and proportional erroe models
## use optim with Brent method if npar=1, setting boundaries
## use regular optim otherwise

########################################################### Compute likelihood, used in main_estep.R and in func_distcond.R
# Input
## args: list, use ind.ioM
## phiM: matrix of parameters
## pres: vector of residual error parameters
## Dargs: list, use transform.par, modeltype, structural.model(), etype.exp, IdM, XM, yM
## DYF: a vector, modified locally but not exported

compute.LLy.old<-function(phiM,args,Dargs,DYF,pres) {
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  if (Dargs$modeltype=="structural"){
    gpred<-error(fpred,pres,Dargs$XM$ytype)
    DYF[args$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
  } else {
    DYF[args$ind.ioM]<- -fpred
  }
  # for(iout in 1:2) {
  #   print(summary(fpred[Dargs$XM$ytype==iout]))
  #   print(summary(gpred[Dargs$XM$ytype==iout]))
  # }
  U<-colSums(DYF)
  return(U)
}

if(FALSE) {
# Input
## phiM: matrix of parameters
## Dargs: list, here we use the elements transform.par, model (structural model), IdM, XM, yM, ind.ioM (a vector of 0/1 indicating which elements of the square matrix DYF have data), list of outcomes (type, error.model, error.function)
## error.parameters: list of the error parameters for the different outcomes in the model (NULL for non-continuous outcome)
## DYF: a matrix of dimension nmax=max(n_i) times (N*nb.chains); the column for subject i has the last nmax-n_i elements set to 0, so that colSums sums on the observations for that subject 
# Output 
## U: vector with either sum(-LL) (continuous models, minus the constant sum(log(1/sqrt(2*pi)))) or sum(-logpdf) (likelihood models)
compute.LLy<-function(phiM, Dargs, DYF, error.parameters) {
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$model(psiM,Dargs$IdM,Dargs$XM)
  lpred<-fpred
  for(iout in 1:length(Dargs$outcome)) {
    idx1<-which(Dargs$XM$ytype==iout)
    # print(summary(fpred[idx1]))
    if(Dargs$outcome[[iout]]@type=="continuous") {
      if(Dargs$outcome[[iout]]@error.model=="exponential") fpred[idx1]<-log(cutoff(fpred[idx1]))
      gpred<-Dargs$outcome[[iout]]@error.function(fpred[idx1], Dargs$outcome[[iout]]@error.parameters)
      lpred[idx1]<-0.5*((Dargs$yM[idx1]-fpred[idx1])/gpred)**2+log(gpred)
      # print(summary(gpred))
    } else lpred[idx1]<- (-fpred[idx1])
  }
  DYF[Dargs$ind.ioM]<-lpred
  U<-colSums(DYF)
  return(U)
}
} else # load final debugged version
  source(file.path(saemixDir,"Rext","func_aux.R"))

################################### Testing
# Creating input for compute.LLy.mat
## from initialise
# using several Markov chains
nchains<-2
chdat<-new(Class="SaemixRepData",data=pkpd.saemix, nb.chains=nchains)
NM<-chdat["NM"]
IdM<-chdat["dataM"]$IdM
yM<-chdat["dataM"]$yM
XM<-chdat["dataM"][,c(pkpd.saemix["name.predictors"],pkpd.saemix["name.cens"],pkpd.saemix["name.mdv"],pkpd.saemix["name.ytype"]),drop=FALSE]
io<-matrix(data=0,nrow=pkpd.saemix["N"],ncol=max(pkpd.saemix["nind.obs"]))
for(i in 1:pkpd.saemix@N)
  io[i,1:pkpd.saemix["nind.obs"][i]]<-1
ioM<-do.call(rbind,rep(list(io),nchains))
ind.ioM <- which(t(ioM)!=0)
DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])

phiM<-psiM<-do.call(rbind, rep(list(pkpd.psi0), pkpd.saemix@N*nchains)) # h=id() so phiM=psiM

# Dargs for initial version
Dargs1<-list(transform.par=c(0,0,0,0), IdM=IdM, XM=XM, yM=yM, modeltype="structural", etype.exp=integer(0), structural.model=pkpd.model@model) # normal distribution for parameters
pres1<-c(pkpd.model@outcome[[1]]@error.parameters,pkpd.model@outcome[[2]]@error.parameters,0) # old version: needs 2 parameters for each error model
xvec1<-compute.LLy.old(phiM,args=list(ind.ioM=ind.ioM),Dargs1,DYF,pres1)

# Dargs for new version
Dargs2<-list(transform.par=c(0,0,0,0), IdM=IdM, XM=XM, yM=yM, ind.ioM=ind.ioM, model=pkpd.model@model, outcome=pkpd.model@outcome)
pres2<-list(pkpd.model@outcome[[1]]@error.parameters, pkpd.model@outcome[[2]]@error.parameters)
xvec2<-compute.LLy(phiM, Dargs2, DYF)
if(max(abs(xvec1-xvec2))>0) cat("Mismatch between old and new versions\n")

### Final change
# error.parameters are already in Dargs$outcome, removed from compute.LLy

#### Eco 12/04/2022 - stopped here

###########################################################
