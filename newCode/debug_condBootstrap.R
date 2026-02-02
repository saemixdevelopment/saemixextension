# Debug conditional bootstrap for TTE models and others...

# Loading libraries
library(xtable)
library(ggplot2)

# Loading saemix
library(saemix)

# Folders
saemixDir<-"/home/eco/work/saemix/saemixextension"
workDir<-file.path(saemixDir, "paperSaemix3")
# setwd(workDir)
figDir <- file.path(workDir, "figs")
nsim<-200 # Number of simulations

###################################################### Data
data(lung.saemix)

# Managing covariates
## create dummy covariates for ECOG=1 and ECOG=2 or 3
## impute missing pat.karno to the median
lung2<-lung.saemix
lung2$ecog1<-ifelse(lung2$ph.ecog==1,1,0)
lung2$ecog23<-ifelse(lung2$ph.ecog>1,1,0)
lung2$pat.karno[is.na(lung2$pat.karno)]<-median(lung2$pat.karno, na.rm=TRUE)

saemix.data<-saemixData(name.data=lung2,header=TRUE,name.group=c("id"),
                        name.predictors=c("time","status","cens"),name.response=c("status"),
                        name.covariates=c("age", "sex", "ecog1","ecog23", "ph.karno", "pat.karno"),
                        units=list(x="days",y="",covariates=c("yr","","-","-","%","%")), verbose=FALSE)

###################################################### Model

weibulltte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  Te <- psi[id,1] # Parameters of the Weibull model
  gamma <- psi[id,2]
  Nj <- length(T)
  
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  hazard <- (gamma/Te)*(T/Te)^(gamma-1) # h
  H <- (T/Te)^gamma # H=ln(S)
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))=ln(S)=-H
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))=ln(S)+ln(h)
  return(logpdf)
}

simulateWeibullTTE <- function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  delta <- xidep[,3] # censoring indicator
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  tmax <- max(T[cens]) # maximum censoring time observed in dataset
  init <- which(T==0)
  Te <- psi[,1] # Parameters of the Weibull model
  gamma <- psi[,2]
  Nj <- length(T)
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  tevent<-T
  Vj<-runif(dim(psi)[1])
  tsim<-Te*(-log(Vj))^(1/gamma) #   events
  tevent[T>0]<-tsim
  tevent[delta==1 & tevent>T] <- T[delta==1 & tevent>T] # subject-specific censoring time
  #  tevent[delta==0 & tevent>tmax] <- tmax # censoring to tmax (for subjects who experienced an event)
  #  tevent[tevent[dead]>tmax] <- tmax # for subjects who initially experienced the event, use maximal censoring time
  return(tevent)
}


saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)

covmodel <- cbind(c(0,1,0,1,0,0),rep(0,6))
weibull.model.cov2<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                                simulate.function = simulateWeibullTTE,
                                psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                                transform.par=c(1,1),covariance.model=matrix(c(0,0,0,1),ncol=2, byrow=TRUE), 
                                covariate.model=covmodel, verbose=FALSE)

###################################################### Fit

weibcov.fit2<-saemix(weibull.model.cov2,saemix.data,saemix.options)

nboot<-2
# essai bootstrap - case
boot.tab <- saemix.bootstrap(weibcov.fit2, method="case", nboot=nboot)

# essai bootstrap - conditional
boot.tab <- saemix.bootstrap(weibcov.fit2, method="conditional", nboot=nboot)


###################################################### Debug bootstrap
weibcov.fit2@model@simulate.function<-simulateWeibullTTE
saemixObject <- weibcov.fit2
nboot<-2
saemix.options<-saemixObject["options"]
saemix.options$directory<-"current"
saemix.options$fix.seed<-FALSE
saemix.options$map<-FALSE   # Only parameter estimates are required for bootstrap
saemix.options$fim<-FALSE
saemix.options$displayProgress<-FALSE 
saemix.options$save.graphs<-FALSE
saemix.options$save<-FALSE
saemix.options$ll.is<-FALSE
saemix.options$print<-FALSE

if(saemixObject@model@modeltype=="structural") idx.eps<-saemixObject@model@indx.res else idx.eps<-integer(0)
idx.iiv<-saemixObject@model@indx.omega
idx.rho<-which(saemixObject@model@covariance.model[lower.tri(saemixObject@model@covariance.model)]==1)
bootstrap.distribution<-failed.runs<-data.frame()
nelements <- length(saemixObject@results@fixed.effects)+length(idx.iiv)+length(idx.rho)+length(idx.eps)
# Starting point: estimates from the fit 
model.boot<-saemixObject["model"]
model.boot@psi0 <- model.boot["betaest.model"]
model.boot@psi0[model.boot["betaest.model"]==1]<-saemixObject@results@fixed.effects

###################################################### Case bootstrap
method<-"case"

# Problem with covariate names, for some reason the wrong covariates are used
saemixObject@data@name.covariates <- rownames(model.boot@covariate.model)

# A bootstrap sample - case bootstrap
data.boot <- dataGen.case(saemixObject)
par(mfrow=c(2,2))
hist(saemixObject@data@data$time[saemixObject@data@data$time>0], breaks=50)
hist(data.boot@data$time[data.boot@data$time>0], breaks=50)

# Fit
fit.boot<-try(saemix(model.boot, data.boot, saemix.options))

###################################################### Conditional bootstrap
# ToDo
method<-"conditional"
nboot<-2
nsamp<-100
saemix.options <- saemixObject@options

# saemix.predict: ok
if(method=="residual" | method=="conditional") {
  saemixObject<-saemix.predict(saemixObject) # estimate individual parameters and compute residuals (currently iwres are needed also for conditional but need to modify this in a further extension to cNP ECO TODO)
}
if(method=="conditional") {
  ndone <- dim(saemixObject@results@phi.samp)
  if(!is.null(ndone)) ndone<-ndone[3] else ndone<-0
  if(ndone<nsamp) {
    if(saemixObject@options$warnings) message("Not enough samples in the object, sampling from the conditional distribution\n")
    saemixObject<-conddist.saemix(saemixObject, nsamp=nsamp) # estimate conditional distributions and sample residuals
  }
  eta.sampc<-centerDist.NPcond(saemixObject, nsamp=nsamp) # Center eta samples from the conditional distribution, to avoid doing this repeatedly
}

# A bootstrap sample - conditional bootstrap
data.boot <- dataGen.NP(saemixObject, nsamp=nsamp,eta.sampc=eta.sampc, conditional=TRUE)
fit.boot<-try(saemix(model.boot, data.boot, saemix.options))

# Checking
