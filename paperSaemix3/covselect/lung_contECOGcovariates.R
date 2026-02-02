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
saveFigures <- FALSE

###################################################### Data
data(lung.saemix)

lung1<-lung.saemix
lung1$pat.karno[is.na(lung1$pat.karno)]<-median(lung1$pat.karno, na.rm=TRUE)

saemix.data.contPH<-saemixData(name.data=lung1,header=TRUE,name.group=c("id"),
                               name.predictors=c("time","status","cens"),name.response=c("status"),
                               name.covariates=c( "sex", "ph.ecog", "ph.karno", "pat.karno", "age"),
                               units=list(x="days",y="",covariates=c("","-","%","%","yr")), verbose=FALSE)

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

###################################################### Fit Weibull model

# Weibull model
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

## Weibull simulations
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

saemix.model<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                          transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
tte.fit<-saemix(saemix.model,saemix.data,saemix.options)

###################################################### Covariate selection (ECOG as continuous variable)

# Stepwise procedure considering ECOG PH as continuous
saemix.model<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                          simulate.function = simulateWeibullTTE,
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                          transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
weibull.fit.cont<-saemix(saemix.model,saemix.data.contPH,saemix.options)
print(weibull.fit.cont)

covtte.fit.contPH <- step.saemix(weibull.fit.cont, direction="both")

covmodelcont <- cbind(c(1,1,0,1,0),c(0,0,1,0,0))
weibull.model.cont<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                                psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                                transform.par=c(1,1),covariance.model=matrix(c(0,0,0,1),ncol=2, byrow=TRUE), 
                                covariate.model=covmodelcont, verbose=FALSE)
weibull.fit.cov<-saemix(weibull.model.cont,saemix.data.contPH,saemix.options)
print(weibull.fit.cov)

## removing patient Karno on Te and ph Karno on gamma, too small
covmodelcont <- cbind(c(1,1,0,0,0),rep(0,5))
weibull.model.cont<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                                psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                                transform.par=c(1,1),covariance.model=matrix(c(0,0,0,1),ncol=2, byrow=TRUE), 
                                covariate.model=covmodelcont, verbose=FALSE)
weibull.fit.cov<-saemix(weibull.model.cont,saemix.data.contPH,saemix.options)
print(weibull.fit.cov)

# same but IIV on Te
weibull.model.cont2<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                                 psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                                 transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE), 
                                 covariate.model=covmodelcont, verbose=FALSE)
weibull.fit.cov2<-saemix(weibull.model.cont2,saemix.data.contPH,saemix.options)
print(weibull.fit.cov2)
