# Folders
workDir<-getwd() 

# @Eco
workDir<-"/home/eco/work/saemix/saemixextension/paperSaemix3"
saemixDir <- "/home/eco/work/saemix/saemixextension"
setwd(workDir)

# Libraries
library(saemix)

# Library survival for KM plot
library(survival)

library(gridExtra)
library(tidyverse)

# Whether to save the plots
saveFigs<-FALSE
figDir <- getwd()

# Number of bootstrap samples
runBootstrap <- FALSE # to read the results from disk
nboot <-10
# nboot <- 200

# Covariate model building (stepwise algorithm, takes a while to run)
runCovariateModelSelection<-FALSE

######################################################################## Data
# Lung cancer data

data(lung.saemix)
# all covariates (but need to manage the missing covariates)
# ECOG status treated as continuous
# missing patient Karnofsky scores set to median (in 3 patients)
# other covariates still have missing values
lung1<-lung.saemix
lung1$pat.karno[is.na(lung1$pat.karno)]<-median(lung1$pat.karno, na.rm=TRUE)

saemix.data.contPH<-saemixData(name.data=lung1,header=TRUE,name.group=c("id"),
                               name.predictors=c("time","status","cens"),name.response=c("status"),
                               name.covariates=c( "sex", "ph.ecog", "ph.karno", "pat.karno", "age"),
                               units=list(x="days",y="",covariates=c("","-","%","%","yr")), verbose=FALSE)

plotDiscreteData(saemix.data, outcome="tte", which.cov="sex")

xplot1<-plotDiscreteData(saemix.data, outcome="tte", which.cov="sex")
xplot2<-plotDiscreteData(saemix.data, outcome="tte", which.cov="ph.ecog")
grid.arrange(grobs=list(xplot1, xplot2), nrow=1, ncol=2)

# Histogram
hist(lung.saemix$time[lung.saemix$status==1])

# Note: missing data in pat.karno, wt.loss and meal.cal
if(FALSE)
  print(summary(lung.saemix))

##### Kaplan-Meier plot
lung.surv<-lung.saemix[lung.saemix$time>0,]
lung.surv$status<-lung.surv$status+1
Surv(lung.surv$time, lung.surv$status) # 1=censored, 2=dead
nonpar.fit <- survfit(Surv(time, status) ~ 1, data = lung.surv)
plot(nonpar.fit)

######################################################################## Weibull model for TTE data
# Model function
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
  H <- (T/Te)^gamma # H= -ln(S)
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))=ln(S)=-H
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))=ln(S)+ln(h)
  return(logpdf)
}

saemix.model<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                          transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)

# Fitting in saemix
tte.fit<-saemix(saemix.model,saemix.data,saemix.options)
plot(tte.fit, plot.type="convergence")
print(tte.fit)

######################################################################## Model evaluation
# Simulation function
simulateWeibullTTE <- function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  Te <- psi[,1] # Parameters of the Weibull model
  gamma <- psi[,2]
  Nj <- length(T)
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  tevent<-T
  Vj<-runif(dim(psi)[1])
  tsim<-Te*(-log(Vj))^(1/gamma) # nsuj events
  tevent[T>0]<-tsim
  tevent[tevent[cens]>T[cens]] <- T[tevent[cens]>T[cens]]
  return(tevent)
}
# Adding simulation function to the model object in the fit
tte.fit@model@simulate.function<-simulateWeibullTTE

### Simulations
simtte.fit <- simulateDiscreteSaemix(tte.fit, nsim=500)

# VPC (KM-VPC)
gpl <- discreteVPC(simtte.fit, outcome="TTE")
plot(gpl)

# Checking the simulation function
xidep1<-saemix.data@data[,saemix.data@name.predictors]
nsuj<-saemix.data@N
psiM<-data.frame(Te=rnorm(nsuj, mean=tte.fit@results@fixed.effects[1], sd=2), gamma=tte.fit@results@fixed.effects[2])
id1<-rep(1:nsuj, each=2)
simtime<-simulateWeibullTTE(psiM, id1, xidep1)

par(mfrow=c(1,2))
hist(saemix.data@data$time[saemix.data@data$time>0], breaks=30, xlim=c(0,1050),xlab="Time", main="Original data")
hist(simtime[simtime>0], breaks=30, xlim=c(0,1050), xlab="Time", main="Simulated data")

# Ignoring the cens column and assuming a common censoring time instead
simulateWeibullTTE.maxcens <- function(psi,id,xidep) {
  etime<-xidep[,1]
  censoringtime <- max(etime)
  Te <- psi[,1]
  gamma <- psi[,2]
  N<-dim(psi)[1]
  Vj<-runif(N)
  T<-Te*(-log(Vj))^(1/gamma)
  T[T>censoringtime]<-censoringtime
  etime[etime>0]<-T
  return(etime)
}
simtime.maxcens<-simulateWeibullTTE.maxcens(psiM, id1, xidep1)

# Compare the two simulation models
par(mfrow=c(1,3))
hist(saemix.data@data$time[saemix.data@data$time>0], breaks=30, xlim=c(0,1050), xlab="Time", main="Original data")
hist(simtime[simtime>0], breaks=30, xlim=c(0,1050), xlab="Time", main="Simulated data")
hist(simtime.maxcens[simtime.maxcens>0], breaks=30, xlim=c(0,1050), xlab="Time", main="Simulated data")

ypred<-predict(tte.fit)

#####################################################
# Use survival package to assess Survival curve 
xtim<-seq(0,max(lung.saemix$time), length.out=200)
estpar<-tte.fit@results@fixed.effects
estse<-tte.fit@results@se.fixed
ypred<-exp(-(xtim/estpar[1])^(estpar[2]))

######################### Not recommended
# Computing SE for the survival curve based on linearised FIM (probably not a good idea) through the delta-method
invfim<-solve(tte.fit@results@fim[1:2,1:2])
xcal<- (xtim/estpar[1])^estpar[2]
dsdgamma<- -log(xtim/estpar[1]) * xcal *exp(-xcal)
dsdalpha<- estpar[2]/estpar[1] * xcal *exp(-xcal)
xmat<-rbind(dsdalpha, dsdgamma)
#    x1<-t(xmat[,1:3]) %*% invfim %*% xmat[,1:3]
sesurv<-rep(0,length(xcal))
for(i in 1:length(xcal))
  sesurv[i]<-sqrt(t(xmat[,i]) %*% invfim %*% xmat[,i])

# Comparison between KM and parametric fit
par(mfrow=c(1,1))
plot(nonpar.fit, xlab = "Days", ylab = "Overall survival probability")
lines(xtim,ypred, col="red",lwd=2)
lines(xtim,ypred+1.96*sesurv, col="red",lwd=1, lty=2)
lines(xtim,ypred-1.96*sesurv, col="red",lwd=1, lty=2)

######################################################################## Comparing models
# Exponential
exptte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  Te <- psi[id,1] # Parameters of the Weibull model
  Nj <- length(T)
  
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  hazard <- (1/Te) # H'
  H <- (T/Te) #  H= -ln(S)
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))=ln(S)=-H
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))=ln(S)+ln(h)
  
  return(logpdf)
}

saemix.model.exp<-saemixModel(model=exptte.model,description="Exponential TTE model",modeltype="likelihood",
                              psi0=matrix(c(1),ncol=1,byrow=TRUE,dimnames=list(NULL,  c("Te"))),
                              transform.par=c(1),covariance.model=matrix(c(1),ncol=1, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
exptte.fit<-saemix(saemix.model.exp,saemix.data,saemix.options)
plot(exptte.fit, plot.type="convergence")
print(exptte.fit)

# Gompertz

gomptte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  Te <- psi[id,1] # Parameters of the Weibull model
  gamma <- psi[id,2]
  teprim <- Te/log(1+log(2)/gamma)
  Nj <- length(T)
  
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  hazard <- (gamma/teprim)*exp(T/teprim) # h
  H <- gamma*(exp(T/teprim)-1) # H
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))=ln(S)=-H
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))=ln(S)+ln(h)
  return(logpdf)
}
saemix.model.gomp<-saemixModel(model=gomptte.model,description="Gompertz TTE model",modeltype="likelihood",
                               psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                               transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
gomptte.fit<-saemix(saemix.model.gomp,saemix.data,saemix.options)
plot(gomptte.fit, plot.type="convergence")
print(gomptte.fit)

# Gamma
# incomplete gamma function for (x,a) : gamma(a) * pgamma(x, a, 1, lower = FALSE)

gammatte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  Te <- psi[id,1] # Parameters of the Weibull model
  lambda <- psi[id,2]
  Nj <- length(T)
  
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  #  hazard <- (lambda/Te) * (lambda*T/Te)^(lambda-1) * exp(lambda*T/Te) / (gamma(lambda) - pgamma(lambda*T/Te,lambda, 1, lower=FALSE))
  hazard <- (T/Te)^(lambda-1) * exp(-T/Te) /gamma(lambda) / Te
  #  H <- pgamma(T/Te, lambda, 1, lower=FALSE) / gamma(lambda)
  H <- pgamma(T/Te, lambda) # incomplete gamma gammainc(x,a)=pgamma(x,a)*gamma(a) and H=gammainc(T/Te, lambda) / gamma(lambda)
  #  H <- (1-pgamma(T/Te, lambda,1, lower=FALSE))
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))
  return(logpdf)
}
saemix.model.gamma<-saemixModel(model=gammatte.model,description="Gamma TTE model",modeltype="likelihood",
                                psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","k"))),
                                transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
gammatte.fit<-try(saemix(saemix.model.gamma,saemix.data,saemix.options))
plot(gammatte.fit, plot.type="convergence")
print(gammatte.fit)


# Log-logistic 
logis.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  Te <- psi[id,1] # Parameters of the Weibull model
  gamma <- psi[id,2]
  Nj <- length(T)
  
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  hazard <- (gamma/Te)*(T/Te)^(gamma-1) /(1+(T/Te)^gamma) # H'
  H <- log(1+(T/Te)^gamma) # H= -ln(S)
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))=ln(S)=-H
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))=ln(S)+ln(h)
  return(logpdf)
}

saemix.model.logis<-saemixModel(model=logis.model,description="Log-logistic TTE model",modeltype="likelihood",
                                psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                                transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
logistte.fit<-saemix(saemix.model.logis,saemix.data,saemix.options)
plot(logistte.fit, plot.type="convergence")
print(logistte.fit)

# Table comparing the models
restte<-data.frame(Model=c("Exponential","Weibull","Gompertz","Gamma","Log-logistic"), 
                   BIC=c(BIC(exptte.fit),BIC(tte.fit), BIC(gomptte.fit), BIC(gammatte.fit), BIC(logistte.fit)))
print(restte)

# Model evaluation comparing Weibull and Gompertz (best fitting models)
# Simulate events based on the observed individual censoring time
simulateGompertzTTE <- function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  delta <- xidep[,3] # censoring indicator
  cens<-which(delta==1) # censoring times (subject specific)
  tmax <- max(T[cens]) # maximum censoring time observed in dataset
  init <- which(T==0)
  Te <- psi[,1] # Parameters of the Weibull model
  gamma <- psi[,2]
  teprim <- Te/log(1+log(2)/gamma)
  
  Nj <- length(T)
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  tevent<-T
  Vj<-runif(dim(psi)[1])
  tsim<-teprim*log(1-log(Vj)/gamma) #   events
  tevent[T>0]<-tsim
  tevent[delta==1 & tevent>T] <- T[delta==1 & tevent>T] # subject-specific censoring time
  #  tevent[delta==0 & tevent>tmax] <- tmax # censoring to tmax (for subjects who experienced an event)
  #  tevent[tevent[dead]>tmax] <- tmax # for subjects who initially experienced the event, use maximal censoring time
  return(tevent)
}

# Checking the simulation function
xidep1<-saemix.data@data[,saemix.data@name.predictors]
nsuj<-saemix.data@N
psiM<-data.frame(Te=rnorm(nsuj, mean=gomptte.fit@results@fixed.effects[1], sd=2), gamma=gomptte.fit@results@fixed.effects[2])
id1<-rep(1:nsuj, each=2)
simtime<-simulateGompertzTTE(psiM, id1, xidep1)

par(mfrow=c(1,2))
hist(saemix.data@data$time[saemix.data@data$time>0], breaks=30, xlim=c(0,1050),xlab="Time", main="Original data")
hist(simtime[simtime>0], breaks=30, xlim=c(0,1050), xlab="Time", main="Simulated data")

gomptte.fit@model@simulate.function<-simulateGompertzTTE
simgomptte.fit <- simulateDiscreteSaemix(gomptte.fit, nsim=500)

gpl2 <- discreteVPC(simgomptte.fit, outcome="TTE")
plot(gpl2)

grid.arrange(gpl,gpl2, nrow=1)

######################################################################## Covariate model
# Stepwise covariate model using the Weibull model

# Toggle to TRUE to run (takes a while)
if(runCovariateModelSelection)
  covtte.fit <- step.saemix(tte.fit, direction="both")

# Covariate model (final step of covtte.fit above)

# Covariate model with only sex and ECOG score

######################################################################## Bootstrap SE **TODO** see fig
# 

######################################################################## RTTE model

# Simulating RTTE data by simulating from U(0,1) and inverting the cdf
simul.rtte.unif<-function(psi) { # xidep, id not important, we only use psi
  censoringtime <- 3
  maxevents <- 30
  Te <- psi[,1]
  gamma <- psi[,2]
  simdat<-NULL
  N<-nrow(psi)
  for(i in 1:N) {
    eventTimes<-c(0)
    T<-0
    Vj<-runif(1)
    #    T <- (-log(Vj)*Te[i])^(gamma[i])
    T<-Te[i]*(-log(Vj))^(1/gamma[i])
    nev<-0
    while (T < censoringtime & nev<maxevents){
      eventTimes <- c(eventTimes, T)  
      nev<-nev+1
      Vj<-runif(1)
      #      T <- T+(-log(Vj)*Te[i])^(gamma[i])
      #      T<-(-log(Vj)*Te[i] + T^(1/gamma[i]))^(gamma[i])
      T<-Te[i]*(-log(Vj) + (T/Te[i])^(gamma[i]))^(1/gamma[i])
    }
    if(nev==maxevents) {
      message("Reached maximum number of events\n")
    }
    eventTimes<-c(eventTimes, censoringtime)
    cens<-rep(1,length(eventTimes))
    cens[1]<-cens[length(cens)]<-0
    simdat<-rbind(simdat,
                  data.frame(id=i, T=eventTimes, status=cens))
  }
  return(simdat)
}

# Subjects
set.seed(12345)
param<-c(2, 1.5, 0.5)
# param<-c(4, 1.2, 0.3)
omega<-c(0.25,0.25)
nsuj<-200
risk<-rep(0,nsuj)
risk[(nsuj/2+1):nsuj]<-1
psiM<-data.frame(Te=param[1]*exp(rnorm(nsuj,sd=omega[1])), gamma=param[2]*exp(param[3]*risk+rnorm(nsuj,sd=omega[2])))
simdat <- simul.rtte.unif(psiM)
simdat$risk<-as.integer(simdat$id>(nsuj/2))

saemix.data<-saemixData(name.data=simdat, name.group=c("id"), name.predictors=c("T"), name.response="status", name.covariates="risk", verbose=FALSE)

rtte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  N <- nrow(psi) # nb of subjects
  Nj <- length(T) # nb of events (including 0 and censoring times)
  # censoringtime = 6
  censoringtime = max(T) # same censoring for everyone
  Te <- psi[id,1]
  gamma <- psi[id,2]
  tinit <- which(T==0) # indices of beginning of observation period
  tcens <- which(T==censoringtime) # indices of censored events 
  tevent <- setdiff(1:Nj, append(tinit,tcens)) # indices of non-censored event times
  hazard <- (gamma/Te)*(T/Te)^(gamma-1)
  H <- (T/Te)^gamma
  logpdf <- rep(0,Nj)
  logpdf[tcens] <- -H[tcens] + H[tcens-1]
  logpdf[tevent] <- -H[tevent] + H[tevent-1] + log(hazard[tevent])
  return(logpdf)
}

saemix.model.base<-saemixModel(model=rtte.model,description="Repeated TTE model",modeltype="likelihood",
                               psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                               transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.model<-saemixModel(model=rtte.model,description="Repeated TTE model",modeltype="likelihood",
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                          transform.par=c(1,1),covariate.model=matrix(c(0,1),ncol=2),
                          covariance.model=matrix(c(1,0,0,1),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, displayProgress=FALSE, print=FALSE)
rtte.fit<-saemix(saemix.model,saemix.data,saemix.options)
plot(rtte.fit, plot.type="convergence")
print(rtte.fit@results)

