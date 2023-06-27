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

library(tidyverse)

# Whether to save the plots
saveFigs<-FALSE
figDir <- getwd()

# Number of bootstrap samples
runBootstrap <- FALSE # to read the results from disk
nboot <-10
# nboot <- 200

# Data
data(lung.saemix)
saemix.data<-saemixData(name.data=lung.saemix,header=TRUE,name.group=c("id"),
                        name.predictors=c("time","status","cens"),name.response=c("status"),
                        name.covariates=c("age", "sex", "ph.ecog", "ph.karno", "pat.karno", "wt.loss","meal.cal"),
                        units=list(x="days",y="",covariates=c("yr","","-","%","%","cal","pounds")), verbose=FALSE)

plotDiscreteData(saemix.data, outcome="tte", which.cov="sex")

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

#### Weibull model for TTE data

weibulltte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  lambda <- psi[id,1] # Parameters of the Weibull model
  beta <- psi[id,2]
  Nj <- length(T)
  
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  hazard <- (beta/lambda)*(T/lambda)^(beta-1) # ln(H')
  H <- (T/lambda)^beta # ln(H)
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))
  return(logpdf)
}

saemix.model<-saemixModel(model=weibulltte.model,description="time model",modeltype="likelihood",
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("lambda","beta"))),
                          transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
tte.fit<-saemix(saemix.model,saemix.data,saemix.options)
plot(tte.fit, plot.type="convergence")
summary(tte.fit)

# Model evaluation

simulateWeibullTTE <- function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  lambda <- psi[,1] # Parameters of the Weibull model
  beta <- psi[,2]
  Nj <- length(T)
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  tevent<-T
  Vj<-runif(dim(psi)[1])
  tsim<-lambda*(-log(Vj))^(1/beta) # nsuj events
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
psiM<-data.frame(lambda=rnorm(nsuj, mean=tte.fit@results@fixed.effects[1], sd=2), beta=tte.fit@results@fixed.effects[2])
id1<-rep(1:nsuj, each=2)
simtime<-simulateWeibullTTE(psiM, id1, xidep1)

par(mfrow=c(1,2))
hist(saemix.data@data$time[saemix.data@data$time>0], breaks=30, xlab="Time", main="Original data")
hist(simtime[simtime>0], breaks=30, xlim=c(0,1000), xlab="Time", main="Simulated data")

# Ignoring the cens column and assuming a common censoring time instead
simulateWeibullTTE.maxcens <- function(psi,id,xidep) {
  etime<-xidep[,1]
  censoringtime <- max(etime)
  lambda <- psi[,1]
  beta <- psi[,2]
  N<-dim(psi)[1]
  Vj<-runif(N)
  T<-lambda*(-log(Vj))^(1/beta)
  T[T>censoringtime]<-censoringtime
  etime[etime>0]<-T
  return(etime)
}
simtime.maxcens<-simulateWeibullTTE.maxcens(psiM, id1, xidep1)

# Compare the two simulation models
par(mfrow=c(1,3))
hist(saemix.data@data$time[saemix.data@data$time>0], breaks=30, xlab="Time", main="Original data")
hist(simtime[simtime>0], breaks=30, xlim=c(0,1000), xlab="Time", main="Simulated data")
hist(simtime.maxcens[simtime.maxcens>0], breaks=30, xlim=c(0,1000), xlab="Time", main="Simulated data")

ypred<-predict(tte.fit)

# Use survival package to assess Survival curve
xtim<-seq(0,max(lung.saemix$time), length.out=200)
estpar<-tte.fit@results@fixed.effects
estse<-tte.fit@results@se.fixed
ypred<-exp(-(xtim/estpar[1])^(estpar[2]))

# Computing SE for the survival curve based on linearised FIM (probably not a good idea) through the delta-method
invfim<-solve(tte.fit@results@fim[1:2,1:2])
xcal<- (xtim/estpar[1])^estpar[2]
dsdbeta<- -log(xtim/estpar[1]) * xcal *exp(-xcal)
dsdalpha<- estpar[2]/estpar[1] * xcal *exp(-xcal)
xmat<-rbind(dsdalpha, dsdbeta)
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

########################################################################
### RTTE model

# Simulating RTTE data by simulating from U(0,1) and inverting the cdf
simul.rtte.unif<-function(psi) { # xidep, id not important, we only use psi
  censoringtime <- 3
  maxevents <- 30
  lambda <- psi[,1]
  beta <- psi[,2]
  simdat<-NULL
  N<-nrow(psi)
  for(i in 1:N) {
    eventTimes<-c(0)
    T<-0
    Vj<-runif(1)
    #    T <- (-log(Vj)*lambda[i])^(beta[i])
    T<-lambda[i]*(-log(Vj))^(1/beta[i])
    nev<-0
    while (T < censoringtime & nev<maxevents){
      eventTimes <- c(eventTimes, T)  
      nev<-nev+1
      Vj<-runif(1)
      #      T <- T+(-log(Vj)*lambda[i])^(beta[i])
      #      T<-(-log(Vj)*lambda[i] + T^(1/beta[i]))^(beta[i])
      T<-lambda[i]*(-log(Vj) + (T/lambda[i])^(beta[i]))^(1/beta[i])
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
psiM<-data.frame(lambda=param[1]*exp(rnorm(nsuj,sd=omega[1])), beta=param[2]*exp(param[3]*risk+rnorm(nsuj,sd=omega[2])))
simdat <- simul.rtte.unif(psiM)
simdat$risk<-as.integer(simdat$id>(nsuj/2))

saemix.data<-saemixData(name.data=simdat, name.group=c("id"), name.predictors=c("T"), name.response="status", name.covariates="risk", verbose=FALSE)

rtte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  N <- nrow(psi) # nb of subjects
  Nj <- length(T) # nb of events (including 0 and censoring times)
  # censoringtime = 6
  censoringtime = max(T) # same censoring for everyone
  lambda <- psi[id,1]
  beta <- psi[id,2]
  tinit <- which(T==0) # indices of beginning of observation period
  tcens <- which(T==censoringtime) # indices of censored events 
  tevent <- setdiff(1:Nj, append(tinit,tcens)) # indices of non-censored event times
  hazard <- (beta/lambda)*(T/lambda)^(beta-1)
  H <- (T/lambda)^beta
  logpdf <- rep(0,Nj)
  logpdf[tcens] <- -H[tcens] + H[tcens-1]
  logpdf[tevent] <- -H[tevent] + H[tevent-1] + log(hazard[tevent])
  return(logpdf)
}

saemix.model.base<-saemixModel(model=rtte.model,description="Repeated TTE model",modeltype="likelihood",
                               psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("lambda","beta"))),
                               transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.model<-saemixModel(model=rtte.model,description="Repeated TTE model",modeltype="likelihood",
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("lambda","beta"))),
                          transform.par=c(1,1),covariate.model=matrix(c(0,1),ncol=2),
                          covariance.model=matrix(c(1,0,0,1),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, displayProgress=FALSE, print=FALSE)
rtte.fit<-saemix(saemix.model,saemix.data,saemix.options)
plot(rtte.fit, plot.type="convergence")
print(rtte.fit@results)

