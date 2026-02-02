# Loading libraries
library(xtable)
library(ggplot2)
library(gridExtra)

# Loading saemix
library(saemix)

# Folders
# workDir<-"/home/eco/work/saemix/saemixextension/paperSaemix3"
workDir <- getwd() # working directory should be paperSaemix3
# setwd(workDir)
figDir <- file.path(workDir, "figs")
nsim<-200 # Number of simulations
saveFigures <- FALSE

###################################################### Data exploration
data(lung.saemix)

# all covariates (but need to manage the missing covariates)
# ECOG status treated as categorical, creating  dummy covariates for ECOG=1 and ECOG=2 or 3 
### note: could do this with the transformation functions
# missing patient Karnofsky scores set to median (in 3 patients)
# other covariates still have missing values
lung2<-lung.saemix
lung2$ecog1<-ifelse(lung2$ph.ecog==1,1,0)
lung2$ecog23<-ifelse(lung2$ph.ecog>1,1,0)
lung2$pat.karno[is.na(lung2$pat.karno)]<-median(lung2$pat.karno, na.rm=TRUE)
lung.data<-saemixData(name.data=lung2,header=TRUE,name.group=c("id"),
                      name.predictors=c("time","status","cens"),name.response=c("status"),
                      name.covariates=c("age", "sex", "ecog1","ecog23", "ph.karno", "pat.karno","ph.ecog"),
                      units=list(x="days",y="",covariates=c("yr","","-","-","%","%")), verbose=FALSE)

plotDiscreteData(lung.data, outcome="tte", which.cov="sex")
xplot1<-plotDiscreteData(lung.data, outcome="tte", which.cov="sex")

# Stratifying on ECOG status
xplot2<-plotDiscreteData(lung.data, outcome="tte", which.cov="ph.ecog")

grid.arrange(grobs=list(xplot1, xplot2), nrow=1, ncol=2)

# Histogram
hist(lung.saemix$time[lung.saemix$status==1])

###################################################### Fit Weibull model

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

saemix.model<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                          transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
tte.fit<-saemix(saemix.model,lung.data,saemix.options)

###################################################### Fit alternative TTE models
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
exptte.fit<-saemix(saemix.model.exp,lung.data,saemix.options)

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
gomptte.fit<-saemix(saemix.model.gomp,lung.data,saemix.options)

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
gammatte.fit<-try(saemix(saemix.model.gamma,lung.data,saemix.options))

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
logistte.fit<-saemix(saemix.model.logis,lung.data,saemix.options)

###################################################### Table with the results, diagnostic plots for the two best fitting models

# Table
restte<-data.frame(Model=c("Exponential","Weibull","Gompertz","Gamma","Log-logistic"), 
                   BIC=c(BIC(exptte.fit),BIC(tte.fit), BIC(gomptte.fit), BIC(gammatte.fit), BIC(logistte.fit)))
print(xtable(restte), only.contents=TRUE, include.rownames=FALSE,  floating=F, sanitize.rownames.function = identity)

# Comparing KM plots for Gompertz and Weibull

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
tte.fit@model@simulate.function <- simulateWeibullTTE
simtte.fit <- simulateDiscreteSaemix(tte.fit, nsim=500)

gpl <- discreteVPC(simtte.fit, outcome="TTE")

## Gompertz simulations
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

gomptte.fit@model@simulate.function<-simulateGompertzTTE
simgomptte.fit <- simulateDiscreteSaemix(gomptte.fit, nsim=500)

gpl2 <- discreteVPC(simgomptte.fit, outcome="TTE", main="Gompertz")

if(saveFigs) {
  namfig<-"lung_compareTTEfits.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  grid.arrange(gpl,gpl2, nrow=1)
  dev.off()
}

###################################################### Covariate model

# Stepwise procedure - selects  sex and ECOG 2-3, with IIV on gamma
# externalised to lung_stepSelection.R

if(FALSE)
  covtte.fit <- step.saemix(tte.fit, direction="both")

covmodel <- cbind(c(0,1,0,1,0,0),rep(0,6))
weibull.model.cov<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                          psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                          transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE), 
                          covariate.model=covmodel, verbose=FALSE)
weibull.model.cov2<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                               psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                               transform.par=c(1,1),covariance.model=matrix(c(0,0,0,1),ncol=2, byrow=TRUE), 
                               covariate.model=covmodel, verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
weibcov.fit<-saemix(weibull.model.cov,lung.data,saemix.options)
weibcov.fit2<-saemix(weibull.model.cov2,lung.data,saemix.options)

print(weibcov.fit)
print(weibcov.fit2)

# large SE on omega(gamma) but not sure we can actually estimate IIV with only 1 point per subject
# effects of sex and ECOG have p-values that are actually NS with a Wald test (but not sure about how accurate the FIM is for TTE models)

covmodelsex <- cbind(c(0,1,0,0,0,0),rep(0,6))
covmodelECOG <- cbind(c(0,0,0,1,0,0),rep(0,6))

weibull.model.sex<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                                psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                                transform.par=c(1,1),covariance.model=matrix(c(0,0,0,1),ncol=2, byrow=TRUE), 
                                covariate.model=covmodelsex, verbose=FALSE)
weibull.model.ECOG<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                               psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                               transform.par=c(1,1),covariance.model=matrix(c(0,0,0,1),ncol=2, byrow=TRUE), 
                               covariate.model=covmodelECOG, verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
weibcov.fit.sex<-saemix(weibull.model.sex,lung.data,saemix.options)
print(weibcov.fit.sex)

weibcov.fit.ECOG<-saemix(weibull.model.ECOG,lung.data,saemix.options)
print(weibcov.fit.ECOG)

# Diagnostic plots 
nsim<-200
weibcov.fit2@model@simulate.function <- simulateWeibullTTE
ysimTTE<-simulateDiscreteSaemix(weibcov.fit2, nsim=nsim)
discreteVPC(ysimTTE)


###################################################### Bootstrap SE

# Final model with covariates - case bootstrap
tab1 <- read.table(file.path(workDir,"covselect", "bootstrapCase_weibullTTEcov.res"), header=F, skip=1)
colnames(tab1)<-c("irep",weibcov.fit2@results@name.fixed,weibcov.fit2@results@name.random)
tab1[,6]<-sqrt(tab1[,6]) # sd instead of var
apply(tab1,2,sd)
x1<-apply(tab1,2,quantile,c(0.025,0.975))
yfit<-weibcov.fit2
df2<-data.frame(Parameter=c(colnames(tab1)[2:5],"omega.gamma"), Estimate=c(yfit@results@fixed.effects, sqrt(yfit@results@omega[2,2])), SDboot=apply(tab1,2,sd)[2:6])
df2<-cbind(df2, t(x1[,2:6]))

par(mfrow=c(2,3))
for(i in 2:6) {
  xvec<-tab1[,i]
  limx <- quantile(xvec,c(0.025,0.975))
  xvec<-xvec[xvec>=limx[1] & xvec<=limx[2]]
  hist(xvec, main=colnames(tab1)[i], xlab=colnames(tab1)[i], xlim=limx, breaks=50)
}

# No covariates - less variability
tab1<-read.table(file.path(workDir, "boostrapRes", "bootstrapCase_weibullTTE.res"), header=F, skip=1)

apply(tab1,2,sd)
x1<-apply(tab1,2,quantile,c(0.025,0.975))
yfit<-tte.fit
df<-data.frame(Parameter=c("Te","gamma"), Estimate=yfit@results@fixed.effects, SDboot=apply(tab1,2,sd)[2:3])
df<-cbind(df, t(x1[,2:3]))

# Covariate model, ECOG as continuous
tab1<-read.table(file.path(workDir, "boostrapRes", "bootstrapCase_weibullTTEcont.res"), header=F, skip=1)
apply(tab1,2,sd)
x1<-apply(tab1,2,quantile,c(0.025,0.975))
yfit<-weibull.fit.cov

df<-data.frame(Parameter=yfit@results@conf.int[-c(5),1], Estimate=yfit@results@conf.int[-c(5),2], SDboot=apply(tab1,2,sd)[-c(1)])
df<-cbind(df, t(x1[,-c(1)]))
print(xtable(df), include.rownames = FALSE)

# Covariate model, ECOG as categorical
weibcov.fit2@model@simulate.function<-simulateWeibullTTE

# TODO debug why not working
# cond.TTE <- saemix.bootstrap(weibcov.fit2, nboot=2)

tab1<-read.table(file.path(workDir,"covselect",  "bootstrapCase_weibullTTEcov.res"), header=F, skip=1)

apply(tab1,2,sd)
x1<-apply(tab1,2,quantile,c(0.025,0.975))
yfit<-weibcov.fit2
df2<-data.frame(Parameter=c("Te","beta_sex","beta_e1","beta_e23","gamma"), Estimate=yfit@results@fixed.effects, SDboot=apply(tab1,2,sd)[2:6])
df2<-cbind(df2, t(x1[,2:6]))



###################################################### Misc
############# 
# Log-logistic- alternate parameterisation
logis.model2<-function(psi,id,xidep) {
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
  # alternate parameterisation, same results
  #  hazard <- (gamma/Te)*(T/Te)^(gamma-1) /(1+(T/Te)^gamma)^2 # H'
  #  H <- 1/(1+(Te/T)^gamma) # H= -ln(S)
  
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))=ln(S)=-H
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))=ln(S)+ln(h)
  return(logpdf)
}

saemix.model.logis2<-saemixModel(model=logis.model2,description="Log-logistic TTE model",modeltype="likelihood",
                                psi0=matrix(c(200,1),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                                transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
try.fit<-saemix(saemix.model.logis,lung.data,saemix.options)
print(try.fit)

############# Not run
# Putting variability on a dummy variable

dumweibtte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  Te <- psi[id,1] # Parameters of the Weibull model
  gamma <- psi[id,2]
  dummy <-psi[id,3]
  Nj <- length(T)
  
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  hazard <- (gamma/Te)*(T/Te)^(gamma-1) # h
  H <- (T/Te)^gamma # H=ln(S)
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))=ln(S)=-H
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))=ln(S)+ln(h)
  return(logpdf)
}

# Currently not possible to put IIV on a (fixed) dummy variable, need to change this (see Alexandra's code)
if(FALSE) {
  saemix.model.dummy<-saemixModel(model=dumweibtte.model,description="Weibull TTE model",modeltype="likelihood",
                                  psi0=matrix(c(1,2,1),ncol=3,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma","dummy"))),
                                  transform.par=c(1,1,1),fixed.estim=c(1,1,0),covariance.model=diag(c(0,0,1)), verbose=FALSE)
  saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
  dumtte.fit<-saemix(saemix.model.dummy,lung.data,saemix.options)
  print(dumtte.fit)
}
