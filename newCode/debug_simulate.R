saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
datDir<-file.path(saemixDir,"data")

library(saemix)
library(ggplot2)

###############################################
# Theophylline example
theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T)
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
                        name.group=c("Id"),name.predictors=c("Dose","Time"),
                        name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
                        units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

# Model with covariate Weight

saemix.model<-saemixModel(model=model1cpt,modeltype="structural",
                          description="One-compartment model with first-order absorption",
                          psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
                          transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
theo.fit1<-saemix(saemix.model,saemix.data,saemix.options)

# Simulation
# Simulate individual parameters
simpsi <- simulateIndividualParameters(theo.fit1, nsim=2)
summary(simpsi)
theo.fit1@results@fixed.effects[1:3] # compare to median, should be close

# Simulate data
simdat<-simulate(theo.fit1, nsim=2)
expect_equal(dim(theo.fit1@sim.data@data)[1],0)
expect_equal(dim(simdat@sim.data@datasim)[1],theo.fit1@data@ntot.obs*2)
par(mfrow=c(1,3))
plot(simdat@data, new=FALSE)
plot(simdat@sim.data, irep=1, new=FALSE)
plot(simdat@sim.data, irep=2, new=FALSE)


################################################ Binary data example 
# toenail data
toenail.saemix<-read.table(file.path(datDir, "toenail.saemix.tab"), header=TRUE)

saemix.data<-saemixData(name.data=toenail.saemix,name.group=c("id"),name.predictors=c("time","y"), name.response="y",
                        name.covariates=c("treatment"),name.X=c("time"))

# saemix model
binary.model<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  #  logpdf<-rep(0,length(tim))
  P.obs = (y==0)*(1-pevent)+(y==1)*pevent
  logpdf <- log(P.obs)
  return(logpdf)
}
simulBinary<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-1/(1+exp(-logit))
  ysim<-rbinom(length(tim),size=1, prob=pevent)
  return(ysim)
}
saemix.model<-saemixModel(model=binary.model,description="Binary model",simulate.function=simulBinary,
                          modeltype="likelihood",
                          psi0=matrix(c(0,-.5,0,0.5),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                          transform.par=c(0,0), covariate.model=c(0,1),covariance.model=matrix(c(1,0,0,1),ncol=2))

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE)

# saemix fit
binary.fit<-saemix(saemix.model,saemix.data,saemix.options)

# Simulation
# Simulate individual parameters
simpsi <- simulateIndividualParameters(binary.fit, nsim=20)
summary(simpsi)
binary.fit@results@fixed.effects

# Simulate data
simdat<-simulate(binary.fit, nsim=2)
plot(binary.fit@data, type="binary")

simdat<-simulateDiscreteSaemix(binary.fit, nsim=2)
plotDiscreteData(binary.fit@data, outcome="binary")

# mirror plot
plotDiscreteDataElement(simdat, outcome="binary", mirror=TRUE)

################################################ Categorical data
# Knee data
knee.saemix<-read.table(file.path(datDir, "knee.saemix.tab"), header=TRUE)

knee.data<-saemixData(name.data=knee.saemix,name.group=c("id"),
                      name.predictors=c("y", "time"), name.X=c("time"),
                      name.covariates = c("Age","Sex","treatment","Age2"),
                      units=list(x="d",y="", covariates=c("yr","-","-","yr2")), verbose=FALSE)

# Fitting PO model 
ordinal.model<-function(psi,id,xidep) {
  y<-xidep[,1]
  time<-xidep[,2]
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  alp4<-psi[id,4]
  beta<-psi[id,5]
  
  logit1<-alp1 + beta*time
  logit2<-logit1+alp2
  logit3<-logit2+alp3
  logit4<-logit3+alp4
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  pge4<-exp(logit4)/(1+exp(logit4))
  pobs = (y==1)*pge1+(y==2)*(pge2 - pge1)+(y==3)*(pge3 - pge2)+(y==4)*(pge4 - pge3)+(y==5)*(1 - pge4)
  logpdf <- log(pobs)
  
  return(logpdf)
}
# simulate function
simulateOrdinal<-function(psi,id,xidep) {
  y<-xidep[,1]
  time<-xidep[,2]
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  alp4<-psi[id,4]
  beta<-psi[id,5]
  
  logit1<-alp1 + beta*time
  logit2<-logit1+alp2
  logit3<-logit2+alp3
  logit4<-logit3+alp4
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  pge4<-exp(logit4)/(1+exp(logit4))
  x<-runif(length(time))
  ysim<-1+as.integer(x>pge1)+as.integer(x>pge2)+as.integer(x>pge3)+as.integer(x>pge4)
  return(ysim)
}

covmodel2<-covmodel1<-matrix(data=0,ncol=5,nrow=4)
covmodel2[3,5]<-covmodel2[4,1]<-1
saemix.model.cov2<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood", 
                               simulate.function=simulateOrdinal,
                               psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                               transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)),
                               covariate.model = covmodel2)

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, nb.chains=10)
ord.fit.cov2<-saemix(saemix.model.cov2,knee.data,saemix.options)

# Simulate data
simord<-simulate(ord.fit.cov2, nsim=2)
plot(ord.fit.cov2@data, type="categorical") # not working
plotDiscreteData(simord@data, outcome="categorical")

# mirror plot
plotDiscreteData(simord@data, outcome="categorical")
plotDiscreteDataElement(simord, outcome="categorical")
plotDiscreteDataElement(simord, outcome="categorical", mirror=TRUE)
plotDiscreteDataElement(simord, outcome="categorical", mirror=TRUE, irep=2)

################################################ Count data
# RAPI data
rapi.saemix<-read.table(file.path(datDir, "rapi.saemix.tab"), header=TRUE)

# Data
rapi.data<-saemixData(name.data=rapi.saemix, name.group=c("id"),
                        name.predictors=c("time","rapi"),name.response=c("rapi"),
                        name.covariates=c("gender"),
                        units=list(x="months",y="",covariates=c("")))

## Models
# Poisson with a time effect
count.poisson<-function(psi,id,xidep) { 
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  logp <- -lambda + y*log(lambda) - log(factorial(y))
  return(logp)
}
saemix.simulatePoisson<-function(psi, id, xidep) {
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  y<-rpois(length(time), lambda=lambda)
  return(y)
}
# Fits
## Poisson
### Gender effect on intercept and slope
saemix.model.poi.cov2<-saemixModel(model=count.poisson,description="Count model Poisson",modeltype="likelihood",   
                                   simulate.function=saemix.simulatePoisson, 
                                   psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                   transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)),
                                   covariance.model =matrix(data=1, ncol=2, nrow=2),
                                   covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, fim=FALSE)

### Fit with saemix
poisson.fit.cov2<-saemix(saemix.model.poi.cov2,rapi.data,saemix.options)

simdat<-simulateDiscreteSaemix(poisson.fit.cov2, nsim=2)

plotDiscreteData(poisson.fit.cov2@data, outcome="count")

plotDiscreteDataElement(poisson.fit.cov2, outcome="count")
plotDiscreteDataElement(simdat, outcome="count", mirror=TRUE)

################################################ TTE 
# Lung cancer data
lung.saemix<-read.table(file.path(datDir, "lung.saemix.tab"), header=TRUE)

lung1<-lung.saemix
lung1$pat.karno[is.na(lung1$pat.karno)]<-median(lung1$pat.karno, na.rm=TRUE)

saemix.data.contPH<-saemixData(name.data=lung1,header=TRUE,name.group=c("id"),
                               name.predictors=c("time","status","cens"),name.response=c("status"),
                               name.covariates=c( "sex", "ph.ecog", "ph.karno", "pat.karno", "age"),
                               units=list(x="days",y="",covariates=c("","-","%","%","yr")), verbose=FALSE)

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
                          simulate.function=simulateWeibullTTE,
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                          transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)
tte.fit<-saemix(saemix.model,saemix.data.contPH,saemix.options)

covmodelcont <- cbind(c(1,1,0,0,0),rep(0,5))
weibull.model.cont<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                                simulate.function=simulateWeibullTTE,
                                psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                                transform.par=c(1,1),covariance.model=matrix(c(0,0,0,1),ncol=2, byrow=TRUE), 
                                covariate.model=covmodelcont, verbose=FALSE)
weibull.fit.cov<-saemix(weibull.model.cont,saemix.data.contPH,saemix.options)

plotDiscreteData(weibull.fit.cov@data, outcome="tte")

# Mirror plot
simdat<-simulateDiscreteSaemix(weibull.fit.cov, nsim=2)
plotDiscreteDataElement(weibull.fit.cov, outcome="tte")
plotDiscreteDataElement(simdat, outcome="tte", mirror=TRUE)

################################################ RTTE 
# TODO
