library(saemix)
saemixDir<-"/home/eco/work/saemix/saemixextension"
workDir<-file.path(saemixDir,"paperSaemix3","bootstrapRes")
nboot <- 200

# source(file.path(saemixDir,"R","func_bootstrap.R"))

data(rapi.saemix)
rapi.saemix$gender <- ifelse(rapi.saemix$gender=="Men",1,0) # Female=reference class as in Atkins

########################################################################################################
# Poisson model

saemix.data<-saemixData(name.data=rapi.saemix, name.group=c("id"),
                        name.predictors=c("time","rapi"),name.response=c("rapi"),
                        name.covariates=c("gender"),
                        units=list(x="months",y="",covariates=c("")), verbose=FALSE)

## Poisson with a time effect
# Model
count.poisson<-function(psi,id,xidep) { 
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  logp <- -lambda + y*log(lambda) - log(factorial(y))
  return(logp)
}
# Simulation function
countsimulate.poisson<-function(psi, id, xidep) {
  time<-xidep[,1]
  y<-xidep[,2]
  ymax<-max(y)
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  y<-rpois(length(time), lambda=lambda)
  y[y>ymax]<-ymax+1 # truncate to maximum observed value to avoid simulating aberrant values
  return(y)
}

saemix.model.poi.cov2<-saemixModel(model=count.poisson,description="Count model Poisson",simulate.function=countsimulate.poisson, 
                                   modeltype="likelihood",   
                                   psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                   transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)),
                                   covariance.model =matrix(data=1, ncol=2, nrow=2),
                                   covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE), verbose=FALSE)

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, fim=FALSE, print=FALSE)

### Fit with saemix
poisson.fit.cov2<-saemix(saemix.model.poi.cov2,saemix.data,saemix.options)
summary(poisson.fit.cov2)

case.count <- saemix.bootstrap(poisson.fit.cov2, method="case", nboot=nboot) 
write.table(case.count,file.path(workDir, "bootstrapCase_rapiCov_Poisson.res"), row.names = FALSE, quote=F)

########################################################################################################
# ZIP Poisson model

## Zero-inflated Poisson model
# Model
count.poissonzip<-function(psi,id,xidep) {
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  p0<-psi[id,3] # Probability of zero's
  lambda<- exp(intercept + slope*time)
  logp <- log(1-p0) -lambda + y*log(lambda) - log(factorial(y)) # Poisson
  logp0 <- log(p0+(1-p0)*exp(-lambda)) # Zeroes
  logp[y==0]<-logp0[y==0]
  return(logp)
}
# Simulation function
countsimulate.poissonzip<-function(psi, id, xidep) {
  time<-xidep[,1]
  y<-xidep[,2]
  ymax<-max(y)
  intercept<-psi[id,1]
  slope<-psi[id,2]
  p0<-psi[id,3] # Probability of zero's
  lambda<- exp(intercept + slope*time)
  prob0<-rbinom(length(time), size=1, prob=p0)
  y<-rpois(length(time), lambda=lambda)
  y[prob0==1]<-0
  y[y>ymax]<-ymax+1 # truncate to maximum observed value to avoid simulating aberrant values
  return(y)
}

saemix.model.zip.cov2<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                                   simulate.function = countsimulate.poissonzip,
                                   psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                                   transform.par=c(0,0,3), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.3,0)),
                                   covariate.model = matrix(c(1,1,0),ncol=3, byrow=TRUE), verbose=FALSE)

zippoisson.fit.cov2<-saemix(saemix.model.zip.cov2,saemix.data,saemix.options)
summary(zippoisson.fit.cov2)

case.count <- saemix.bootstrap(zippoisson.fit.cov2, method="case", nboot=nboot) 
write.table(case.count,file.path(workDir, "bootstrapCase_rapiCov_ZIP.res"), row.names = FALSE, quote=F)

# with correlation
mat1<-diag(c(1,1,0))
mat1[1,2]<-mat1[2,1]<-1
saemix.model.zip.corr<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                                   simulate.function = countsimulate.poissonzip,
                                   psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                                   transform.par=c(0,0,3), covariance.model=mat1, omega.init=diag(c(0.5,0.3,0)),
                                   covariate.model = matrix(c(1,1,0),ncol=3, byrow=TRUE), verbose=FALSE)

zippoisson.fit.corr<-saemix(saemix.model.zip.corr,saemix.data,saemix.options)
summary(zippoisson.fit.corr)

case.count <- saemix.bootstrap(zippoisson.fit.corr, method="case", nboot=nboot) 
write.table(case.count,file.path(workDir, "bootstrapCase_rapiCov_ZIP_corr.res"), row.names = FALSE, quote=F)

cond.count <- saemix.bootstrap(zippoisson.fit.corr, method="conditional", nboot=nboot) 


########################################################################################################
# Hurdle model - without correlation

## Hurdle - 2 models 
saemix.data1<-saemixData(name.data=rapi.saemix[rapi.saemix$rapi>0,], name.group=c("id"),
                         name.predictors=c("time","rapi"),name.response=c("rapi"),
                         name.covariates=c("gender"),
                         units=list(x="week",y="",covariates=c("")), verbose=FALSE)

rapi.saemix$y0<-as.integer(rapi.saemix$rapi>0)
saemix.data0<-saemixData(name.data=rapi.saemix, name.group=c("id"),
                         name.predictors=c("time","y0"),name.response=c("y0"),
                         name.covariates=c("gender"),
                         units=list(x="week",y="",covariates=c("")), verbose=FALSE)

# Fit Binomial model to saemix.data0
binary.model<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  pobs = (y==0)*(1-pevent)+(y==1)*pevent
  logpdf <- log(pobs)
  return(logpdf)
}
# Associated simulation function
simulBinary<-function(psi,id,xidep) {
  tim<-xidep[,1]
  #  y<-xidep[,2] # not used for simulation
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  ysim<-rbinom(length(tim),size=1, prob=pevent)
  return(ysim)
}

# Fit truncated Poisson model to saemix.data1
## Poisson with a time effect
# Model
truncated.poisson<-function(psi,id,xidep) { 
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  logp <- y*log(lambda) - log(factorial(y)) - log(exp(lambda)-1)
  return(logp)
}
# simulating from a truncated Poisson using the inverse function
# Peter Dalgaard on R-help list: https://stat.ethz.ch/pipermail/r-help/2005-May/070680.html
rtpois <- function(N, lambda) {
  qpois(runif(N, dpois(0, lambda), 1), lambda)
}

# Simulation function
countsimulate.truncatedpoisson<-function(psi, id, xidep) {
  # left: truncate to 0
  # right: truncate to maximum observed value to avoid simulating aberrant values
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  #  y<-rtruncPois(length(time), y, lambda=lambda) 
  y<-rtpois(length(time), lambda=lambda)
  return(y)
}

########################## No correlation between parameters

saemix.hurdle0<-saemixModel(model=binary.model,description="Binary model",
                            modeltype="likelihood",simulate.function=simulBinary,
                            psi0=matrix(c(0.5,-.1,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                            transform.par=c(0,0), covariate.model=c(1,1),
                            covariance.model =matrix(data=1, ncol=2, nrow=2), omega.init=diag(c(1,0.3)), verbose=FALSE)

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE, displayProgress=FALSE, print=FALSE)

hurdlefit0<-saemix(saemix.hurdle0,saemix.data0,saemix.options)
summary(hurdlefit0)

saemix.hurdle1.cov2<-saemixModel(model=truncated.poisson,description="Count model Poisson",modeltype="likelihood",   
                                 simulate.function = countsimulate.truncatedpoisson,
                                 psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                 transform.par=c(0,0), omega.init=diag(c(0.5, 0.01)),
                                 covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)

hurdlefit1<-saemix(saemix.hurdle1.cov2,saemix.data1,saemix.options)
summary(hurdlefit1)

case.count <- saemix.bootstrap(hurdlefit0, method="case", nboot=nboot) 
write.table(case.count,file.path(workDir, "bootstrapCase_rapiCov_Hurdle0.res"), row.names = FALSE, quote=F)

case.count <- saemix.bootstrap(hurdlefit1, method="case", nboot=nboot) 
write.table(case.count,file.path(workDir, "bootstrapCase_rapiCov_Hurdle1.res"), row.names = FALSE, quote=F)

########################## With correlation between parameters

saemix.hurdle0.corr<-saemixModel(model=binary.model,description="Binary model",
                            modeltype="likelihood",simulate.function=simulBinary,
                            psi0=matrix(c(0.5,-.1,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                            covariance.model =matrix(data=1, ncol=2, nrow=2),
                            transform.par=c(0,0), covariate.model=c(1,1),
                            omega.init=diag(c(1,0.3)), verbose=FALSE)

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE, displayProgress=FALSE, print=FALSE)

hurdlefit0.corr<-saemix(saemix.hurdle0.corr,saemix.data0,saemix.options)
summary(hurdlefit0.corr)

saemix.hurdle1.corr<-saemixModel(model=truncated.poisson,description="Count model Poisson",modeltype="likelihood",   
                                 simulate.function = countsimulate.truncatedpoisson,
                                 psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                 transform.par=c(0,0), omega.init=diag(c(0.5, 0.01)),
                                 covariance.model =matrix(data=1, ncol=2, nrow=2),
                                 covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)

hurdlefit1.corr<-saemix(saemix.hurdle1.corr,saemix.data1,saemix.options)
summary(hurdlefit1.corr)

case.count <- saemix.bootstrap(hurdlefit0.corr, method="case", nboot=nboot) 
write.table(case.count,file.path(workDir, "bootstrapCase_rapiCov_Hurdle0_corr.res"), row.names = FALSE, quote=F)

case.count <- saemix.bootstrap(hurdlefit1.corr, method="case", nboot=nboot) 
write.table(case.count,file.path(workDir, "bootstrapCase_rapiCov_Hurdle1_corr.res"), row.names = FALSE, quote=F)

########################################################################################################
