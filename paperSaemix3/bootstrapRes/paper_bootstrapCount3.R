library(saemix)
saemixDir<-"/home/eco/work/saemix/saemixextension"
workDir<-file.path(saemixDir,"paperSaemix3","bootstrapRes")
nboot <- 200

# source(file.path(saemixDir,"R","func_bootstrap.R"))

data(rapi.saemix)
rapi.saemix$gender <- ifelse(rapi.saemix$gender=="Men",1,0) # Female=reference class as in Atkins

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
                            omega.init=diag(c(1,0.3)), verbose=FALSE)

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE, displayProgress=FALSE, print=FALSE)

hurdlefit0<-saemix(saemix.hurdle0,saemix.data0,saemix.options)
summary(hurdlefit0)

if(FALSE){
saemix.hurdle1.cov2<-saemixModel(model=truncated.poisson,description="Count model Poisson",modeltype="likelihood",   
                                 simulate.function = countsimulate.truncatedpoisson,
                                 psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                 transform.par=c(0,0), omega.init=diag(c(0.5, 0.01)),
                                 covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)

hurdlefit1<-saemix(saemix.hurdle1.cov2,saemix.data1,saemix.options)
summary(hurdlefit1)
  case.count <- saemix.bootstrap(hurdlefit1, method="case", nboot=nboot) 
  write.table(case.count,file.path(workDir, "bootstrapCase_rapiCov_Hurdle1.res"), row.names = FALSE, quote=F)
cond.count <- saemix.bootstrap(hurdlefit1, method="conditional", nboot=nboot) 
write.table(cond.count,file.path(workDir, "bootstrapCond_rapiCov_Hurdle1.res"), row.names = FALSE, quote=F)
}
case.count <- saemix.bootstrap(hurdlefit0, method="case", nboot=nboot) 
write.table(case.count,file.path(workDir, "bootstrapCase_rapiCov_Hurdle0.res"), row.names = FALSE, quote=F)
cond.count <- saemix.bootstrap(hurdlefit0, method="conditional", nboot=nboot) 
write.table(cond.count,file.path(workDir, "bootstrapCond_rapiCov_Hurdle0.res"), row.names = FALSE, quote=F)

########################################################################################################
