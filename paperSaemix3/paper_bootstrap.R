library(saemix)
saemixDir<-"/home/eco/work/saemix/saemixextension"
workDir<-file.path(saemixDir,"paperSaemix3")
nboot <- 5

# Sections to run (set to TRUE)

# binary data (toenail)
runmodelBin <- FALSE

# categorical data (knee)
runmodelCat <- TRUE
runmodelCatCov <- FALSE

# count data (rapi)
runmodelCount <- FALSE

# TTE model (lung)
runbasemodelTTE <- FALSE
runcovmodelContTTE<-FALSE
runcovmodelCatTTE<-FALSE

########################################################################################################
#################################################### binary data: toenail

########################################################################################################
#################################################### categorical data: knee

if(runmodelCat || runmodelCatCov) {
    # Data, restricted to the 2 covariates used in the model
    data(knee.saemix)
    ordknee.data<-saemixData(name.data=knee.saemix,name.group=c("id"),
                            name.predictors=c("y", "time"), name.X=c("time"),
                            name.covariates = c("treatment","Age2"),
                            units=list(x="d",y="", covariates=c("-","yr2")), verbose=FALSE)
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
}
if(runmodelCat) {
    # Model without covariate
    saemix.model<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                    simulate.function=simulateOrdinal, psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5, byrow=TRUE, 
                    dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), transform.par=c(0,1,1,1,1),
                    omega.init=diag(c(100, 1, 1, 1, 1)), covariance.model = diag(c(1,0,0,0,1)), verbose=FALSE)
    saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, nb.chains=10, nbiter.saemix=c(600,100), print=FALSE)
    ord.fit<-saemix(saemix.model,ordknee.data,saemix.options)            
    case.ordinal <- saemix.bootstrap(ord.fit, method="case", nboot=nboot) 
    write.table(case.ordinal,file.path(workDir, "bootstrapCase_knee.res"), row.names = FALSE, quote=F)
    cond.ordinal <- saemix.bootstrap(ord.fit, method="conditional", nboot=nboot) 
    write.table(cond.ordinal,file.path(workDir, "bootstrapCond_kneeCov.res"), row.names = FALSE, quote=F)
    }
if(runmodelCatCov) {
    # Model with covariate
    covariate.model <- matrix(data=0, nrow=2, ncol=5)
    covariate.model[1,2]<-covariate.model[1,5]<-covariate.model[2,1]<-1
    ordmodel.cov<-saemixModel(model=ordinal.model,description="Ordinal categorical model",
        modeltype="likelihood",simulate.function=simulateOrdinal, 
        psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5, byrow=TRUE, dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), transform.par=c(0,1,1,1,1),
        omega.init=diag(c(100, 1, 1, 1, 1)), covariate.model=covariate.model, 
        covariance.model = diag(c(1,1,1,1,0)), verbose=FALSE)
    ord.fit.cov<-saemix(ordmodel.cov,ordknee.data,saemix.options)
    
    # Bootstrap for both models 
    case.ordinal <- saemix.bootstrap(ord.fit.cov, method="case", nboot=nboot) 
    write.table(case.ordinal,file.path(workDir, "bootstrapCase_kneeCov.res"), row.names = FALSE, quote=F)
    cond.ordinal <- saemix.bootstrap(ord.fit.cov, method="conditional", nboot=nboot) 
    write.table(cond.ordinal,file.path(workDir, "bootstrapCond_kneeCov.res"), row.names = FALSE, quote=F)
}
########################################################################################################
#################################################### count data: rapi


########################################################################################################
#################################################### TTE analysis
# Data
data(lung.saemix)

# ECOG considered as continuous
# other covariates still have missing values
lung1<-lung.saemix
lung1$pat.karno[is.na(lung1$pat.karno)]<-median(lung1$pat.karno, na.rm=TRUE)

saemix.data.contPH<-saemixData(name.data=lung1,header=TRUE,name.group=c("id"),
                               name.predictors=c("time","status","cens"),name.response=c("status"),
                               name.covariates=c( "sex", "ph.ecog", "ph.karno", "pat.karno", "age"),
                               units=list(x="days",y="",covariates=c("","-","%","%","yr")), verbose=FALSE)

# ECOG considered as categorical
## Managing covariates - creating dummy covariates for ECOG=1 and ECOG=2 or 3, setting missing pat.karno to the median
lung2<-lung.saemix
lung2$ecog1<-ifelse(lung2$ph.ecog==1,1,0)
lung2$ecog23<-ifelse(lung2$ph.ecog>1,1,0)
lung2$pat.karno[is.na(lung2$pat.karno)]<-median(lung2$pat.karno, na.rm=TRUE)
saemix.data<-saemixData(name.data=lung2,header=TRUE,name.group=c("id"),
                        name.predictors=c("time","status","cens"),name.response=c("status"),
                        name.covariates=c("age", "sex", "ecog1","ecog23", "ph.karno", "pat.karno"),
                        units=list(x="days",y="",covariates=c("yr","","-","-","%","%")), verbose=FALSE)

####################################################
# Model
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

####################################################
# Without covariate
if(runbasemodelTTE) {
  weibull.model<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                             simulate.function = simulateWeibullTTE,
                             psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                             transform.par=c(1,1),covariance.model=matrix(c(0,0,0,1),ncol=2, byrow=TRUE), verbose=FALSE)
  weibtte.fit<-saemix(weibull.model,saemix.data,saemix.options)
  case.TTE <- saemix.bootstrap(weibtte.fit, method = "case")
  write.table(case.TTE,file.path(workDir, "bootstrapCase_weibullTTE.res"), row.names = FALSE, quote=F)
  cond.TTE <- saemix.bootstrap(weibtte.fit)
  write.table(cond.TTE,file.path(workDir, "bootstrapCond_weibullTTE.res"), row.names = FALSE, quote=F)
}

####################################################
# Covariate model, ECOG as continuous
if(runcovmodelContTTE) {
  covmodelcont <- cbind(c(1,1,0,0,0),rep(0,5))
  weibull.model.cont<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                                  simulate.function = simulateWeibullTTE,
                                  psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                                  transform.par=c(1,1),covariance.model=matrix(c(0,0,0,1),ncol=2, byrow=TRUE), 
                                  covariate.model=covmodelcont, verbose=FALSE)
  weibull.fit.cov<-saemix(weibull.model.cont,saemix.data.contPH,saemix.options)
  case.TTE <- saemix.bootstrap(weibull.fit.cov, method = "case")
  write.table(case.TTE,file.path(workDir, "bootstrapCase_weibullTTEcont.res"), row.names = FALSE, quote=F)
  cond.TTE <- saemix.bootstrap(weibull.fit.cov)
  write.table(cond.TTE,file.path(workDir, "bootstrapCond_weibullTTEcont.res"), row.names = FALSE, quote=F)
}

####################################################
# Covariate model, ECOG as categorical
if(runcovmodelCatTTE) {
  covmodel <- cbind(c(0,1,1,1,0,0),rep(0,6))
  weibull.model.cov2<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                                  simulate.function = simulateWeibullTTE,
                                  psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                                  transform.par=c(1,1),covariance.model=matrix(c(0,0,0,1),ncol=2, byrow=TRUE), 
                                  covariate.model=covmodel, verbose=FALSE)
  weibcov.fit2<-saemix(weibull.model.cov2,saemix.data,saemix.options)
  print(weibcov.fit2)
  case.TTE <- saemix.bootstrap(weibcov.fit2, method = "case")
  write.table(case.TTE,file.path(workDir, "bootstrapCase_weibullTTEcov.res"), row.names = FALSE, quote=F)
  cond.TTE <- saemix.bootstrap(weibcov.fit2, method = "conditional")
  write.table(cond.TTE,file.path(workDir, "bootstrapCond_weibullTTEcov.res"), row.names = FALSE, quote=F)
}
