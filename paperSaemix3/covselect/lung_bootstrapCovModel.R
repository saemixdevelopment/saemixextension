# Loading libraries
library(xtable)
library(ggplot2)

# Loading saemix
library(saemix)
source("/home/eco/work/saemix/saemixextension/R/func_bootstrap.R")

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

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)

covmodel <- cbind(c(0,1,0,1,0,0),rep(0,6))
weibull.model.cov2<-saemixModel(model=weibulltte.model,description="Weibull TTE model",modeltype="likelihood",
                                psi0=matrix(c(300,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("Te","gamma"))),
                                transform.par=c(1,1),covariance.model=matrix(c(0,0,0,1),ncol=2, byrow=TRUE), 
                                covariate.model=covmodel, verbose=FALSE)

###################################################### Fit

weibcov.fit2<-saemix(weibull.model.cov2,saemix.data,saemix.options)

###################################################### Bootstrap
nboot<-1000

case.TTE <- saemix.bootstrap(weibcov.fit2, method = "case", nboot=nboot)
write.table(case.TTE,file.path(workDir,"covselect", "bootstrapCase_weibullTTEcov.res"), row.names = FALSE, quote=F)

cond.TTE <- saemix.bootstrap(weibcov.fit2, nboot=nboot)
write.table(cond.TTE,file.path(workDir, "covselect", "bootstrapCond_weibullTTEcov.res"), row.names = FALSE, quote=F)

