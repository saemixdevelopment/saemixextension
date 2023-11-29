# Loading libraries
library(xtable)
library(ggplot2)
library(tidyr)

# Loading saemix
#library(saemix)

# Folders
workDir<-"/home/eco/work/saemix/discreteEval"
# setwd(workDir)

saemixDir<-"/home/eco/work/saemix/saemixextension"
figDir <- file.path(workDir, "figs")
nsim<-200 # Number of simulations
saveFigures <- FALSE

###################################################### Data exploration

data(rapi.saemix) 
rapi.saemix$gender <- ifelse(rapi.saemix$gender=="Men",1,0)  # Female=reference class as in Atkins

saemix.data<-saemixData(name.data=rapi.saemix, name.group=c("id"),
                        name.predictors=c("time","rapi"),name.response=c("rapi"),
                        name.covariates=c("gender"),
                        units=list(x="months",y="",covariates=c(""), verbose=FALSE)

xpl <- plotDiscreteData(saemix.data, outcome="count", which.cov="gender", breaks=c(0:9, 16, 25,80))
print(xpl)

###################################################### Fit Poisson model
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
  y<-xidep[,2] # used to set an upper bound to counts
  ymax<-max(y)
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  y<-rpois(length(time), lambda=lambda)
  y[y>ymax]<-ymax+1 # truncate to maximum observed value to avoid simulating aberrant values
  return(y)
}

### Gender effect on intercept and slope
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

# Proportion of zeroes
nsim<-100
yfit1<-simulateDiscreteSaemix(poisson.fit.cov2, nsim=nsim)
obs.prop0 <-  length(yfit1@data@data$rapi[yfit1@data@data$rapi==0])/yfit1@data@ntot.obs
poiss.prop0 <-length(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim==0])/length(yfit1@sim.data@datasim$ysim)
cat("Observed proportion of 0's",obs.prop0,"\n")
cat("      Poisson model, p=",poiss.prop0,"\n")

###################################################### Fit ZI-Poisson model
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

### ZIP Poisson with gender on both intercept and slope
saemix.model.zip.cov2<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                                   simulate.function = countsimulate.poissonzip,
                                   psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                                   transform.par=c(0,0,3), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.3,0)),
                                   covariate.model = matrix(c(1,1,0),ncol=3, byrow=TRUE), verbose=FALSE)

zippoisson.fit.cov2<-saemix(saemix.model.zip.cov2,saemix.data,saemix.options)

# Proportion of zeroes
ysim.zip<-simulateDiscreteSaemix(zippoisson.fit.cov2, nsim=nsim)
zip.prop0 <-length(ysim.zip@sim.data@datasim$ysim[ysim.zip@sim.data@datasim$ysim==0])/length(ysim.zip@sim.data@datasim$ysim)
cat("      ZI-Poisson model, p=",zip.prop0,"\n")

###################################################### Fit hurdle model from Atkins
# Separate in two datasets
saemix.data1<-saemixData(name.data=rapi.saemix[rapi.saemix$rapi>0,], name.group=c("id"),
                         name.predictors=c("time","rapi"),name.response=c("rapi"),
                         name.covariates=c("gender"),
                         units=list(x="week",y="",covariates=c("")))

rapi.saemix$y0<-as.integer(rapi.saemix$rapi>0)
saemix.data0<-saemixData(name.data=rapi.saemix, name.group=c("id"),
                         name.predictors=c("time","y0"),name.response=c("y0"),
                         name.covariates=c("gender"),
                         units=list(x="week",y="",covariates=c("")))

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
saemix.hurdle0<-saemixModel(model=binary.model,description="Binary model",
                            modeltype="likelihood",simulate.function=simulBinary,
                            psi0=matrix(c(0.5,-.1,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                            transform.par=c(0,0), covariate.model=c(1,1),
                            covariance.model=matrix(c(1,0,0,1),ncol=2), omega.init=diag(c(1,0.3)))

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE, displayProgress=FALSE)

hurdlefit0<-saemix(saemix.hurdle0,saemix.data0,saemix.options)

# Fit Poisson model to saemix.data1
saemix.hurdle1.cov2<-saemixModel(model=count.poisson,description="Count model Poisson",modeltype="likelihood",   
                                 simulate.function = countsimulate.poisson,
                                 psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                 transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)),
                                 covariance.model =matrix(data=1, ncol=2, nrow=2),
                                 covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE))
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)

hurdlefit1<-saemix(saemix.hurdle1.cov2,saemix.data1,saemix.options)

summary(hurdlefit0)
summary(hurdlefit1)

# VPC
yhurdle0<-simulateDiscreteSaemix(hurdlefit0, nsim=nsim)
discreteVPC(yhurdle0)
hurd.prop0 <-length(yhurdle0@sim.data@datasim$ysim[yhurdle0@sim.data@datasim$ysim==0])/length(yhurdle0@sim.data@datasim$ysim)

###################################################### Results in a table
# Comparing models
xtab1 <- cbind("Poisson"=c(format((poisson.fit.cov2@results@fixed.effects), digits=1),"",format(BIC(poisson.fit.cov2), digits=1, nsmall=1),format(poiss.prop0, digits=2)),
      "ZIP"=c(format((zippoisson.fit.cov2@results@fixed.effects)[1:4], digits=1),format(zippoisson.fit.cov2@results@fixed.effects[5], digits=1),format(BIC(zippoisson.fit.cov2), digits=1, nsmall=1),format(zip.prop0, digits=2)),
      "Hurdle"=c(format((hurdlefit1@results@fixed.effects), digits=1),"","-",format(hurd.prop0, digits=2))
)
rownames(xtab1)<-c("$\\alpha_0$","$\\beta_{male, int}$","$\\alpha_1$ (month$^{-1}$)","$\\beta_{male,time}$","$p_0$","BIC","Prop(Y=0)")

print(xtable(xtab1), only.contents=TRUE, include.rownames=T, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity)

# Bootstrap SE
case.count <- read.table(file.path(saemixDir,"bootstrap","results","rapi_caseBootstrap.res"), header=T)
cond.count <- read.table(file.path(saemixDir,"bootstrap","results","rapi_condBootstrap.res"), header=T)
nboot<-dim(case.count)[1]
case.count <- case.count[!is.na(case.count[,2]),]
cond.count <- cond.count[!is.na(cond.count[,2]),]

par.estim<-c(zippoisson.fit.cov2@results@fixed.effects,sqrt(diag(zippoisson.fit.cov2@results@omega)[zippoisson.fit.cov2@results@indx.omega]))
df2<-data.frame(parameter=colnames(case.count)[-c(1)], saemix=format(par.estim, digits = 1, nsmall=2))
for(i in 1:2) {
  if(i==1) {
    resboot1<-case.count
    namboot<-"case"
  } else {
    resboot1<-cond.count
    namboot <-"cNP"
  }
  for(icol in 7:8) resboot1[,icol]<-sqrt(resboot1[,icol])
  mean.bootDist<-apply(resboot1, 2, mean)[-c(1)]
  sd.bootDist<-apply(resboot1, 2, sd)[-c(1)]
  quant.bootDist<-apply(resboot1[-c(1)], 2, quantile, c(0.025, 0.975))
  l1<-paste0(format(mean.bootDist, digits=1, nsmall=2)," (",format(sd.bootDist,digits=1, nsmall=2, trim=T),")")
  l2<-paste0("[",format(quant.bootDist[1,], digits=1, nsmall=2),", ",format(quant.bootDist[2,],digits=1, nsmall=2, trim=T),"]")
  df2<-cbind(df2, l1, l2)
  i1<-3+2*(i-1)
  colnames(df2)[i1:(i1+1)]<-paste0(namboot,".",c("estimate","CI"))
}
print(df2)

rownames(df2)<-c("$\\alpha_0$","$\\beta_{male, int}$","$\\alpha_1$  (month$^{-1}$)","$\\beta_{male,time}$","$p_0$","$\\omega_0$","$\\omega_1$")
print(xtable(df2[,2:4]), only.contents=TRUE, include.rownames=T, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity)


###################################################### Diagnostics

discreteVPC(ysim.zip, outcome="count", breaks=c(0:9,16,25,80), which.cov="gender")

###################################################### Quick check
if(FALSE) {
  tab.men <- table(rapi.saemix$rapi[rapi.saemix$time==0 & rapi.saemix$gender==1] == 0)
  tab.women <- table(rapi.saemix$rapi[rapi.saemix$time==0 & rapi.saemix$gender==0] == 0)
  cat("Observed proportion of 0's at time 0 (in women):",tab.women[2]/sum(tab.women),"\n")
  
  xcal<-rnorm(1000, mean=hurdlefit0@results@fixed.effects[1], sd=sqrt(hurdlefit0@results@omega[1,1]))
  cat("Expected proportion of 0's at time 0 (in women):",mean(1-1/(1+exp(-xcal))),"\n")

xcal1<-rnorm(1000, mean=1.56, sd=sqrt(2.30766)) # omega not reported in Atkins et al.
cat("Expected proportion of 0's at time 0 (in women) from Atkins, assuming same variance:",mean(1-1/(1+exp(-xcal1))),"\n")
xcal2<-rnorm(1000, mean=1.56, sd=sqrt(0.01)) # omega not reported in Atkins et al.
mean(1-1/(1+exp(-xcal2)))
}
