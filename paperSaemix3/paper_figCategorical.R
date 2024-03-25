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

###################################################### Data exploration
# Data
data(knee.saemix)

ordknee.data<-saemixData(name.data=knee.saemix,name.group=c("id"),
                        name.predictors=c("y", "time"), name.X=c("time"),
                        name.covariates = c("Age","Sex","treatment","Age2"),
                        units=list(x="d",y="", covariates=c("yr","-","-","yr2")), verbose=FALSE)

plotDiscreteData(ordknee.data, outcome="categorical", which.cov="treatment")

plotDiscreteData(ordknee.data, outcome="categorical", which.cov="Sex")

# ERROR TODO: should be time (d) on the X-axis...not y

if(saveFigs) {
  namfig<-"knee_rawDataPropTime.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plotDiscreteData(ordknee.data, outcome="categorical", which.cov="treatment")
  dev.off()
}

###################################################### Fit proportional odds model

# Model for ordinal responses
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

# Fitting
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, nb.chains=10, nbiter.saemix=c(600,100), print=FALSE)
#saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, nb.chains=10, fim=FALSE)

# Saemix model without covariates
saemix.model<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                            simulate.function=simulateOrdinal, psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5, byrow=TRUE, 
                                                                           dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), transform.par=c(0,1,1,1,1),
                            omega.init=diag(c(100, 1, 1, 1, 1)), covariance.model = diag(c(1,0,0,0,1)), verbose=FALSE)
ord.fit<-saemix(saemix.model,ordknee.data,saemix.options)
summary(ord.fit)

# Model with covariates resulting from the stepwise algorithm
covariate.model <- matrix(data=0, nrow=2, ncol=5)
covariate.model[1,2]<-covariate.model[1,5]<-covariate.model[2,1]<-1
ordmodel.cov<-saemixModel(model=ordinal.model,description="Ordinal categorical model",
                          modeltype="likelihood",simulate.function=simulateOrdinal, 
                          psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5, byrow=TRUE, dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), transform.par=c(0,1,1,1,1),
                          omega.init=diag(c(100, 1, 1, 1, 1)), covariate.model=covariate.model, 
                          covariance.model = diag(c(1,1,1,1,0)), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, nb.chains=10, nbiter.saemix=c(600,100), print=FALSE)
ord.covfit<-saemix(ordmodel.cov,ordknee.data,saemix.options)            

###################################################### Covariate models

###################################################### SE - redo for covariate model
if(FALSE) {
  # Redone after changes to the bootstrap code (namely in starting values)
  case.ordinal2 <- read.table(file.path(saemixDir,"paperSaemix3","bootstrapCase_knee.res"), header=T)
  cond.ordinal2 <- read.table(file.path(saemixDir,"paperSaemix3","bootstrapCond_knee.res"), header=T)
  
  # previous implementation
  case.ordinal3 <- read.table(file.path(saemixDir,"bootstrap","results","knee_caseBootstrap.res"), header=T)
  cond.ordinal3 <- read.table(file.path(saemixDir,"bootstrap","results","knee_condBootstrap.res"), header=T)
  
  # comparing both fits
  summary(case.ordinal2)
  summary(case.ordinal3)
  # comparing with conditional bootstrap
  summary(cond.ordinal2)
}

case.ordinal <- read.table(file.path(saemixDir,"paperSaemix3","bootstrapCase_knee.res"), header=T)
cond.ordinal <- read.table(file.path(saemixDir,"paperSaemix3","bootstrapCond_knee.res"), header=T)
nboot<-dim(case.ordinal)[1]
case.ordinal <- case.ordinal[!is.na(case.ordinal[,2]),]
cond.ordinal <- cond.ordinal[!is.na(cond.ordinal[,2]),]

par.estim<-format(c(ord.fit@results@fixed.effects,diag(ord.fit@results@omega)[ord.fit@results@indx.omega]), digits=2, nsmall=1)
df2<-data.frame(parameter=colnames(case.ordinal)[-c(1)], saemix=par.estim)
for(i in 1:2) {
  if(i==1) {
    resboot1<-case.ordinal
    namboot<-"case"
  } else {
    resboot1<-cond.ordinal
    namboot <-"cNP"
  }
  mean.bootDist<-apply(resboot1, 2, mean, na.rm=T)[-c(1)]
  sd.bootDist<-apply(resboot1, 2, sd, na.rm=T)[-c(1)]
  quant.bootDist<-apply(resboot1[-c(1)], 2, quantile, c(0.025, 0.975), na.rm=T)
  l1<-paste0(format(mean.bootDist, digits=1, scientific=FALSE)," (",format(sd.bootDist,digits=1, nsmall=0, trim=T, scientific=FALSE),")")
  l2<-paste0("[",format(quant.bootDist[1,], digits=1, trim=T, scientific=FALSE),", ",format(quant.bootDist[2,],digits=1, trim=T, scientific=FALSE),"]")
  df2<-cbind(df2, l1, l2)
  i1<-3+2*(i-1)
  colnames(df2)[i1:(i1+1)]<-paste0(namboot,".",c("estimate","CI"))
}
print(df2)
df2.nocov<-df2


case.ordinal <- read.table(file.path(saemixDir,"paperSaemix3","bootstrapCase_kneeCov.res"), header=T)
cond.ordinal <- read.table(file.path(saemixDir,"paperSaemix3","bootstrapCond_kneeCov.res"), header=T)
nboot<-dim(case.ordinal)[1]
case.ordinal <- case.ordinal[!is.na(case.ordinal[,2]),]
cond.ordinal <- cond.ordinal[!is.na(cond.ordinal[,2]),]

par.estim<-format(c(ord.covfit@results@fixed.effects,diag(ord.covfit@results@omega)[ord.covfit@results@indx.omega]), digits=2, nsmall=1)
df2<-data.frame(parameter=colnames(case.ordinal)[-c(1)], saemix=par.estim)
for(i in 1:2) {
  if(i==1) {
    resboot1<-case.ordinal
    namboot<-"case"
  } else {
    resboot1<-cond.ordinal
    namboot <-"cNP"
  }
  mean.bootDist<-apply(resboot1, 2, mean, na.rm=T)[-c(1)]
  sd.bootDist<-apply(resboot1, 2, sd, na.rm=T)[-c(1)]
  quant.bootDist<-apply(resboot1[-c(1)], 2, quantile, c(0.025, 0.975), na.rm=T)
  l1<-paste0(format(mean.bootDist, digits=1, scientific=FALSE)," (",format(sd.bootDist,digits=1, nsmall=0, trim=T, scientific=FALSE),")")
  l2<-paste0("[",format(quant.bootDist[1,], digits=1, trim=T, scientific=FALSE),", ",format(quant.bootDist[2,],digits=1, trim=T, scientific=FALSE),"]")
  df2<-cbind(df2, l1, l2)
  i1<-3+2*(i-1)
  colnames(df2)[i1:(i1+1)]<-paste0(namboot,".",c("estimate","CI"))
}
print(df2)
df2.cov<-df2


df3<-df2.cov[,c(2,6)] # conditional bootstrap
df3<-rbind(df3, c("-", "-"))
rownames(df3)<-c("$\\alpha_1$","$\\beta_{Age^2,\\alpha_1}$","$\\alpha_2$","$\\beta_{Trt,\\alpha_2}$","$\\alpha_3$","$\\alpha_4$","$\\beta$", "$\\beta_{Trt,\\beta}$", "$\\omega_{\\alpha_1}$","$\\omega_{\\alpha_2}$", "$\\omega_{\\alpha_3}$", "$\\omega_{\\alpha_4}$", "$\\omega_{\\beta}$")
df4<-df2.nocov[,c(2,6)] # model without covariates
df3.nocov<-rbind(df4[1,], c("-", "-"), df4[2,],c("-", "-"),df4[3:5,],c("-", "-"),df4[6,],c("-", "-"),c("-", "-"),c("-", "-"),df4[7,])
df3<-cbind(df3.nocov, df3)
rownames(df3)<-c("$\\alpha_1$","$\\beta_{Age^2,\\alpha_1}$","$\\alpha_2$","$\\beta_{Trt,\\alpha_2}$","$\\alpha_3$","$\\alpha_4$","$\\beta$", "$\\beta_{Trt,\\beta}$", "$\\omega_{\\alpha_1}$","$\\omega_{\\alpha_2}$", "$\\omega_{\\alpha_3}$", "$\\omega_{\\alpha_4}$", "$\\omega_{\\beta}$")

print(xtable(df3), only.contents=TRUE, include.rownames=TRUE,  floating=F, sanitize.rownames.function = identity)

