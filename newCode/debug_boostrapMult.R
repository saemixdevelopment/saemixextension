# Debugging bootstrap for multiple response models
####################################################################
# Loading Lucas' functions

saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-"/home/eco/work/saemix/versions/lucas/prov/saemixDev/R"
datDir<-file.path(saemixDir,"data")
datDir2 <- "/home/eco/work/saemix/workshop26/data"

# Libraries
library(ggplot2)
library(MASS)
library(rlang)
library(npde)

source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixData-methods.R"))
source(file.path(progDir,"SaemixData-methods_covariates.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"main.R"))
source(file.path(progDir,"func_aux.R"))
source(file.path(progDir,"main_initializeMainAlgo.R")) # ...
source(file.path(progDir,"main_estep.R"))
source(file.path(progDir,"main_mstep.R"))
source(file.path(progDir,"func_FIM.R"))
source(file.path(progDir,"func_map.R"))
source(file.path(progDir,"func_plots.R"))
source(file.path(progDir,"func_distcond.R"))
source(file.path(progDir,"func_simulations.R"))
source(file.path(progDir,"compute_LL.R"))
source(file.path(progDir,"func_npde.R"))
source(file.path(progDir,"func_estimParam.R"))
source(file.path(progDir,"backward.R"))
source(file.path(progDir,"forward.R"))
source(file.path(progDir,"func_stepwise.R"))
source(file.path(progDir,"stepwise.R"))
source(file.path(progDir,"func_compare.R"))
#source(file.path(progDir,"func_bootstrap.R"))
source(file.path(progDir,"func_exploreData.R"))
source(file.path(progDir,"func_discreteVPC.R"))
source(file.path(progDir,"progressBar.R"))

# Bootstrap function, modified
source(file.path(saemixDir,"R","func_bootstrap_mult.R"))

####################################################################
# Joint model
prothro.saemix <- read.table(file.path(datDir2,"prothro.saemix.tab"), header=TRUE)

joint.data<-saemixData(name.data=prothro.saemix, name.group=c("id"), 
                       name.predictors=c("time","obs"),
                       name.response="obs",name.ytype = "ytype")

# Joint model function
JMmodel<-function(psi,id,xidep) {
  time <- xidep$time
  ytype<-xidep$ytype  # type of response (1: continuous, 2: event)
  b0 <- psi[id,1] 
  b1 <- psi[id,2] 
  h0 <- psi[id,3]
  alpha <- psi[id,4] 
  # Longitudinal response
  ypred <- b0+b1*time # predictions for the longitudinal part
  # Event outcome
  T<-time[ytype==2] # vector of times (survival response)
  Nj <- length(T)
  ev <- xidep$obs[ytype==2] # vector of observations (survival response)
  cens <- which(ev==0)       # with censored ones
  ind <- which(ev==1)      # and event ones 
  # Creating vectors of the same length of T to compute likelihood of the survival part 
  b0b <- b0[ytype==2] # to have vectors of the same length as T 
  b1b <- b1[ytype==2]
  h0b <- h0[ytype==2]
  alphab <- alpha[ytype==2]
  
  haz <- h0b*exp(alphab*(b0b+b1b*T)) # instantaneous hazard
  # cumulative hazard (explicit expression in our case)
  H <- (h0b/(alphab*b1b))*exp((b0b+b1b*T)*alphab)-(h0b/(alphab*b1b))*exp(alphab*b0b) 
  
  logpdf <- rep(0,Nj)
  logpdf[cens] <- -H[cens] # likelihood contributions for censored observations
  logpdf[ind] <- -H[ind] + log(haz[ind]) # likelihood contributions for event observations 
  
  ypred[ytype==2] = logpdf
  return(ypred)
}

# Initial parameters
param<-c(73,1.25,0.6,0.0001) 
omega.sim<-c(18, 3, 0.05, 0.01)
sigma.sim <- 17

### saemix Model 
joint.model<-saemixModel(model=JMmodel,description="JM LMEM-TTE constant baseline hazard", 
                         modeltype=c("structural","likelihood"),
                         psi0=matrix(param,ncol=4,byrow=TRUE,
                                     dimnames=list(NULL, c("b0","b1","h0","alpha"))),
                         transform.par=c(0,0,1,0), covariance.model=diag(c(1,1,0,0)),
                         fixed.estim = c(1,1,1,1),error.model = "constant",
                         omega.init = diag(omega.sim))

saemix.options<-saemixControl(seed=12345, map=T, fim=T, ll.is=TRUE, save.graphs = F) 
# please, specify save.graphs=F (currently not extended)
fit.joint <- saemix(joint.model, joint.data, saemix.options)

###################################
# Bootstrap for joint model
nboot <- 2
case.count <- saemix.bootstrap(fit.joint, method="case", nboot=nboot) 
case.count

# Debugging Case bootstrap for joint models
## only case for the moment (difficult to simulate from a joint model)
####################################################################
# Multiple response model - 2 Gaussian responses
load(file.path(datDir2,"sd_iv_rich_pkpd.rda"))
tab1 <- sd_iv_rich_pkpd[,-c(4,10)]
tab2 <- sd_iv_rich_pkpd[,-c(3,10)]
colnames(tab1)<-colnames(tab2)<-c("Id","Time","Obs","Weight","Age","Dose","Sex","Race")
tab1$ytype <- 1
tab2$ytype <- 2
datPKPD <- rbind(tab1,tab2)
datPKPD <- datPKPD[order(datPKPD$Id, datPKPD$Time, datPKPD$ytype),]
pkpd.data<-saemixData(name.data=datPKPD,
                      name.group=c("Id"),name.predictors=c("Time","Dose"),name.ytype="ytype",
                      name.response=c("Obs"),name.covariates=c("Weight","Age","Sex"),
                      units=list(x="hr",y=c("mg/L","-"),covariates=c("kg","yr","-")))
# Centering covariates
pkpd.data<-transformContCov(pkpd.data,Weight,centering=70,transformation=log,verbose=FALSE, newCovName = "lWeight")
pkpd.data<-transformContCov(pkpd.data,Age,centering="median",transformation=log,verbose=FALSE, newCovName = "lAge")
# Recreating pkpd.data with our new covariates
pkpd.data<-saemixData(name.data=pkpd.data@data,
                      name.group=c("Id"),name.predictors=c("Time","Dose"),name.ytype="ytype",
                      name.response=c("Obs"),name.covariates=c("lWeight","lAge","Sex"),
                      units=list(x="hr",y=c("mg/L","-"),covariates=c("logkg","logyr","-")))

datPKPD.transcov <- pkpd.data@data[,c("Id","Time","Obs","lWeight","lAge","Dose","Sex","ytype")]

pkpdmodel <- function(psi,id,xidep) { 
  tim<-xidep[,1]  
  dose<-xidep[,2]
  ytype<-xidep$ytype
  CL<-psi[id,1]
  V<-psi[id,2]
  Emax<-psi[id,3]
  EC50<-psi[id,4]
  k<-CL/V
  ypred<-dose/V*(exp(-k*tim))
  ypred2 <- Emax*ypred/(ypred+EC50)
  ypred[ytype==2] <- ypred2[ytype==2]
  return(ypred)
}
psi0 <- c(0.9, 9, 9.5, 1.2)
pkpd.model<-saemixModel(model=pkpdmodel,modeltype=c("structural","structural"),
                        description="PK/PD model with direct response model", 
                        psi0=matrix(c(psi0,0.75,1,0,0),ncol=4, byrow=TRUE,
                                    dimnames=list(NULL, c("CL","V","Emax","EC50"))),transform.par=c(1,1,1,1),
                        error.model=c("combined","combined"))
saemix.options<-list(seed=1234,save=FALSE,save.graphs=FALSE)
pkpd.fit <- saemix(pkpd.model, pkpd.data, saemix.options)


###################################
# Bootstrap for multiple responses (Gaussian)
nboot <- 2
case.count <- saemix.bootstrap(pkpd.fit, method="case", nboot=nboot) 
case.count

pkpd.fit@options$warnings<-TRUE
cond.count <- saemix.bootstrap(pkpd.fit, method="conditional", nboot=nboot) 
cond.count

####################################################################
# Multiple response model - 1 Gaussian + 1 binary response
head(datPKPD.transcov)
datPKPD.bin <- datPKPD.transcov
datPKPD.bin$Obs[datPKPD.bin$ytype==2]<-ifelse(datPKPD.bin$Obs[datPKPD.bin$ytype==2]>3,1,0)

pkpdmodel.bin <- function(psi,id,xidep) { 
  tim<-xidep[,1]  
  dose<-xidep[,2]
  y<-xidep[,3]
  ytype<-xidep$ytype
  y<-y[ytype==2]
  CL<-psi[id,1]
  V<-psi[id,2]
  inter<-psi[id,3]
  slope<-psi[id,4]
  k<-CL/V
  ypred<-dose/V*(exp(-k*tim))
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  logpdf<-rep(0,length(tim[ytype==2]))
  P.obs = (y==0)*(1-pevent[ytype==2])+(y==1)*pevent[ytype==2]
  ypred[ytype==2] <- log(P.obs)
  return(ypred)
}


pkpdmodel.simbin <- function(psi,id,xidep) { 
  tim<-xidep[,1]  
  dose<-xidep[,2]
  ytype<-xidep$ytype
  CL<-psi[id,1]
  V<-psi[id,2]
  inter<-psi[id,3]
  slope<-psi[id,4]
  k<-CL/V
  ypred<-dose/V*(exp(-k*tim))
  logit<-inter+slope*tim
  pevent<-exp(logit[ytype==2])/(1+exp(logit[ytype==2]))
  ysim<-rbinom(length(pevent),size=1, prob=pevent)
  ypred[ytype==2] <- ysim
  return(ypred)
}



pkpd.bindata<-saemixData(name.data=datPKPD.bin,
                      name.group=c("Id"),name.predictors=c("Time","Dose","Obs"),name.ytype="ytype",
                      name.response=c("Obs"),name.covariates=c("lWeight","lAge","Sex"),
                      units=list(x="hr",y=c("mg/L","-"),covariates=c("logkg","logyr","-")))

pkpd.binmodel<-saemixModel(model=pkpdmodel.bin,modeltype=c("structural","likelihood"),
                        description="PK/PD model with binary response", 
                        psi0=matrix(c(psi0[1:2],1,0,0.75,1,0,0),ncol=4, byrow=TRUE,
                                    dimnames=list(NULL, c("CL","V","intercept","slope"))),transform.par=c(1,1,0,0),
                        error.model=c("combined"))
saemix.options<-list(seed=1234,save=FALSE,save.graphs=FALSE)
pkpd.binfit <- saemix(pkpd.binmodel, pkpd.bindata, saemix.options)


###################################
# Bootstrap for multiple responses (Gaussian)
case.count <- saemix.bootstrap(pkpd.binfit, method="case", nboot=nboot) 
case.count

pkpd.binfit@options$warnings<-TRUE
pkpd.binfit@model@simulate.function<-pkpdmodel.simbin
cond.count <- saemix.bootstrap(pkpd.binfit, method="conditional", nboot=nboot) 
cond.count

####################################################################
