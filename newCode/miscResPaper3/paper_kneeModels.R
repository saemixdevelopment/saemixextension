workDir<-"/home/eco/work/saemix/saemixextension/paperSaemix3"
saemixDir <- "/home/eco/work/saemix/saemixextension"
setwd(workDir)

# Libraries
library(saemix)

data(knee.saemix)

knee.data<-saemixData(name.data=knee.saemix,name.group=c("id"),
                        name.predictors=c("y", "time"), name.X=c("time"),
                        name.covariates = c("Age","Sex","treatment","Age2"),
                        units=list(x="d",y="", covariates=c("yr","-","-","yr2")), verbose=FALSE)


# Transition matrices - nobody worsens (odd ?)
ncat<-length(unique(knee.saemix$y))
transmat <- matrix(data=0, nrow=ncat, ncol=ncat)
for(i in 2:dim(knee.saemix)[1]) {
  if(knee.saemix$id[i]==knee.saemix$id[i-1]) { # same subject
    transmat[knee.saemix$y[i-1],knee.saemix$y[i]] <- 1+transmat[knee.saemix$y[i-1],knee.saemix$y[i]]
  }
}
transmat

diag(transmat)/rowSums(transmat)
table(knee.saemix$y)
table(knee.saemix$y[knee.saemix$time==0])/length(knee.saemix$y[knee.saemix$time==0])
table(knee.saemix$y[knee.saemix$time==10])/length(knee.saemix$y[knee.saemix$time==0])

################################################################## Categorical models
# currently only PO works, the other ones won't start

#################################
# Fitting PO model - ok but fit not very good (see diagnostics)
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

# Saemix model
saemix.model.PO<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                          simulate.function=simulateOrdinal, 
                          psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5, byrow=TRUE, 
                                     dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), transform.par=c(0,1,1,1,1),
                          omega.init=diag(c(100, 1, 1, 1, 1)), covariance.model = diag(c(1,0,0,0,1)), verbose=FALSE)

# Fitting
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, nb.chains=10, nbiter.saemix=c(600,100), print=FALSE)
#saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, nb.chains=10, fim=FALSE)

ord.fit<-saemix(saemix.model.PO,knee.data,saemix.options)
summary(ord.fit)

#################################
# Multinomial model

ordinal.model.multinom<-function(psi,id,xidep) {
  y<-xidep[,1]
  time<-xidep[,2]
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  alp4<-psi[id,4]
  beta1<-psi[id,5]
  beta2<-psi[id,6]
  beta3<-psi[id,7]
  beta4<-psi[id,8]
  
  logit1<-alp1 + beta1*time
  logit2<-alp2 + beta2*time
  logit3<-alp3 + beta3*time
  logit4<-alp4 + beta4*time
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  pge4<-exp(logit4)/(1+exp(logit4))
  pobs = (y==1)*pge1+(y==2)*(pge2 - pge1)+(y==3)*(pge3 - pge2)+(y==4)*(pge4 - pge3)+(y==5)*(1 - pge4)
  logpdf <- log(pobs)
  
  return(logpdf)
}

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, nb.chains=10, nbiter.saemix=c(600,100), print=FALSE)

# Saemix model - not working, tried with many different combinations of psi0 and omega.init => check initialisation ??
xpsi <- c(-1.8, 0.4, 0.7, -0.03, 0.2,-0.05, -0.1, -0.2)
xpsi <- c(-1.8, 0.4, 0.7, -0.03, -0.1,-0.1, -0.1, -0.2)
xpsi <- c(-1.8, 0.4, 0.7, -0.03, 0.1,0.1, 0.1, 0.2)
xpsi <- c(1.8, 0.4, 0.7, 0.3, 0.3,-0.1, -0.1, -0.2)
xpsi <- c(-1.8, -3, -1.8, -0.3, 0.1,0.2,0.1,-0.1)
xpsi <- c(-28, -4, -1.6, 0.5, 0.05, 0.35, 0.89, 0.23) # Monolix

pcalc <- c(xpsi[1:4], xpsi[1:4]+xpsi[5:8]*10)
cat(1/(1+exp(-pcalc)),"\n")
multinom.model<-saemixModel(model=ordinal.model.multinom,description="Multinomial model",modeltype="likelihood",
                            psi0=matrix(xpsi,ncol=8, byrow=TRUE, dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta1","beta2","beta3", "beta4"))), 
                            transform.par=c(rep(0,4), rep(1,4)),
                            omega.init=diag(c(10, 10,10,10, 0.1,0.1,0.1,0.1)), covariance.model = diag(c(1,0,0,0,1,1,0,0)), verbose=FALSE)
multinom.model<-saemixModel(model=ordinal.model.multinom,description="Multinomial model",modeltype="likelihood",
                            psi0=matrix(xpsi,ncol=8, byrow=TRUE, dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta1","beta2","beta3", "beta4"))), 
                            transform.par=c(rep(0,4), rep(1,4)),
                            omega.init=diag(c(10, 10,10,10, 1,1,1,1)), covariance.model = diag(c(1,0,0,0,1,1,0,0)), verbose=FALSE)
multinom.model<-saemixModel(model=ordinal.model.multinom,description="Multinomial model",modeltype="likelihood",
                            psi0=matrix(xpsi,ncol=8, byrow=TRUE, dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta1","beta2","beta3", "beta4"))), 
                            transform.par=c(rep(0,4), rep(1,4)),
                            omega.init=diag(c(1,1,1,1, 1,1,1,1)), covariance.model = diag(c(1,0,0,0,1,1,0,0)), verbose=FALSE)

# Fitting
multinom.fit<-saemix(multinom.model,knee.data,saemix.options)
summary(multinom.fit)

#################################
# Fitting non-PO model
ordinal.model.nonPO<-function(psi,id,xidep) {
  y<-xidep[,1]
  time<-xidep[,2]
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  alp4<-psi[id,4]
  beta1<-psi[id,5]
  beta2<-psi[id,6]

  logit1<-alp1 + beta1*time
  logit2<-alp1 +alp2 + (beta1+beta2)*time
  logit3<-alp1 +alp2 + alp3 + (beta1+beta2)*time
  logit4<-alp1 +alp2 + alp3+ alp4 + (beta1+beta2)*time
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  pge4<-exp(logit4)/(1+exp(logit4))
  pobs = (y==1)*pge1+(y==2)*(pge2 - pge1)+(y==3)*(pge3 - pge2)+(y==4)*(pge4 - pge3)+(y==5)*(1 - pge4)
  logpdf <- log(pobs)
  
  return(logpdf)
}

# Saemix model - not working either
nonPO.model<-saemixModel(model=ordinal.model.nonPO,description="Ordinal categorical model",modeltype="likelihood",
                         psi0=matrix(c(-15,6,8,12, 0.5,0.1),ncol=6, byrow=TRUE, dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta1","beta2"))),
                         transform.par=c(0,1,1,1,1,0),
                         omega.init=diag(c(100, 1, 1, 1, 1,1)), covariance.model = diag(c(1,0,0,0,1,1)), verbose=FALSE)

subnonPO.fit<-saemix(nonPO.model,knee.data,saemix.options)
summary(subnonPO.fit)

################################################################## Binary models
################################# 
# Binarised response for each category, to assess the PO assumption
# note: would need eg bootstrap estimates of the SE to test whether equality of betas is a decent assumption (range between -1.2 and -1.66)

# saemix model
binary.model<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  logpdf<-rep(0,length(tim))
  P.obs = (y==0)*(1-pevent)+(y==1)*pevent
  logpdf <- log(P.obs)
  return(logpdf)
}

# simulation function (used for diagnostics)
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

saemix.binmodel<-saemixModel(model=binary.model,description="Binary model",simulate.function=simulBinary, modeltype="likelihood",
                          psi0=matrix(c(1.4,-0.2,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("inter","slope"))),
                          transform.par=c(0,0), covariance.model=matrix(c(1,0,0,0),ncol=2), omega.init=diag(c(2,0.5)), verbose=FALSE)

for (icat in 1:4) {
  knee1<-knee.saemix
  knee1$y<-ifelse(knee1$y>icat,1,0)
  
  saemix.data<-saemixData(name.data=knee1,name.group=c("id"),
                          name.predictors=c("time","y"), name.X=c("time"),
                          name.covariates = c("Age","Sex","treatment","Age2"),
                          units=list(x="d",y="", covariates=c("yr","-","-","yr2")), verbose=FALSE)
  
  saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE, print=FALSE)
  binary.fit<-saemix(saemix.binmodel,saemix.data,saemix.options)
  cat("-------------------------\n")
  cat("Model for Y>",icat,"\n")
  print(binary.fit)
}

################################# Simple models (non mixed), suggest multinomial better than PO
# MASS to fit polr
library(MASS)
fit.polr <- polr(as.factor(y)~time, data=knee.saemix)

# nnet to fit multinomial (non PO) models
library(nnet)
mlm <- multinom(as.factor(y)~time, data=knee.saemix, Hess=TRUE)

M1 <- logLik(fit.polr)
M2 <- logLik(mlm)
(G <- -2*(M1[1] - M2[1]))
pchisq(G,3,lower.tail = FALSE)
