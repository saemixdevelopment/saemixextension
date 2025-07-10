library(saemix)
library(ggplot2)

saemixDir <- "/home/eco/work/saemix/saemixextension"
source(file.path(saemixDir,"newCode","checkInitialFixedEffects.R"))

############################
# Data and model function
data(knee.saemix)

knee.data<-saemixData(name.data=knee.saemix,name.group=c("id"),
                      name.predictors=c("y", "time"), name.X=c("time"),
                      name.covariates = c("Age","Sex","treatment","Age2"),
                      units=list(x="d",y="", covariates=c("yr","-","-","yr2")), verbose=FALSE)

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
  pobs = (y==1)*pge1 + (y==2)*(pge2 - pge1)+(y==3)*(pge3 - pge2)+(y==4)*(pge4 - pge3)+(y==5)*(1 - pge4)
  logpdf <- log(pobs)
  
  return(logpdf)
}
# Multinomial model with cutoff for negative probabilities
ordinal.model.multinom2<-function(psi,id,xidep) {
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
  # to avoid negative probabilities
  p2<-cutoff.eps(pge2 - pge1)
  p3<-cutoff.eps(pge3 - pge2)
  p4<-cutoff.eps(pge4 - pge3)
  p5<-cutoff.eps(1 - pge4)
  pobs = (y==1)*cutoff.eps(pge1) + (y==2)*p2+(y==3)*p3+(y==4)*p4+(y==5)*p5
  logpdf <- log(pobs)
  
  return(logpdf)
}

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, nb.chains=10, nbiter.saemix=c(600,100), print=FALSE)

############################
# Saemix model - not working, tried with many different combinations of psi0 and omega.init => check initialisation ??
xpsi <- c(-1.8, 0.4, 0.7, -0.03, 0.2,-0.05, -0.1, -0.2)
xpsi <- c(-1.8, 0.4, 0.7, -0.03, -0.1,-0.1, -0.1, -0.2)
xpsi <- c(-1.8, 0.4, 0.7, -0.03, 0.1,0.1, 0.1, 0.2)
xpsi <- c(1.8, 0.4, 0.7, 0.3, 0.3,-0.1, -0.1, -0.2)
xpsi <- c(-1.8, -3, -1.8, -0.3, 0.1,0.2,0.1,-0.1)
xpsi <- c(-28, -4, -1.6, 0.5, 0.05, 0.35, 0.89, 0.23) # Monolix, logit=NA for time=7 and 10 when y=4
xpsi <- c(-28, -4, -1.6, 0.5, 0.05, 0.35, 0.4, 0.23) # Monolix, changing beta3 to 0.4 instead of 0.89, now all probabilities are non negative for the population values, but the model still won't start fitting

pcalc <- c(xpsi[1:4], xpsi[1:4]+xpsi[5:8]*10)
cat(1/(1+exp(-pcalc)),"\n")

# With the first version, some probabilities are negative at baseline or the run can't start
multinom.model<-saemixModel(model=ordinal.model.multinom,description="Multinomial model",modeltype="likelihood",
                            psi0=matrix(xpsi,ncol=8, byrow=TRUE, dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta1","beta2","beta3", "beta4"))), 
                            transform.par=c(rep(0,4), rep(1,4)),
                            omega.init=diag(c(10, 10,10,10, 1,1,1,1)), covariance.model = diag(c(1,0,0,0,1,1,0,0)), verbose=FALSE)
multinom.fit<-saemix(multinom.model,knee.data,saemix.options)

# With the second version (tolerance), the run starts but LL can't be computed
multinom.model2<-saemixModel(model=ordinal.model.multinom2,description="Multinomial model",modeltype="likelihood",
                            psi0=matrix(xpsi,ncol=8, byrow=TRUE, dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta1","beta2","beta3", "beta4"))), 
                            transform.par=c(rep(0,4), rep(1,4)),
                            omega.init=diag(c(10, 10,10,10, 1,1,1,1)), covariance.model = diag(c(1,0,0,0,1,1,0,0)), verbose=FALSE)
multinom.fit<-saemix(multinom.model2,knee.data,saemix.options)
yfit <- conddist.saemix(multinom.fit)
yfit<-saemix.predict(yfit)

# can't compute LL by importance sampling
# Checking predictions - ok
summary(multinom.fit@results@predictions)

# Computing probabilities for the estimated population parameters: -Inf/NA
lfit<-checkInitialFixedEffects(multinom.model, knee.data, psi=multinom.fit@results@fixed.effects)
table(knee.data@data$y[is.na(lfit)], knee.data@data$time[is.na(lfit)])

###################################
# Checking values returned for the population parameters

ypred<-checkInitialFixedEffects(multinom.model, knee.data)

# Finding out where pb is 
knee.data@data[is.na(ypred),]
table(knee.data@data$y[!is.na(ypred)], knee.data@data$time[!is.na(ypred)])
table(knee.data@data$y[is.na(ypred)], knee.data@data$time[is.na(ypred)])

mydat <- knee.data
xidep <- mydat@data[,mydat@name.predictors]
psi<-do.call(rbind,rep(list(xpsi),mydat@N))
ypred1<-ordinal.model.multinom(psi, mydat@data[,"index"], xidep)

xidep1<-xidep[1:4,]
psi1 <- psi[1,,drop=FALSE]
ypred1<-ordinal.model.multinom(psi1, rep(1,4), xidep1)

# By hand => some probabilities (y=4 at time 7 and 10) become negative
## we can code differently to let them go down to a tiny number (tolerance)
xtim<-c(0,3,7,10)
logit1<-xpsi[1] + xpsi[5]*xtim
logit2<-xpsi[2] + xpsi[6]*xtim
logit3<-xpsi[3] + xpsi[7]*xtim
logit4<-xpsi[4] + xpsi[8]*xtim
exp(logit3)/(1+exp(logit3))
exp(logit4)/(1+exp(logit4))
data.frame(xtim, p1=exp(logit1)/(1+exp(logit1)), p2=exp(logit2)/(1+exp(logit2))-exp(logit1)/(1+exp(logit1)),
           p3=exp(logit3)/(1+exp(logit3))-exp(logit2)/(1+exp(logit2)), p4=exp(logit4)/(1+exp(logit4)) - exp(logit3)/(1+exp(logit3)), 
           p5=1-exp(logit4)/(1+exp(logit4)))

###################################
# Other initial values for omega
multinom.model<-saemixModel(model=ordinal.model.multinom,description="Multinomial model",modeltype="likelihood",
                            psi0=matrix(xpsi,ncol=8, byrow=TRUE, dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta1","beta2","beta3", "beta4"))), 
                            transform.par=c(rep(0,4), rep(1,4)),
                            omega.init=diag(c(10, 10,10,10, 0.1,0.1,0.1,0.1)), covariance.model = diag(c(1,0,0,0,1,1,0,0)), verbose=FALSE)
multinom.model<-saemixModel(model=ordinal.model.multinom,description="Multinomial model",modeltype="likelihood",
                            psi0=matrix(xpsi,ncol=8, byrow=TRUE, dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta1","beta2","beta3", "beta4"))), 
                            transform.par=c(rep(0,4), rep(1,4)),
                            omega.init=diag(c(1,1,1,1, 1,1,1,1)), covariance.model = diag(c(1,0,0,0,1,1,0,0)), verbose=FALSE)

