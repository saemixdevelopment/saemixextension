#### Joint model example : prothro dataset in R package JM
library(JM)
library(pracma)
library(ggplot2)
data("prothro")
data("prothros")

gp = ggplot(data=prothro[which(prothro$id %in% 1:20),], aes(x=time, y=pro, group = id))+geom_point(lwd=1.5)+geom_line(col="#CC0033",lwd=0.8)+theme_classic()+
  ylab("Marker observations")+xlab("Days")+theme(axis.text = element_text(size=14),
                                                 axis.title = element_text(size=16))
gp

table(prothros$death)


######### estimation with saemix 

saemixDir <- "C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension"
workDir <- file.path(saemixDir, "joint")
setwd(workDir)

library(Cairo)
library("viridis")
library(rlang)

# Loading of extended functions of the library 
progDir<-file.path(saemixDir, "R")
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"func_plots.R"))

source(file.path(workDir,"multi_aux.R"))
source(file.path(workDir,"multi_initializeMainAlgo.R"))
source(file.path(workDir,"multi_estep.R"))
source(file.path(workDir,"multi_mstep.R"))
source(file.path(workDir,"multi_main.R"))
source(file.path(workDir,"multi_map.R"))
source(file.path(workDir,"compute_LL_multi.R"))

############################################## Data and model
# Formatting and creating data 
d1 = prothro[,c(1,2,3)]
d1$ytype=1
colnames(d1)[2] = "obs"
d2 = prothros[,c(1,3,2)]
d2$ytype = 2
colnames(d2)[2] = "obs"
colnames(d2)[3] = "time"
data_joint = rbind(d1,d2)
dataJM<-saemixData(name.data=data_joint, name.group=c("id"), name.predictors=c("time","obs"),
                   name.response="obs",name.ytype = "ytype")



#### initializing parameters

param<-c(73,1.25,0.6,0.0001)
omega.sim<-c(18, 3, 0.05, 0.01)
sigma.sim <- 17

# Model definition 

JMmodel<-function(psi,id,xidep) {
  ytype<-xidep$ytype  # type of response (1: continuous, 2: event)
  b0 <- psi[id,1] 
  b1 <- psi[id,2] 
  h0 <- psi[id,3]
  alpha <- psi[id,4] 
  
  ypred <- b0+b1*xidep[,1]   
  
  T<-xidep[ytype==2,1]# vector of times (survival part)
  Nj <- length(T)
  ev = xidep$obs[ytype==2]
  cens<-which(ev==0) 
  ind <- which(ev==1)
  b0b = b0[ytype==2] # to have vectors of the same length as T 
  b1b = b1[ytype==2]
  h0b = h0[ytype==2]
  alphab = alpha[ytype==2]
  
  haz <- h0b*exp(alphab*(b0b+b1b*T)) # instantaneous hazard
  H <- (h0b/(alphab*b1b))*exp((b0b+b1b*T)*alphab)-(h0b/(alphab*b1b))*exp(alphab*b0b) # cumulative hazard
  
  logpdf <- rep(0,Nj)
  logpdf[cens] <- -H[cens] 
  logpdf[ind] <- -H[ind] + log(haz[ind]) 
  
  ypred[ytype==2] = logpdf
  return(ypred)
}


jointTTE<-saemixModel(model=JMmodel,description="JM LMEM-TTE (prothro data)",modeltype=c("structural","likelihood"),
                      psi0=matrix(param,ncol=4,byrow=TRUE,dimnames=list(NULL, c("b0","b1","h0","alpha"))),
                      transform.par=c(0,0,1,0), covariance.model=diag(c(1,1,0,0)),
                      fixed.estim = c(1,1,1,1),error.model = "constant",
                      omega.init = diag(omega.sim))

saemix.data<-dataJM
saemix.model<-jointTTE
saemix.options<-saemixControl(seed=12345, map=T, fim=T, ll.is=TRUE, save.graphs = F) # please, specify save.graphs=F (currently not extended)
yfit <- saemix.multi(saemix.model, saemix.data, saemix.options)

summary(yfit) # parameter estimates + likelihood 
yfit@results@fim # inverse of the Fisher Information Matrix 
sqrt(diag(yfit@results@fim)) # SE of parameter estimates 


################### WITH a nonlinear mixed-effects model ###################
# time consuming because the model as no explicit likelihood expression...
# use of numerical integration to compute joint likelihood 

param<-c(20,2,0.6,0.5,2,-0.039)
omega.sim<-c(5, 1, 0.6,0.5,0.02, 0.03)
sigma.sim <- 17

# model
JMmodel_nl<-function(psi,id,xidep) {
  ytype<-xidep$ytype  
  b0 <- psi[id,1] 
  b1 <- psi[id,2] 
  b2 <- psi[id,3]
  a <- psi[id,4]
  h0 <- psi[id,5]
  alpha <- psi[id,6]
  
  ypred <- b0+a*(exp(-b1*xidep[,1])-exp(-b2*xidep[,1]))  
  
  T<-xidep[ytype==2,1]
  Nj <- length(T)
  ev = xidep$obs[ytype==2]
  cens<-which(ev==0)  
  ind <- which(ev==1)
  b0b = b0[ytype==2]
  b1b = b1[ytype==2]
  b2b = b2[ytype==2]
  ab = a[ytype==2]
  h0b = h0[ytype==2]
  alphab = alpha[ytype==2]
  
  f=function(x) seq(0,x,length.out=100)
  tab = mapply(f,T)
  tab = t(tab)
  pas = tab[,2]-tab[,1]
  
  haz <- h0b*exp(alphab*(b0b+ab*(exp(-b1b*T)-exp(-b2b*T))))
  hazt <- h0b*exp(alphab*(b0b+ab*(exp(-b1b*tab)-exp(-b2b*tab))))
  H = apply(hazt,1,sum)*pas   # no close form... need to integrate numerically 
  
  logpdf <- rep(0,Nj)
  logpdf[cens] <- -H[cens] 
  logpdf[ind] <- -H[ind] + log(haz[ind])
  
  ypred[ytype==2] = logpdf
  return(ypred)
}

jointTTE_nl<-saemixModel(model=JMmodel_nl,description="JM lin longi one tte",modeltype=c("structural","likelihood"),
                      psi0=matrix(param,ncol=6,byrow=TRUE,dimnames=list(NULL, c("b0","b1","b2","a","h0","alpha"))),
                      transform.par=c(0,1,1,1,0,0), covariance.model=diag(c(1,1,1,1,0,0)),
                      fixed.estim = c(1,1,1,1,1,1),error.model = "constant",
                      omega.init = diag(omega.sim))

saemix.model_nl<-jointTTE_nl
yfit_nl <- saemix.multi(saemix.model_nl, saemix.data, saemix.options)

yfit_nl@results@fim
sqrt(diag(yfit_nl@results@fim))
