library(saemix)
nboot<-5

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
  # Model without covariate
  saemix.model<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                            simulate.function=simulateOrdinal, psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5, byrow=TRUE, 
                                                                           dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), transform.par=c(0,1,1,1,1),
                            omega.init=diag(c(100, 1, 1, 1, 1)), covariance.model = diag(c(1,0,0,0,1)), verbose=FALSE)
  saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, nb.chains=10, nbiter.saemix=c(600,100), print=FALSE)
  ord.fit<-saemix(saemix.model,ordknee.data,saemix.options)            
  case.ordinal <- saemix.bootstrap(ord.fit, method="case", nboot=nboot) 
  case.ordinal <- saemix.bootstrap2(ord.fit, method="case", nboot=nboot) 
  
  # Running model again from final estimates: failed to find a valid initial parameter guess :-/
  # so something weird going on...
  saemix.model2<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood", simulate.function=simulateOrdinal, 
                             psi0=matrix(ord.fit@results@fixed.effects*.5,ncol=5, byrow=TRUE, dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), 
                             transform.par=c(0,1,1,1,1),omega.init=diag(c(100, 1, 1, 1, 1)), covariance.model = diag(c(1,0,0,0,1)), verbose=FALSE)
  ord.fit2<-saemix(saemix.model2,ordknee.data,saemix.options)        
  
  # 
  x1<-table(knee.saemix$y[knee.saemix$time==0])
  x1[1]/sum(x1)
  