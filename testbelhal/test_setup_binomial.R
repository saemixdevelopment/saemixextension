###################################################################################
cat("Simulating binary example with covariate\n")

nsuj<-1000
xtim<-c(0:3)
parnam<-c("Intercept","beta.time","beta.RF")
param<-c(0,-0.37,0.85)
omega<-c(.21,.1,0)

partab<-as.data.frame(matrix(data=0,nrow=nsuj,ncol=3,dimnames=list(NULL,parnam)))
partab$beta.RF[(nsuj/2+1):nsuj]<-param[3]
for(i in 1:2) partab[,i]<-rnorm(nsuj,mean=param[i],sd=omega[i])

psim<-data.frame()
for(itim in xtim) {
  logit.sim<-partab[,1]+partab[,2]*itim+partab[,3]
  xtab<-exp(logit.sim)/(1+exp(logit.sim))
  psim<-rbind(psim,xtab)
}
datsim<-data.frame(id=rep(1:nsuj,each=length(xtim)),time=rep(xtim,nsuj),psim=unlist(psim))
rownames(datsim)<-NULL
ysim<-rbinom(nsuj*length(xtim),size=1,prob=datsim$psim)
summary(datsim)
datsim$y<-ysim
datsim$risk<-ifelse(datsim$id>500,1,0)

if(FALSE) {
  namfich<-file.path(datDir,"simulatedBinary.txt")
  write.table(datsim,namfich,row.names=F,quote=F)
  namfich<-file.path(datDir,"simulatedBinary_par.txt")
  write.table(partab,namfich,row.names=T,quote=F)
}

# Running saemix
saemix.data<-saemixData(name.data=file.path(datDir,"simulatedBinary.txt"),
      name.group=c("id"),name.predictors=c("time","y"), 
      name.covariates=c("risk"),name.X=c("time"))

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

saemix.model<-saemixModel(model=binary.model,description="Binary model",
        modeltype="likelihood",
        psi0=matrix(c(0,-.5,0.5,0),ncol=2,byrow=TRUE,dimnames=list(NULL,parnam[1:2])),
        transform.par=c(0,0),covariance.model=matrix(c(1,0,0,1),ncol=2))

model.cov<-saemixModel(model=binary.model,description="Binary model",
        modeltype="likelihood",
        psi0=matrix(c(0,-.5,0.5,0),ncol=2,byrow=TRUE,dimnames=list(NULL,parnam[1:2])),covariate.model=matrix(c(1,0),ncol=2),
        transform.par=c(0,0),covariance.model=matrix(c(1,0,0,1),ncol=2))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
# saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

binary.fit<-saemix(model.cov,saemix.data,saemix.options)

if(FALSE) {
# Check with contingency table that predictions match observations
  vec<-predict(binary.fit)
  ypred<-as.numeric(exp(vec)>=0.5)
  table(ypred,datsim$y)
  fisher.test(table(ypred,datsim$y))
  # not stellar :-(

# Compare predicted and simulated individual parameters
  parsim<-read.table(file.path(datDir,"simulatedBinary_par.txt"),header=T)
  head(binary.fit@results@map.psi)
  head(parsim)
  par(mfrow=c(1,2))
  for(i in 1:2) {
    vec<-parsim[,i]
    if(i==1) vec<-vec+parsim[,3]*param[3]
    plot(vec,binary.fit@results@map.psi[,i],xlab="Simulated",ylab="Estimated", main=parnam[i],pch=20)
    abline(0,1)
    print(cor.test(binary.fit@results@map.psi[,i],vec))
  }
  par(mfrow=c(1,2))
  for(i in 1:2) {
    vec<-parsim[,i]
    if(i==1) vec<-vec+parsim[,3]*param[3]
    plot(vec,binary.fit@results@cond.mean.psi[,i],xlab="Simulated",ylab="Estimated", main=parnam[i],pch=20)
    abline(0,1)
    print(cor.test(binary.fit@results@cond.mean.psi[,i],vec))
  }
  par(mfrow=c(1,2))
  for(i in 1:2) {
    plot(binary.fit@results@map.psi[,i],binary.fit@results@cond.mean.psi[,i],xlab="Simulated",ylab="Estimated", main=parnam[i],pch=20)
    abline(0,1)
    print(cor.test(binary.fit@results@cond.mean.psi[,i],binary.fit@results@map.psi[,i]))
  }
}
