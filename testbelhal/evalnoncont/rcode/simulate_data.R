# Input
# parpop, omega, respar - population fixed and random effects
# xidep - design
# datDir, namsimdat - names of directory and scenario

if(iscenar %in% c(1:3)) {
  for(isim in 1:nsim) {
    partab<-as.data.frame(matrix(data=0,nrow=nsuj,ncol=2,dimnames=list(NULL,parnam)))
    for(i in 1:2) partab[,i]<-rnorm(nsuj,mean=param[i],sd=omega[i])

    psiM<-data.frame()
    for(itim in xtim) {
      logit.sim<-partab[,1]+partab[,2]*itim
      xtab<-exp(logit.sim)/(1+exp(logit.sim))
      psiM<-rbind(psiM,xtab)
    }
    datsim<-data.frame(id=rep(1:nsuj,each=length(xtim)),time=rep(xtim,nsuj),psiM=unlist(psiM))
    xpar<-data.frame(id=1:nsuj,psiM)
    rownames(datsim)<-NULL
    ysim<-rbinom(nsuj*length(xtim),size=1,prob=datsim$psiM)
    summary(datsim)
    datsim$y<-ysim
    datsim$risk<-ifelse(datsim$id>500,1,0)
    
    namfich<-paste('data_',namsimdat,isim,".tab",sep="")
    write.table(datsim,file.path(datDir,namfich),quote=FALSE,row.names=FALSE)
    namfich<-paste('param_',namsimdat,isim,".tab",sep="")
    write.table(xpar,file.path(datDir,namfich),quote=FALSE,row.names=FALSE)
  }
}

    