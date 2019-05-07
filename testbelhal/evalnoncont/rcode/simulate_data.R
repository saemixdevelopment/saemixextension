# Input
# parpop, omega, respar - population fixed and random effects
# xidep - design
# datDir, namsimdat - names of directory and scenario

if(iscenar %in% c(1:9)) {
  library(mvtnorm)
  for(isim in 1:nsim) {
    etas<-rmvnorm(nsuj, mean=c(0,0), sigma=omega[2:3,2:3])
    etas<-cbind(rnorm(nsuj,mean=0,sd=sqrt(omega[1,1])),etas)
    parsuj<-cbind(etas,rep(parpop[4],nsuj))
    for(i in 1:3) parsuj[,i]<-parpop[i]*exp(parsuj[,i])
    if(iscenar==1) {
      psiM<-parsuj[,1:3]
      nampar<-nampar[1:3]
      } else psiM<-parsuj[,1:4]
    idM<-xidep[,1]
    xiM<-xidep[,2,drop=FALSE]
    fpred<-modfun(psiM,idM,xiM)
    gpred<-fpred*rnorm(length(idM),mean=0,sd=respar)
    xsim<-data.frame(id=idM,dose=xiM,y=fpred+gpred)
    xpar<-data.frame(id=1:nsuj,psiM)
    colnames(xpar)[-c(1)]<-nampar
    
    namfich<-paste('data_',namsimdat,isim,".tab",sep="")
    write.table(xsim,file.path(datDir,namfich),quote=FALSE,row.names=FALSE)
    namfich<-paste('param_',namsimdat,isim,".tab",sep="")
    write.table(xpar,file.path(datDir,namfich),quote=FALSE,row.names=FALSE)
  }
}
