# Input
# resDir, namsimdat - names of directory and scenario
# parpop, omega, respar - population fixed and random effects

namOrig<-file.path(resDir,paste("parpop_",namsimdat,".res",sep=""))
namSE<-file.path(resDir,paste("sepop_",namsimdat,".res",sep=""))

partab<-read.table(namOrig,header=T)
parSE<-read.table(namSE,header=T)
npar<-dim(parSE)[2]-1

simfix<-parpop
vec<-unlist(omega[lower.tri(omega)])
simomega<-c(diag(omega),vec[vec!=0])
simse<-respar
popsim<-c(simfix,simse,simomega)

empSE<-apply(partab,2,sd)
empSE<-empSE[-c(1,npar+2,npar+3)]

idx0<-which(parpop==0)
parsim0<-popsim
parsim0[idx0]<-1

nsim<-dim(partab)[1]
tab.comp<-do.call(rbind,rep(list(popsim),nsim))
tab.comp0<-do.call(rbind,rep(list(parsim0),nsim))

diffpar<-(partab[-c(1,npar+2,npar+3)]-tab.comp)/tab.comp0
rbias<-apply(diffpar,2,mean)
rmse<-apply(diffpar,2,sd)

tabempSE<-do.call(rbind,rep(list(empSE),nsim))
tabSE<-(parSE[-c(1)]-tabempSE)/tabempSE
rbiasSE<-apply(tabSE,2,mean)
rmseSE<-apply(tabSE,2,sd)
