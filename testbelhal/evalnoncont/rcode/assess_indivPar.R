# Input
# resDir, namsimdat - names of directory and scenario
# parpop, omega, respar - population fixed and random effects
# nsim - number of simulations

namOrig<-file.path(resDir,paste("parpop_",namsimdat,".res",sep=""))
partab<-read.table(namOrig,header=T)

tab<-tab.shr<-tab.mpe<-NULL
for(isim in 1:nsim) {
  namfich<-paste('param_',namsimdat,isim,".tab",sep="")
  parsim<-read.table(file.path(datDir,namfich),header=TRUE)

  namfich<-paste('indpar_',namsimdat,isim,".tab",sep="")
  parest<-read.table(file.path(resDir,namfich),header=TRUE)
  
  # Estimating mean prediction error and RMSE on concentrations
  namfich<-paste('data_',namsimdat,isim,".tab",sep="")
  saemix.data<-saemixData(file.path(datDir,namfich),header=T,name.group="id",name.predictors="dose",name.response="y",units=list(x="-",y="-"), verbose=FALSE)
  xiM<-saemix.data@data[,saemix.data@name.predictors,drop=FALSE]
  idM<-saemix.data@data[,"index"]
  yobs<-saemix.data@data[,saemix.data@name.response]
  
  parpop<-partab[isim,]
  npar<-length(nampar)
  
  # MAP
  tabs<-parsim[,2:(npar+1)]
  tabe<-parest[,2:(npar+1)]
  tabc<-(tabe-tabs)/tabs
  rbias<-apply(tabc,2,mean)
  rmse<-apply(tabc,2,sd)
  l1<-c(isim,rbias,rmse)
  fpred<-saemix.model@model(tabe,idM,xiM)
#  lobs<-c(isim,mean(yobs-fpred),sd(yobs-fpred),mean((yobs-fpred)/yobs),sd((yobs-fpred)/yobs))
  lobs<-c(isim,mean(yobs-fpred),sd(yobs-fpred),mean((yobs-fpred)/fpred),sd((yobs-fpred)/fpred))
  
  # Conditional mean after estimation
  tabe<-parest[,(npar+2):(2*npar+1)]
  tabc<-(tabe-tabs)/tabs
  rbias<-apply(tabc,2,mean)
  rmse<-apply(tabc,2,sd)
  l1<-c(l1,rbias,rmse)
  fpred<-saemix.model@model(tabe,idM,xiM)
#  lobs<-c(lobs,mean(yobs-fpred),sd(yobs-fpred),mean((yobs-fpred)/yobs),sd((yobs-fpred)/yobs))
  lobs<-c(lobs,mean(yobs-fpred),sd(yobs-fpred),mean((yobs-fpred)/fpred),sd((yobs-fpred)/fpred))
  
  # Conditional mean after conditional distribution
  tabe<-parest[,(2*npar+2):(3*npar+1)]
  tabc<-(tabe-tabs)/tabs
  rbias<-apply(tabc,2,mean)
  rmse<-apply(tabc,2,sd)
  l1<-c(l1,rbias,rmse)
  fpred<-saemix.model@model(tabe,idM,xiM)
  lobs<-c(lobs,mean(yobs-fpred),sd(yobs-fpred),mean((yobs-fpred)/fpred),sd((yobs-fpred)/fpred))

  tab<-rbind(tab,l1)
  tab.mpe<-rbind(tab.mpe,lobs)
  
  # Estimating shrinkage - only works ad hoc, for models without covariates and LN distribution
  tabs<-log(tabs)
  shr.sim<-1-apply(tabs,2,var)/parpop[(npar+3):(npar+2+dim(omega)[1])]
  l2<-c(isim,shr.sim)
  tabe<-parest[,2:(npar+1)]
  tabe<-log(tabe)
  varebe<-apply(tabe,2,var)
  shr.est<-1-varebe/parpop[(npar+3):(npar+2+dim(omega)[1])]
  l2<-c(l2,shr.est)
  tabe<-parest[,(npar+2):(2*npar+1)]
  tabe<-log(tabe)
  varebe<-apply(tabe,2,var)
  shr.est<-1-varebe/parpop[(npar+3):(npar+2+dim(omega)[1])]
  l2<-c(l2,shr.est)
  tabe<-parest[,(2*npar+2):(3*npar+1)]
  tabe<-log(tabe)
  varebe<-apply(tabe,2,var)
  shr.est<-1-varebe/parpop[(npar+3):(npar+2+dim(omega)[1])]
  l2<-c(l2,shr.est)
  l2<-unlist(l2)
  
  tab.shr<-rbind(tab.shr,l2)

}
xmean.bias<-apply(tab,2,mean)
xmean.sd<-apply(tab,2,sd)
y.bias<-apply(tab.mpe,2,mean)
y.sd<-apply(tab.mpe,2,sd)


tab.comp<-tab
# Plot
