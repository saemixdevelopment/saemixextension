xarr<-array(1:24,dim=c(3,4,2))

xarr
apply(xarr,c(1,3),mean)

#####################################################
library(testthat)
library(ggplot2)
saemixDir<-"/home/eco/work/saemix/saemixextension"
setwd(saemixDir)
source(file.path(saemixDir,"testeco","helper-source.R"))
datDir<-file.path(saemixDir,"data")

source(file.path(saemixDir,"testeco","test_setup_cont.R"))

# Debugging conditional distribution function
source('~/work/saemix/newLib/saemix/R/func_distcond.R') # old version, conddist1
source('~/work/saemix/saemixextension/R/func_distcond.R') # new version

# Graphs
theo.fit1@options$displayProgress<-TRUE
theo.fit1@options$warnings<-TRUE
dcond1<-conddist1.saemix(theo.fit1,nsamp = 200)
dcond2<-conddist.saemix(theo.fit1,nsamp = 200)

dim(dcond1@results@psi.samp)
dim(dcond2@results@psi.samp)

ypl<-NULL
for(isuj in 1:3) {
  ypl0<-data.frame(id=isuj,ysamp=c(dcond1@results@psi.samp[isuj,1,]),param="ka")
  ypl1<-data.frame(id=isuj,ysamp=c(dcond1@results@psi.samp[isuj,2,]),param="V")
  ypl2<-data.frame(id=isuj,ysamp=c(dcond1@results@psi.samp[isuj,3,]),param="CL")
  ypl4<-rbind(ypl0,ypl1,ypl2)
  ypl4<-cbind(ypl4,method="original")
  ypl<-rbind(ypl,ypl4)
  
  ypl0<-data.frame(id=isuj,ysamp=c(dcond2@results@psi.samp[isuj,1,]),param="ka")
  ypl1<-data.frame(id=isuj,ysamp=c(dcond2@results@psi.samp[isuj,2,]),param="V")
  ypl2<-data.frame(id=isuj,ysamp=c(dcond2@results@psi.samp[isuj,3,]),param="CL")
  ypl3<-rbind(ypl0,ypl1,ypl2)
  ypl3<-cbind(ypl3,method="revised")
  
  ypl<-rbind(ypl,ypl3)
}

ggplot(ypl,aes(ysamp,fill=as.factor(method),colour=as.factor(method)))+geom_density(alpha=0.4) + facet_wrap(id~param, scales="free")  + theme(legend.position = "top") 

# Debugging
saemixObject<-theo.fit1
nsamp<-100
max.iter<-NULL


# k final, on refait le calcul avec les 50 itérations précédentes => ne change rien
eik0<-eik*0
varik0<-eik0

kk0<-0
for(kk in (k-50):(k-1)) {
  kk0<-kk0+1
  ijk<-(k-1)*N
  eik01<-eik0
  eik0<-eik0*(kk0-1)/kk0+econd[,(ijk+1):(ijk+N),]/kk0
  varik0<-varik0*(kk0-1)/kk0+((econd[,(ijk+1):(ijk+N),])**2)/kk0 + (eik01**2)
}
summary(c(eik0-eik))
