xarr<-array(1:24,dim=c(3,4,2))

xarr
apply(xarr,c(1,3),mean)

##################################################### 09/2024
library(testthat)
library(ggplot2)
library(rlang)
saemixDir<-"/home/eco/work/saemix/saemixextension"
setwd(saemixDir)
library(saemix)

################## Theophylline example
source(file.path(saemixDir,"testeco","helper-source.R"))
datDir<-file.path(saemixDir,"data")

source(file.path(saemixDir,"testeco","test_setup_cont.R"))

# Graphs
theo.fit1@options$displayProgress<-TRUE
theo.fit1@options$warnings<-TRUE
dcond2<-conddist.saemix(theo.fit1,nsamp = 200)

dim(dcond2@results@psi.samp)

ypl<-NULL
for(isuj in 1:3) {
  ypl0<-data.frame(id=isuj,ysamp=c(dcond2@results@psi.samp[isuj,1,]),param="ka")
  ypl1<-data.frame(id=isuj,ysamp=c(dcond2@results@psi.samp[isuj,2,]),param="V")
  ypl2<-data.frame(id=isuj,ysamp=c(dcond2@results@psi.samp[isuj,3,]),param="CL")
  ypl3<-rbind(ypl0,ypl1,ypl2)
  ypl3<-cbind(ypl3,method="revised")
}

ggplot(ypl3,aes(ysamp,fill=as.factor(method),colour=as.factor(method)))+geom_density(alpha=0.4) + facet_wrap(id~param, scales="free")  + theme(legend.position = "top") 


# Debugging
saemixObject<-theo.fit1

################## Leo example
leoDir <- "/home/eco/work/saemix/hotline/leoShiny2407"

source(file.path(leoDir,"leo_predictNewdata.R")) # modifié pour autoriser un seul sujet

# returns lemaitre.fit: fit on simulated data
source(file.path(leoDir,"leo_levofloxacin.R"))

# saemix.fit: replace population parameters with the initial ones to use as prior for the EBE estimates on the Test Patient
saemix.fit<-lemaitre.fit
psi0.pop<-c(psi.lemaitre[1:4],beta.cov[1:2],psi.lemaitre[5:6])
saemix.fit@results@fixed.effects<-saemix.fit@results@fixed.psi<-psi0.pop
saemix.fit@results@omega<-omega.lemaitre
saemix.fit@results@respar<-sig.lemaitre
saemix.fit@results@betas<-t(c(log(psi.lemaitre[1:4]),beta.cov[1:2],log(psi.lemaitre[5:6])))
for (icol in 1:length(psi.lemaitre)) saemix.fit@results@mean.phi[,icol]<-psi.lemaitre[icol]
saemix.fit@options$warnings<-TRUE # slightly more info, in particular the nb of iterations

# Test Patient
ID<-c(1,1)
DOSE<-c(500,500)
TIME<-c(12,3)
CONCENTRATION<-c(6,40)
DFG<-c(120,120)
SEX<-c(1,1)
newdata<-data.frame(Id=ID,Dose=DOSE,Time=TIME,Concentration=CONCENTRATION,DFG=DFG,Sex=SEX)
newdata<-cbind(newdata,tau=12,lDFG=log(newdata$DFG/100))

# Predict
xpred<-saemixPredictNewdata(saemix.fit, newdata, type=c("ipred", "ypred", "ppred", "icpred"),nsamp=200)

# Conditional distribution on the simulated data
lemaitre.fit@options$displayProgress<-TRUE
lemaitre.fit@options$warnings<-TRUE
dcond1 <- conddist.saemix(lemaitre.fit, nsamp=200)
dim(dcond1@results@psi.samp)

# Conditional distribution with the population parameters
saemix.fit@options$displayProgress<-TRUE
saemix.fit@options$warnings<-TRUE
dcond2 <- conddist.saemix(saemix.fit, nsamp=200)

dim(dcond2@results@psi.samp)

dcond3 <- conddist.saemix(saemix.fit, nsamp=200, max.iter=200)

ypl<-NULL
for(isuj in 1:dim(dcond1@results@psi.samp)[1]) {
  ypl3<-NULL
  for(ipar in 1:length(saemix.fit@model@name.modpar)) {
    ypl0<-data.frame(id=isuj,ysamp=c(dcond1@results@psi.samp[isuj,ipar,]),param=saemix.fit@model@name.modpar[ipar])
    ypl3<-rbind(ypl3, ypl0)
  }
  ypl3<-cbind(ypl3,method="reestimated parameters")
  ypl<-ypl3
  ypl3<-NULL
  for(ipar in 1:length(saemix.fit@model@name.modpar)) {
    ypl0<-data.frame(id=isuj,ysamp=c(dcond2@results@psi.samp[isuj,ipar,]),param=saemix.fit@model@name.modpar[ipar])
    ypl3<-rbind(ypl3, ypl0)
  }
  ypl3<-cbind(ypl3,method="population parameters")
  ypl<-rbind(ypl, ypl3)
  ypl3<-NULL
  for(ipar in 1:length(saemix.fit@model@name.modpar)) {
    ypl0<-data.frame(id=isuj,ysamp=c(dcond3@results@psi.samp[isuj,ipar,]),param=saemix.fit@model@name.modpar[ipar])
    ypl3<-rbind(ypl3, ypl0)
  }
  ypl3<-cbind(ypl3,method="population parameters, maxiter=200")
  ypl<-rbind(ypl, ypl3)
}

ggplot(ypl,aes(ysamp,fill=as.factor(method),colour=as.factor(method)))+geom_density(alpha=0.4) + facet_wrap(id~param, scales="free")  + theme(legend.position = "top")
# As expected, parameter distributions closer to priors when we reset them
# no difference with more iterations, but changed the code to have at least 50

# Debugging distcond
saemixObject<-saemix.fit
# Resultat:
# quand on le fait pour les 3 sujets simulés, on a des estimations de paramètres qui restent dans le range => ok

# Debugging predictNewdata
saemixObject<-saemix.fit
newdata
type<-c("ipred", "ypred", "ppred", "icpred")
nsamp<-200

################## Checking parameters and their range after distcond
nsamp<-100
max.iter<-NULL
displayPlot <- TRUE
# executer code dans la fonction distcond

# k final, on refait le calcul avec les 50 itérations précédentes 
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
# ne change rien (avec econd)
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

# Moyenne des phi.samp
apply(phi.samp,c(1,2),mean)
cond.mean.phi

# Check that all the estimates are contained within the corresponding quantiles => ok
qlist <- apply(phi.samp,c(1,2),quantile,c(0.025,0.975))
nampar<-saemix.fit@model@name.modpar
nsuj<-dim(cond.mean.phi)[1]
qpl<-data.frame(id=rep(1:nsuj, length(nampar)), value=c(cond.mean.phi), param=rep(nampar, each=nsuj))
qpl1<-NULL
for(i in 1:length(nampar)) qpl1 <-rbind(qpl1,
              data.frame(q025=qlist[,,i][1,],q975=qlist[,,i][2,]))
qpl<-cbind(qpl, qpl1)

ggplot(qpl, aes(x=id, y=value)) + geom_point() + geom_segment(aes(x=id, xend=id, y=q025, yend=q975)) + facet_wrap(vars(param), scales="free")

eik1 <- eik*0

for(kk in (k-50):(k-1)) {
  eik1 <- eik1
  kk0<-kk0+1
  ijk<-(k-1)*N
  eik01<-eik0
  eik0<-eik0*(kk0-1)/kk0+econd[,(ijk+1):(ijk+N),]/kk0
  varik0<-varik0*(kk0-1)/kk0+((econd[,(ijk+1):(ijk+N),])**2)/kk0 + (eik01**2)
}
summary(c(eik0-eik))

################## Checking parameters and their range after estimateIndividualParametersNewdata
qlist <- apply(saemixObject@results@psi.samp,c(1,2),quantile,c(0.025,0.5,0.975))
xtab <- data.frame(param=saemixObject@model@name.modpar, cond.mean=c(saemixObject@results@cond.mean.psi), mean.samp=c(apply(saemixObject@results@psi.samp,c(1,2),mean)))
xt1<-NULL
for (i in 1:dim(xtab)[1])
  xt1<-rbind(xt1, qlist[,,i])
xtab<-cbind(xtab, xt1)
print(xtab)
# cond.mean different from conditional means of samples
#  conditional means of samples is close to median of samples BUT still below 2.5% for Tlag !

tlag200 <- saemixObject@results@psi.samp[1,2,c(1:200)]
hist(tlag200, breaks=50)

##################################################### 03/2020
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
