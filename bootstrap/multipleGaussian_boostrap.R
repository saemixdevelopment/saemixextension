####################################################################
# Bootstrap for multiple responses - two continuous responses
####################################################################
saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-"/home/eco/work/saemix/versions/lucas/prov/saemixDev/R"
datDir<-file.path(saemixDir,"data")
datDir2 <- "/home/eco/work/saemix/workshop26/data"

# Libraries
library(ggplot2)
library(MASS)
library(rlang)
library(npde)

source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixData-methods.R"))
source(file.path(progDir,"SaemixData-methods_covariates.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"main.R"))
source(file.path(progDir,"func_aux.R"))
source(file.path(progDir,"main_initializeMainAlgo.R")) # ...
source(file.path(progDir,"main_estep.R"))
source(file.path(progDir,"main_mstep.R"))
source(file.path(progDir,"func_FIM.R"))
source(file.path(progDir,"func_map.R"))
source(file.path(progDir,"func_plots.R"))
source(file.path(progDir,"func_distcond.R"))
source(file.path(progDir,"func_simulations.R"))
source(file.path(progDir,"compute_LL.R"))
source(file.path(progDir,"func_npde.R"))
source(file.path(progDir,"func_estimParam.R"))
source(file.path(progDir,"backward.R"))
source(file.path(progDir,"forward.R"))
source(file.path(progDir,"func_stepwise.R"))
source(file.path(progDir,"stepwise.R"))
source(file.path(progDir,"func_compare.R"))
#source(file.path(progDir,"func_bootstrap.R"))
source(file.path(progDir,"func_exploreData.R"))
source(file.path(progDir,"func_discreteVPC.R"))
source(file.path(progDir,"progressBar.R"))

####################################################################
# Bootstrap function, modified
source(file.path(saemixDir,"R","func_bootstrap_mult.R"))
set.seed(42919)

# Number of bootstrap samples
nboot <- 500

####################################################################
# PKPD example

load(file.path(datDir2,"sd_iv_rich_pkpd.rda"))
tab1 <- sd_iv_rich_pkpd[,-c(4,10)]
tab2 <- sd_iv_rich_pkpd[,-c(3,10)]
colnames(tab1)<-colnames(tab2)<-c("Id","Time","Obs","Weight","Age","Dose","Sex","Race")
tab1$ytype <- 1
tab2$ytype <- 2
datPKPD <- rbind(tab1,tab2)
datPKPD <- datPKPD[order(datPKPD$Id, datPKPD$Time, datPKPD$ytype),]
pkpd.data<-saemixData(name.data=datPKPD,
                      name.group=c("Id"),name.predictors=c("Time","Dose"),name.ytype="ytype",
                      name.response=c("Obs"),name.covariates=c("Weight","Age","Sex"),
                      units=list(x="hr",y=c("mg/L","-"),covariates=c("kg","yr","-")))
# Centering covariates
pkpd.data<-transformContCov(pkpd.data,Weight,centering=70,transformation=log,verbose=FALSE, newCovName = "lWeight")
pkpd.data<-transformContCov(pkpd.data,Age,centering="median",transformation=log,verbose=FALSE, newCovName = "lAge")
# Recreating pkpd.data with our new covariates
pkpd.data<-saemixData(name.data=pkpd.data@data,
                      name.group=c("Id"),name.predictors=c("Time","Dose"),name.ytype="ytype",
                      name.response=c("Obs"),name.covariates=c("lWeight","lAge","Sex"),
                      units=list(x="hr",y=c("mg/L","-"),covariates=c("logkg","logyr","-")))

datPKPD.transcov <- pkpd.data@data[,c("Id","Time","Obs","lWeight","lAge","Dose","Sex","ytype")]

pkpdmodel <- function(psi,id,xidep) { 
  tim<-xidep[,1]  
  dose<-xidep[,2]
  ytype<-xidep$ytype
  CL<-psi[id,1]
  V<-psi[id,2]
  Emax<-psi[id,3]
  EC50<-psi[id,4]
  k<-CL/V
  ypred<-dose/V*(exp(-k*tim))
  ypred2 <- Emax*ypred/(ypred+EC50)
  ypred[ytype==2] <- ypred2[ytype==2]
  return(ypred)
}
psi0 <- c(0.9, 9, 9.5, 1.2)
pkpd.model<-saemixModel(model=pkpdmodel,modeltype=c("structural","structural"),
                        description="PK/PD model with direct response model", 
                        psi0=matrix(c(psi0,0.75,1,0,0),ncol=4, byrow=TRUE,
                                    dimnames=list(NULL, c("CL","V","Emax","EC50"))),transform.par=c(1,1,1,1),
                        error.model=c("combined","combined"))
saemix.options<-list(seed=1234,save=FALSE,save.graphs=FALSE)
pkpd.fit <- saemix(pkpd.model, pkpd.data, saemix.options)

###################################
# Bootstrap for multiple responses (Gaussian)
casepkpd <- saemix.bootstrap(pkpd.fit, method="case", nboot=nboot) 
write.table(casepkpd, file.path(saemixDir, "bootstrap", "pkpd_caseBootstrap.res"), quote=F, row.names=FALSE)

pkpd.fit@options$warnings<-TRUE
condpkpd <- saemix.bootstrap(pkpd.fit, method="conditional", nboot=nboot) 
write.table(condpkpd, file.path(saemixDir, "bootstrap", "pkpd_condBootstrap.res"), quote=F, row.names=FALSE)

###################################
if(FALSE) {
  casepkpd <- read.table(file.path(saemixDir, "bootstrap", "pkpd_caseBootstrap.res"), header=T)
  condpkpd<- read.table(file.path(saemixDir, "bootstrap", "pkpd_condBootstrap.res"), header=T)
  
  dat1<-cbind(pkpd.fit@results@conf.int[1:12,1:3],case.est=apply(casepkpd[,-c(1)],2,mean), case.sd=apply(casepkpd[,-c(1)],2,sd), 
        cond.est=apply(condpkpd[,-c(1)],2,mean), cond.sd=apply(condpkpd[,-c(1)],2,sd))
}


