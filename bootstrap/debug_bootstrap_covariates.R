###################################################
saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
nielsDir <- "/home/eco/work/theses/niels/bootstrap/bootstrap_niels"

# Libraries
library(ggplot2)
library(MASS)
library(rlang)
library(npde)

# Sourcing saemix functions
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"main.R"))
source(file.path(progDir,"func_aux.R"))
source(file.path(progDir,"main_initialiseMainAlgo.R"))
source(file.path(progDir,"main_estep.R"))
source(file.path(progDir,"main_mstep.R"))
source(file.path(progDir,"func_FIM.R"))
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
source(file.path(progDir,"func_bootstrap.R"))
source(file.path(progDir,"func_exploreData.R"))
source(file.path(progDir,"func_discreteVPC.R"))

###################################################
setwd(nielsDir)
# ARSACS fit
growth.log_param<-function(psi,id,xidep) {
  x<-xidep[,1]
  s_0<-psi[id,1]
  s_max<-psi[id,2]
  alpha <- psi[id,3]
  t50 <- psi[id,4]
  f <- s_max+(s_0-s_max)*1/(1+(exp(alpha*x)-1)/(exp(alpha*t50)-1))
  return(f)
}

i<-1
sara_temp <- read.csv(paste0("sara_mice_",i,".csv"))

# tidyverse
#med_INAS <- (sara_temp%>%filter(!duplicated(ID))%>%summarize(med=median(INAS,na.rm=T)))
#sara_temp$INAS <-  unlist(log(sara_temp$INAS)-log(med_INAS[,"med"]))
#med_BMI <- (sara_temp%>%filter(!duplicated(ID))%>%summarize(med=median(BMI,na.rm=T)))
#sara_temp$BMI <-  unlist(log(sara_temp$BMI)-log(med_BMI[,"med"]))
# covmodel <- matrix(0,3,4)%>%replace(c(12,8,5,1,10),1)

# same without tidyverse
saracov <- sara_temp[!duplicated(sara_temp$ID),]
sara_temp$INAS <- log(sara_temp$INAS)-log(median(saracov$INAS, na.rm=T))
sara_temp$BMI <- log(sara_temp$BMI)-log(median(saracov$BMI, na.rm=T))
covmodel <- rep(0,12)
covmodel[c(1,5,8,10,12)]<-1
covmodel<-matrix(covmodel,nrow=3)

saraARSACS_TSO <- saemixData(name.data=sara_temp, name.group = "ID", name.predictors = c("TSO"),name.response="DV"
                             , units=list(x="yr",y="-"),name.covariates=c("SEX","INAS","AOO_15_"))
saemix.model.log.p<-saemixModel(model=growth.log_param,description="Logistic model",
                                psi0=matrix(c(5,40,0.1,30),ncol=4,byrow=TRUE,dimnames=list(NULL,c("s_0","s_max","alpha","t50"))),
                                transform.par=c(1,1,1,1),fixed.estim=c(1,1,1,1),covariance.model=diag(c(1,0,1,1)),
                                covariate.model = covmodel,
                                omega.init=matrix(0.1,4,4),
                                error.model="constant",)
saemix.options<-list(print=T,nb.chains=10,seed=200182,directory="",nmc.is = 10000,
                     save.graphs=FALSE,displayProgress=FALSE,nbiter.saemix = c(600, 200))

fit.arsacs<-saemix(saemix.model.log.p,saraARSACS_TSO, control=saemix.options)

# Results from Niels
boot.temp <- read.csv(paste0("results_bootstrap_",i,".csv"))
summary(boot.temp)

##############################################
# Bootstraps without simulated annealing
fit.arsacs2 <- fit.arsacs
fit.arsacs2@options$nbiter.sa<-0

####################
# Case bootstrap
boot.arsacs.case <- saemix.bootstrap(fit.arsacs2,method="case",nboot=1)

# check à la main
set.seed(12345)
data.boot.case <- dataGen.case(fit.arsacs2)

####################
# Parametric Bootstrap
fit.arsacs2 <- fit.arsacs
fit.arsacs2@options$nbiter.sa<-0
boot.arsacs.param <- saemix.bootstrap(fit.arsacs2,method="parametric",nboot=1)

set.seed(12345)
data.boot.par <- dataGen.Par(fit.arsacs2)

####################
# Bootstrap with conditional distribution
boot.arsacs.cNP <- saemix.bootstrap(fit.arsacs2,method="conditional",nboot=1)

set.seed(12345)
data.boot.cNP <- dataGen.NP(fit.arsacs2, conditional=TRUE)

# Fails with message
# Error in d1.omega[Uargs$ind.fix11, ] * (t(Uargs$COV1) %*% (suffStat$statphi1 -  : 
#                                                             tableaux de tailles inadéquates
# ? because SA removed maybe ?

##############################################
# Debug parametric and conditional in the presence of covariates...
saemixObject <- fit.arsacs
nsamp<-100

# Main function
saemix.options<-saemixObject["options"]
saemix.options$directory<-"current"
saemix.options$fix.seed<-FALSE
saemix.options$map<-FALSE   # Only parameter estimates are required for bootstrap
saemix.options$fim<-FALSE
saemix.options$displayProgress<-FALSE 
saemix.options$save.graphs<-FALSE
saemix.options$ll.is<-FALSE
saemix.options$print<-FALSE

saemixObject<-saemix.predict(saemixObject)
saemixObject<-conddist.saemix(saemixObject, nsamp=nsamp)
eta.sampc<-centerDist.NPcond(saemixObject, nsamp=nsamp)

if(saemixObject@model@modeltype=="structural") idx.eps<-saemixObject@model@indx.res else idx.eps<-integer(0)
idx.iiv<-saemixObject@model@indx.omega
idx.rho<-which(saemixObject@model@covariance.model[lower.tri(saemixObject@model@covariance.model)]==1)
bootstrap.distribution<-failed.runs<-data.frame()
nelements <- length(saemixObject@results@fixed.effects)+length(idx.iiv)+length(idx.rho)+length(idx.eps)
model.boot<-saemixObject["model"]
model.boot@psi0 <- model.boot["betaest.model"]
model.boot@psi0[model.boot["betaest.model"]==1]<-saemixObject@results@fixed.effects

# Generate dataset
data.boot <- dataGen.NP(saemixObject, nsamp=nsamp,eta.sampc=eta.sampc, conditional=TRUE)

fit.boot<-try(saemix(model.boot, data.boot, saemix.options))

# in a loop - works :-/
for (iboot in 1:5) {
  data.boot <- dataGen.NP(saemixObject, nsamp=nsamp,eta.sampc=eta.sampc, conditional=TRUE)
  fit.boot<-try(saemix(model.boot, data.boot, saemix.options))
  res<-fit.boot@results
  l1<-c(iboot,res@fixed.effects, diag(res@omega)[idx.iiv])
  if(length(idx.rho)>0) l1<-c(l1,res@omega[lower.tri(res@omega)][idx.rho])
  if(length(idx.eps)>0) l1<-c(l1, res@respar[idx.eps])
  if(length(res@ll.lin)>0) l1<-c(l1, res@ll.lin)
  bootstrap.distribution<-rbind(bootstrap.distribution,l1) 
}

# Trying to isolate the issue
boot.cnp <- saemix.bootstrapcNP(saemixObject, nboot=10)


saemix.bootstrapcNP<-function(saemixObject, method="conditional", nboot=200, nsamp=100, saemix.options=NULL) {
    saemixObject<-saemix.predict(saemixObject) # estimate individual parameters and compute residuals (currently iwres are needed also for     ndone <- dim(saemixObject@results@phi.samp)
    ndone <- dim(saemixObject@results@phi.samp)
    if(!is.null(ndone)) ndone<-ndone[3] else ndone<-0
    if(ndone<nsamp) {
      if(saemixObject@options$warnings) message("Not enough samples in the object, sampling from the conditional distribution\n")
      saemixObject<-conddist.saemix(saemixObject, nsamp=nsamp) # estimate conditional distributions and sample residuals
    }
    eta.sampc<-centerDist.NPcond(saemixObject, nsamp=nsamp) # Center eta samples from the conditional distribution, to avoid doing this repeatedly
  if(is.null(saemix.options)) {
    #      saemix.options<-list(directory="current",fix.seed=FALSE,map=FALSE,ll.is=FALSE,displayProgress=FALSE,save.graphs=FALSE,print=FALSE)
    saemix.options<-saemixObject["options"]
    saemix.options$directory<-"current"
    saemix.options$fix.seed<-FALSE
    saemix.options$map<-FALSE   # Only parameter estimates are required for bootstrap
    saemix.options$fim<-FALSE
    saemix.options$displayProgress<-FALSE 
    saemix.options$save.graphs<-FALSE
    saemix.options$ll.is<-FALSE
    saemix.options$print<-FALSE
  }
  verbose <- saemix.options$warnings
  if(saemixObject@model@modeltype=="structural") idx.eps<-saemixObject@model@indx.res else idx.eps<-integer(0)
  idx.iiv<-saemixObject@model@indx.omega
  idx.rho<-which(saemixObject@model@covariance.model[lower.tri(saemixObject@model@covariance.model)]==1)
  bootstrap.distribution<-failed.runs<-data.frame()
  nelements <- length(saemixObject@results@fixed.effects)+length(idx.iiv)+length(idx.rho)+length(idx.eps)
  model.boot<-saemixObject["model"]
  model.boot@psi0 <- model.boot["betaest.model"]
  model.boot@psi0[model.boot["betaest.model"]==1]<-saemixObject@results@fixed.effects
  for(iboot in 1:nboot) {
    if(method=="case")  
      data.boot <- dataGen.case(saemixObject)
    if(method=="residual")
      data.boot <- dataGen.NP(saemixObject,conditional=FALSE)
    if(method=="conditional")
      data.boot <- dataGen.NP(saemixObject, nsamp=nsamp,eta.sampc=eta.sampc, conditional=TRUE)
    if(method=="parametric")
      data.boot <- dataGen.Par(saemixObject)
    fit.boot<-try(saemix(model.boot, data.boot, saemix.options))
    if(is(fit.boot,"try-error")) {
      l1<-c(iboot,rep(NA,nelements))
      failed.runs <- rbind(failed.runs, c(iboot, fit.boot))
    } else {
      res<-fit.boot@results
      l1<-c(iboot,res@fixed.effects, diag(res@omega)[idx.iiv])
      if(length(idx.rho)>0) l1<-c(l1,res@omega[lower.tri(res@omega)][idx.rho])
      if(length(idx.eps)>0) l1<-c(l1, res@respar[idx.eps])
      if(length(res@ll.lin)>0) l1<-c(l1, res@ll.lin)
      
    }
    #    l1<-c(iboot,res@fixed.effects, diag(res@omega)[idx.iiv],res@omega[lower.tri(res@omega)][idx.rho],res@respar[idx.eps], res@se.fixed, res@se.omega[idx.iiv],res@se.cov[lower.tri(res@se.cov)][idx.rho], res@se.respar[idx.eps],res@ll.lin)
    bootstrap.distribution<-rbind(bootstrap.distribution,l1) 
  }
  # Names
  nampar<-colnames(saemixObject@model@covariance.model)
  namcol<-c(saemixObject@results@name.fixed, saemixObject@results@name.random)
  if(length(idx.rho)>0) {
    for(i in 1:(length(nampar)-1)) {
      for(j in (i+1):length(nampar)) {
        if(saemixObject@model@covariance.model[i,j]==1) {
          namcol<-c(namcol,paste("cov.",nampar[i],nampar[j],sep=""))
        }
      }
    }
  }
  if(length(idx.eps)>0) namcol<-c(namcol,saemixObject@model@name.sigma[idx.eps])
  if(length(res@ll.lin)>0) namcol<-c(namcol,"LL.lin")
  #  namcol<-c(saemixObject@results@name.fixed,saemixObject@results@name.random,namcol, saemixObject@results@name.sigma[saemixObject@results@indx.res])
  #  colnames(bootstrap.distribution)<-c("Replicate",namcol,paste("SE",namcol,sep="."),"LL.lin")
  colnames(bootstrap.distribution)<-c("Replicate",namcol)
  if(verbose && dim(failed.runs)[1]>0) {
    cat(dim(failed.runs)[1],"failed:\n")
    print(head(failed.runs))
  }
  return(bootstrap.distribution)
}

