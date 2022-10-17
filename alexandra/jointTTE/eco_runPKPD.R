#################################################

saemixDir <- "/home/eco/work/saemix/saemixextension"
workDir <- file.path(saemixDir, "alexandra","jointTTE")
setwd(workDir)

# library(saemix)
library(ggplot2)
library("viridis")  

simPKPD <- FALSE

# Chargement des fonctions originelles de la librairie
progDir<-file.path(saemixDir, "R")
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixModel.R"))
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
source(file.path(progDir,"func_estimParam.R"))
source(file.path(progDir,"func_npde.R"))
source(file.path(progDir,"backward.R"))
source(file.path(progDir,"forward.R"))
source(file.path(progDir,"stepwise.R"))
source(file.path(progDir,"func_stepwise.R"))
source(file.path(progDir,"func_compare.R"))
source(file.path(progDir,"func_bootstrap.R"))
source(file.path(progDir,"func_exploreData.R"))
source(file.path(progDir,"func_discreteVPC.R"))

#################################################
############# SIMULATION DES DONNEES ############
################### PKPD ########################
#################################################

N<-200
# 4 doses possibles 
dose = sample(c(50,75,100,125),N,replace = T)

# ka_pop = 2
# V_pop = 8
# Cl_pop = 2 # 0.15 trop proche de flip-flop
# Emax_pop = 100
# EC50_pop = 5

param<-c(2,8,2,100,5)
omega<-c(1.4, 0.15, 0.3, 0.15, 0.2)
# omega_V = 0.15
# omega_Cl = 0.3

# omega_Emax = 0.15   # mettre omega_Emax et omega_EC50 = 0 pour que ça marche 
# omega_EC50 = 0.2

sigma_prop = 0.2
sigma_addpk = 0.5
sigma_addpd = 5

iparam<-param
if(simPKPD) {
  pk = function(x, iparam) (D/iparam[2])*(iparam[1]/(iparam[1]-iparam[3]/iparam[2]))*(exp(-(iparam[3]/iparam[2])*x)-exp(-iparam[1]*x))
  pd = function(x, iparam) iparam[4]*pk(x, iparam)/(pk(x, iparam)+iparam[5])
  data_pkpd.add = data.frame(id=rep(1:N,each=20),time = rep(1:10,2*N),obs=NA,ytype=rep(c(rep(1,10),rep(2,10)),N),dose=NA, ka=NA, V=NA, CL=NA, Emax=NA, EC50=NA)
  data_pkpd.prop = data.frame(id=rep(1:N,each=20),time = rep(1:10,2*N),obs=NA,ytype=rep(c(rep(1,10),rep(2,10)),N),dose=NA, ka=NA, V=NA, CL=NA, Emax=NA, EC50=NA)
  for (i in 1:N){
    for(ipar in 1:5) {
      iparam[ipar]<-param[ipar]*exp(rnorm(1,0,omega[ipar]))
      data_pkpd.add[data_pkpd.add$id==i,(ipar+5)] = iparam[ipar]
      data_pkpd.prop[data_pkpd.prop$id==i,(ipar+5)] = iparam[ipar]
    }
    # ka = ka_pop*exp(rnorm(1,0,omega_ka))
    # V = V_pop*exp(rnorm(1,0,omega_V))
    # Cl = Cl_pop*exp(rnorm(1,0,omega_Cl))
    # Emax = Emax_pop*exp(rnorm(1,0,omega_Emax))
    # EC50 = EC50_pop*exp(rnorm(1,0,omega_EC50))
    
    D = dose[i]
    # Additive model
    data_pkpd.add$obs[data_pkpd.add$id==i & data_pkpd.add$ytype==1] = pk(1:10, iparam)+rnorm(10,0,sigma_addpk)
    data_pkpd.add$obs[data_pkpd.add$id==i & data_pkpd.add$ytype==2] = pd(1:10, iparam)+rnorm(10,0,sigma_addpd)
    data_pkpd.add$dose[data_pkpd.add$id==i] = D
    # Proportional model
    data_pkpd.prop$obs[data_pkpd.prop$id==i & data_pkpd.prop$ytype==1] = pk(1:10, iparam)*(1+rnorm(10,0,sigma_prop))
    data_pkpd.prop$obs[data_pkpd.prop$id==i & data_pkpd.prop$ytype==2] = pd(1:10, iparam)*(1+rnorm(10,0,sigma_prop))
    data_pkpd.prop$dose[data_pkpd.prop$id==i] = D
  }
  data_pkpd.add <- data_pkpd.add[order(data_pkpd.add$id, data_pkpd.add$time, data_pkpd.add$ytype),]
  data_pkpd.prop <- data_pkpd.prop[order(data_pkpd.prop$id, data_pkpd.prop$time, data_pkpd.prop$ytype),]
  write.csv(data_pkpd.add, file.path(workDir, "simAlex_pkpd_additive.csv"), quote=F, row.names = FALSE)
  write.csv(data_pkpd.prop, file.path(workDir, "simAlex_pkpd_proportional.csv"), quote=F, row.names = FALSE)
  # write.csv(data_pkpd.add[data_pkpd.add$ytype==1,], file.path(workDir, "simAlex_pk_additive.csv"), quote=F, row.names = FALSE)
  # write.csv(data_pkpd.add[data_pkpd.add$ytype==2,], file.path(workDir, "simAlex_pd_additive.csv"), quote=F, row.names = FALSE)
  # write.csv(data_pkpd.prop[data_pkpd.prop$ytype==1,], file.path(workDir, "simAlex_pk_proportional.csv"), quote=F, row.names = FALSE)
  # write.csv(data_pkpd.prop[data_pkpd.prop$ytype==2,], file.path(workDir, "simAlex_pd_proportional.csv"), quote=F, row.names = FALSE)
} else {
  data_pkpd.add <- read.csv(file.path(workDir, "simAlex_pkpd_additive.csv"), header=TRUE)
  data_pkpd.prop <- read.csv(file.path(workDir, "simAlex_pkpd_proportional.csv"), header=TRUE)
}

data_pkpd <- data_pkpd.add
ggplot(data_pkpd, aes(x=time, y=obs, group=id, colour=as.factor(dose))) + geom_line() + theme_minimal() + scale_color_viridis(discrete=TRUE) + facet_wrap(.~ytype, scales="free") + ylab("Outcome") +  guides(colour=guide_legend(title="Dose"))

data_pkpd <- data_pkpd.prop
ggplot(data_pkpd, aes(x=time, y=obs, group=id, colour=as.factor(dose))) + geom_line() + theme_minimal() + scale_color_viridis(discrete=TRUE) + facet_wrap(.~ytype, scales="free") + ylab("Outcome") +  guides(colour=guide_legend(title="Dose"))

################### Fitting only PK ########################
source("copy_estep_print.R")

# PK model
modelPK<-function(psi,id,xidep) {
  ka <- psi[id,1] 
  V <- psi[id,2]
  Cl <- psi[id,3]
  tim<-xidep[,1] 
  D <- xidep[,2]
  ypred <- (D/V)*(ka/(ka-Cl/V))*(exp(-(Cl/V)*tim)-exp(-ka*tim))
  return(ypred)
}

data_pkpd <- data_pkpd.add[data_pkpd.add$ytype==1,]
saemix.data.add<-saemixData(name.data=data_pkpd, name.group=c("id"), name.predictors=c("time","dose"), name.response="obs")
data_pkpd <- data_pkpd.prop[data_pkpd.prop$ytype==1,]
saemix.data.prop<-saemixData(name.data=data_pkpd, name.group=c("id"), name.predictors=c("time","dose"), name.response="obs")

pkmodel.add<-saemixModel(model=modelPK,description="Only PK",modeltype="structural",
                           psi0=matrix(c(2, 8, 0.15),ncol=3,byrow=TRUE,dimnames=list(NULL, c("ka","V","Cl"))),
                           transform.par=c(1,1,1), covariance.model=diag(c(1,1,1)),
                           omega.init = diag(c(1.96,0.02,0.09)),error.model = c("constant"),error.init = c(1,0))
pkmodel.prop<-saemixModel(model=modelPK,description="Only PK",modeltype="structural",
                         psi0=matrix(c(2, 8, 0.15),ncol=3,byrow=TRUE,dimnames=list(NULL, c("ka","V","Cl"))),
                         transform.par=c(1,1,1), covariance.model=diag(c(1,1,1)),
                         omega.init = diag(c(1,1,1)),error.model = c("proportional"),error.init = c(0,1))
#                         omega.init = diag(c(1.96,0.02,0.09)),error.model = c("proportional"),error.init = c(0,1))
saemix.options<-list(seed=12345,save=FALSE,save.graphs=FALSE, fim=FALSE, displayProgress=FALSE, nbiter.saemix=c(50,20))
pkfit.add <- saemix(pkmodel.add,saemix.data.add,saemix.options)
pkfit.prop <- saemix(pkmodel.prop,saemix.data.prop,saemix.options)

################### Fitting only PD ########################
# PD with individual parameters for PK
modelPD<-function(psi,id,xidep) {
  Emax <- psi[id,1]
  EC50 <- psi[id,2]
  tim<-xidep[,1] 
  D <- xidep[,2]
  ka <- xidep[,3]
  V <- xidep[,4]
  Cl <- xidep[,5]
  ypred <- (D/V)*(ka/(ka-Cl/V))*(exp(-(Cl/V)*tim)-exp(-ka*tim))
  ypred <- Emax*ypred/(ypred+EC50)
  return(ypred)
}

data_pd <- data_pkpd.add[data_pkpd.add$ytype==2,]
saemix.data.add<-saemixData(name.data=data_pd, name.group=c("id"), name.predictors=c("time","dose", "ka","V","CL"), name.response="obs")
data_pd <- data_pkpd.prop[data_pkpd.prop$ytype==2,]
saemix.data.prop<-saemixData(name.data=data_pd, name.group=c("id"), name.predictors=c("time","dose", "ka","V","CL"), name.response="obs")

pdmodel.add<-saemixModel(model=modelPD,description="Only PD",modeltype="structural",
                         psi0=matrix(c(100, 5),ncol=2,byrow=TRUE,dimnames=list(NULL, c("Emax","EC50"))),
                         transform.par=c(1,1), covariance.model=diag(c(1,1)),
                         omega.init = diag(c(0.02,0.04)),error.model = c("constant"),error.init = c(1,0))
pdmodel.prop<-saemixModel(model=modelPD,description="Only PD",modeltype="structural",
                          psi0=matrix(c(100, 5),ncol=2,byrow=TRUE,dimnames=list(NULL, c("Emax","EC50"))),
                          transform.par=c(1,1), covariance.model=diag(c(1,1)),
                          omega.init = diag(c(0.02,0.04)),error.model = c("proportional"),error.init = c(0,1))
saemix.options<-list(seed=12345,save=FALSE,save.graphs=FALSE, fim=FALSE, displayProgress=FALSE, nbiter.saemix=c(50,20))
pdfit.add <- saemix(pdmodel.add,saemix.data.add,saemix.options)
pdfit.prop <- saemix(pdmodel.prop,saemix.data.prop,saemix.options)

################### Fitting PK/PD with Alexandra's code ########################

## LANCER LES FONCTIONS MODIFIEES SUIVANTES 
# modified functions: saemix, initialiseMainAlgo, estep, mstep, compute.LLy, compute.Uy => added .alex
# all other functions are taken from saemix
source("functions_alex.R")
source("main_2longi.R")

# Refait sans la fonction pk à l'intérieur
LONGIjointmodel<-function(psi,id,xidep) {
  ytype<-xidep$ytype
  
  ka <- psi[id,1] 
  V <- psi[id,2]
  Cl <- psi[id,3]
  Emax <- psi[id,4]
  EC50 <- psi[id,5]
  
  tim<-xidep[,1] 
  D <- xidep[,2]
  
  ypred <- (D/V)*(ka/(ka-Cl/V))*(exp(-(Cl/V)*tim)-exp(-ka*tim))
  ypred[ytype==2] <- Emax*ypred[ytype==2]/(ypred[ytype==2]+EC50)
  
  return(ypred)
}

# Proportional error model
data_pkpd <- data_pkpd.prop
saemix.dataPKPD.prop<-saemixData(name.data=data_pkpd, name.group=c("id"), name.predictors=c("time","dose"), name.response="obs",name.ytype = "ytype")
data_pkpd <- data_pkpd.add
saemix.dataPKPD.add<-saemixData(name.data=data_pkpd, name.group=c("id"), name.predictors=c("time","dose"), name.response="obs",name.ytype = "ytype")

pkpdmodel.add<-saemixModel(model=LONGIjointmodel,description="joint pkpd",modeltype="structural",
                          psi0=matrix(param,ncol=5,byrow=TRUE,dimnames=list(NULL, c("ka","V","Cl","Emax","EC50"))),
                          transform.par=c(1,1,1,1,1), covariance.model=diag(c(1,1,1,1,1)),fixed.estim = c(1,1,1,1,1),
                          omega.init = diag(rep(0.5,5)),error.model = c("constant","constant"),error.init = c(1,0,5,0),
                          name.sigma = c("a.1","b.1","a.2","b.2"))

pkpdmodel.prop<-saemixModel(model=LONGIjointmodel,description="joint pkpd",modeltype="structural",
                          psi0=matrix(param,ncol=5,byrow=TRUE,dimnames=list(NULL, c("ka","V","Cl","Emax","EC50"))),
                          transform.par=c(1,1,1,1,1), covariance.model=diag(c(1,1,1,1,1)),fixed.estim = c(1,1,1,1,1),
#                          omega.init = diag(c(1.96,0.02,0.09,0.02,0.04)),error.model = c("proportional","proportional"),error.init = c(0,1,0,1),
                          omega.init = diag(rep(0.5,5)),error.model = c("proportional","proportional"),error.init = c(0,1,0,1),
                          name.sigma = c("a.1","b.1","a.2","b.2"))

saemix.options<-list(seed=12345,save=FALSE,save.graphs=FALSE, fim=FALSE, map=FALSE, displayProgress=FALSE)
saemix.options<-list(seed=12345,save=FALSE,save.graphs=FALSE, fim=FALSE, map=FALSE, displayProgress=FALSE, nbiter.saemix=c(50,20))

# Creating a temporary file
namUyfile<-file.path(workDir,"resUy.res")
system(paste("touch",namUyfile))

pkpdfit.add <- saemix.alex(pkpdmodel.add,saemix.dataPKPD.add,saemix.options)

namUyfile<-file.path(workDir,"resUy.res")
system(paste("touch",namUyfile))
pkpdfit.prop <- saemix.alex(pkpdmodel.prop,saemix.dataPKPD.prop,saemix.options)

head(pkpdfit.prop@results@cond.mean.phi)


# debug
model <- pkpdmodel.prop
data <- saemix.dataPKPD.prop
control <- saemix.options

# Running kiter by steps 
## for estep.alex.outcome - same names for Uargs, Dargs, opt, mean.phi, varList, DYF, phiM
# xmcmc<-estep.alex.outcome(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)

## for compute.LLy.outcome
phiM <- phiM
args<-Uargs
pres<-varList$pres


# Can Alex's code fit single responses ? No (n_mark controls the nb of responses and is set)
# data_pd <- data_pkpd.prop[data_pkpd.prop$ytype==2,]
# saemix.data.prop<-saemixData(name.data=data_pd, name.group=c("id"), name.predictors=c("time","dose", "ka","V","CL"), name.response="obs")
# pdmodel.prop<-saemixModel(model=modelPD,description="Only PD",modeltype="structural",
#                           psi0=matrix(c(100, 5),ncol=2,byrow=TRUE,dimnames=list(NULL, c("Emax","EC50"))),
#                           transform.par=c(1,1), covariance.model=diag(c(1,1)),
#                           omega.init = diag(c(0.02,0.04)),error.model = c("proportional"),error.init = c(0,1))
# alexpdfit.prop <- saemix.alex(pdmodel.prop,saemix.data.prop,saemix.options)

########################################################
# Debug

yfit <- pkpdfit.prop
phiM <- yfit$phiM
Uargs<-args<-yfit$Uargs
Dargs <- yfit$Dargs
DYF<-yfit$DYF
varList<-yfit$varList
pres<-yfit$varList$pres
opt<-yfit$opt
suffStat <- yfit$suffStat
phi<-yfit$phi
betas<-yfit$betas

structural.model <-pkpdmodel.prop@model

############################
## Calcul des vraisemblances
Uy.50 <- compute.LLy.alex(yfit$phiM, yfit$Uargs, yfit$Dargs, yfit$DYF, yfit$varList$pres)
Uy.mean50 <- compute.LLy.alex(yfit$mean.phi, yfit$Uargs, yfit$Dargs, yfit$DYF, yfit$varList$pres)
head(exp(yfit$phiM))
exp(yfit$mean.phi[1,])


compute.LLy.split <- function (phiM, args, Dargs, DYF, pres) {
  psiM <- transphi(phiM, Dargs$transform.par)
  fpred <- Dargs$structural.model(psiM, Dargs$IdM, Dargs$XM)
  ytype = Dargs[["XM"]][["ytype"]]
  for (ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype == ityp] <- log(cutoff(fpred[Dargs$XM$ytype == ityp]))
  gpred <- error(fpred, pres, Dargs$XM$ytype)
  DYF[args$ind.ioM[ytype==1]] <- 0.5 * ((Dargs$yM[ytype==1] - fpred[ytype==1])/gpred[ytype==1])^2 + log(gpred[ytype==1])
  DYF[args$ind.ioM[ytype==2]] <- 0.5 * ((Dargs$yM[ytype==2] - fpred[ytype==2])/gpred[ytype==2])^2 + log(gpred[ytype==2])

  #DYF[Uargs$ind.ioM] <- 0.5 * ((Dargs$yM - fpred)/gpred)^2 + log(gpred)
  # ca doit faire la même chose ça  # oui
  # DYF1<-DYF
  # DYF1[Uargs$ind.ioM] <- 0.5 * ((Dargs$yM - fpred)/gpred)^2 + log(gpred)
  # summary(c(DYF-DYF1))
  
  # Separating both Uy
  DYF.rep1<-DYF[args$ind.ioM[ytype==1]]
  id.rep1<-Dargs$IdM[args$ind.ioM[ytype==1]]
  DYF.rep2<-DYF[args$ind.ioM[ytype==2]]
  id.rep2<-Dargs$IdM[args$ind.ioM[ytype==2]]
  Uy.rep1 <- tapply(DYF.rep1, id.rep1, sum)
  Uy.rep2 <- tapply(DYF.rep2, id.rep2, sum)
  
  U <- colSums(DYF)
  return(list(Uy=list(rep1=Uy.rep1, rep2=Uy.rep2), U=U))
}

Uy.50.split <- compute.LLy.split(yfit$phiM, yfit$Uargs, yfit$Dargs, yfit$DYF, yfit$varList$pres)
Uy.mean50.split <- compute.LLy.split(yfit$mean.phi, yfit$Uargs, yfit$Dargs, yfit$DYF, yfit$varList$pres)

# Should be ~0
summary(Uy.50.split$Uy$rep1+Uy.50.split$Uy$rep2-Uy.50.split$U)

head(Uy.50.split$Uy$rep1)
head(Uy.50.split$Uy$rep2)
head(Uy.mean50.split$Uy$rep1)
head(Uy.mean50.split$Uy$rep2)
head(Uy.mean50.split$Uy$rep1 - Uy.50.split$Uy$rep1)
head(Uy.mean50.split$Uy$rep2 - Uy.50.split$Uy$rep2)

############################
## Calcul des vraisemblances séparées pour chacune des réponses

## Lecture fichier temporaire
res<-scan(namUyfile)
res<-matrix(res, byrow=T, ncol=202)
colnames(res)<-c('kiter',"ytype",1:200)
# idx.meanphi and idx.meanphi+1: value of U.y at the beginning of the E-step
# next series (16 or 32 lines in all) are values obtained during the 3 kernels
idx.meanphi.rep1 <- which(!duplicated(res[,1]))
idx.meanphi.rep2<-idx.meanphi.rep1+1

# kernel 1: 2 iterations [4 lines]
## delta=delta(LL) between new samples and mean.phi
# kernel 2: 2*nb.etas iterations (here, 10) [20 lines]
## delta=delta(LL) between new samples and samples after kernel 1 + delta(1/2 eta.Omega.eta) between new etas and etas after kernel 1
# kernel 3: 8 or 16 iterations here (depends on kiter %% something) [16 or 32 lines]
## same computation as for kernel 2

## Computing deltaU
idx.meanphi.rep1 <- c(idx.meanphi.rep1, dim(res)[1]+1, dim(res)[1]+2)
nsuj<-dim(res)[2]-2
deltaU <- NULL
for(irep in 1:(length(idx.meanphi.rep1)-2)) {
  for(iout in 1:2) {
    # initial value from mean.phi
    lbase <- res[idx.meanphi.rep1[irep]+(iout-1),]
    # kernel 1
    idx1<-idx.meanphi.rep1[irep]+(iout-1)+4
    l1<-res[idx1,]
    tab<-data.frame(kiter=irep, kernel=1, id=1:nsuj, outcome=ifelse(iout==1,"PK","PD"), delta=l1[-c(1:2)]-lbase[-c(1:2)])
    deltaU<-rbind(deltaU, tab)
    # kernel 2
    idx1<-idx.meanphi.rep1[irep]+(iout-1)+24
    l1<-res[idx1,]
    tab<-data.frame(kiter=irep, kernel=2, id=1:nsuj, outcome=ifelse(iout==1,"PK","PD"), delta=l1[-c(1:2)]-lbase[-c(1:2)])
    deltaU<-rbind(deltaU, tab)
    # kernel 3
    idx1<-idx.meanphi.rep1[irep+1]+(iout-1)-2
    l1<-res[idx1,]
    tab<-data.frame(kiter=irep, kernel=3, id=1:nsuj, outcome=ifelse(iout==1,"PK","PD"), delta=l1[-c(1:2)]-lbase[-c(1:2)])
    deltaU<-rbind(deltaU, tab)
  }
}

summary(deltaU)
# outliers
deltaU1<-deltaU[deltaU$delta<1000,]
ggplot(deltaU1, aes(x=kiter, y=delta, group=as.factor(kernel), colour=as.factor(kernel))) + geom_point() + geom_hline(yintercept=-log(0.5)) + geom_hline(yintercept = 7) + facet_wrap(.~outcome)

xll <- cut(deltaU$delta, breaks=c(min(deltaU$delta), 0,7,100, 1000, max(deltaU$delta)), include.lowest=TRUE)
# slightly larger amount of high delta values for PK (makes sense, bad PK parameters have an impact on both PK and PD while bad PD estimates only affect PD)
table(xll, deltaU$outcome)
# Higher values of delta associated with kernel 1 (makes sense, sample all eta's so bound to have a chance of hitting really bad values)
table(xll, deltaU$kernel)

# No glaring difference between delta's for PK and PD, maybe slightly more variability (SD) in PK and slightly higher values on average but nothing completely unbalanced
library(tidyverse)
deltaU %>%
  group_by(outcome) %>%
  summarise(mean(delta), sd(delta), min(delta), max(delta))
############################
# Calcul des stats suffisantes

summary(phiM)

############################
## Prédictions individuelles



############################
# Step by step
model=saemix.model
data=saemix.data
control=saemix.options
