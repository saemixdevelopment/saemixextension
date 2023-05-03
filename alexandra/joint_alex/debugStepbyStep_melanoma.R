#################################################
saemixDir <- "C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension"
workDir <- file.path(saemixDir, "alexandra","joint_alex")
setwd(workDir)

library(ggplot2)
library(Cairo)
library("viridis")  
library(rlang)
library(gsl)

# Chargement des fonctions originelles de la librairie
progDir<-file.path(saemixDir, "R")
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"func_plots.R")) # for saemix.plot.setoptions

################################################# Data and model (original files)
# Creating data and model objects

# data
data_joint <- read.csv("C:/Users/AlexandraLAVALLEY/Documents/Code_saemix/essais/melanoma_CMPR.csv", header=TRUE)
dataJM<-saemixData(name.data=data_joint, name.group=c("id"), name.predictors=c("time","obs"), name.covariates = c("age", "sex"),
                   name.response="obs")

#model JM lin + 2 competing risks
JMmodel<-function(psi,id,xidep) {

  h01 <- psi[id,1] 
  h02 <- psi[id,2]
  hcov1_age <- psi[id,3]
  hcov1_sex <- psi[id,4]
  hcov2_age <- psi[id,5]
  hcov2_sex <- psi[id,6]
  
  T<-xidep[,1] # vector of times partie survie ev 1
  ev = xidep$obs
  Nj <- length(T)
  cens<-which(ev==0)  # indices of censored observations
  ind1 <- which(ev==1) # indices of event 1
  ind2 <- which(ev==2) # indices of event 2
  
  
  haz1 = h01*exp(hcov1_age+hcov1_sex)
  H1 = h01*exp(hcov1_age+hcov1_sex)*T
  #F1 = function(t) p*(1-exp(-t))
  
  haz2 = h02*exp(hcov2_age+hcov2_sex)
  H2 = h02*exp(hcov2_age+hcov2_sex)*T
  #haz2 = rb*(1-F1(1000))*exp(-rb*T)/(1-(1-F1(1000))*(1-exp(-rb*T)))
  #H2 = -log(1-(1-F1(1000))*(1-exp(-rb*T)))
  
  logpdf <- rep(0,Nj)
  #logpdf[cens] <- log(exp(-H1[cens])+exp(-H2[cens])-1)
  logpdf[cens] <- -H1[cens]-H2[cens]
  logpdf[ind1] <- -H1[ind1] + log(haz1[ind1]) 
  logpdf[ind2] <- -H2[ind2] + log(haz2[ind2]) 
  
  
  ypred = logpdf
  
  return(ypred)
}


# joint TTE  
param<-c(0.00001,0.00001,0,0,0,0,0) 
jointTTE<-saemixModel(model=JMmodel,description="JM with CR",modeltype="likelihood",
                      psi0=matrix(param,ncol=7,byrow=TRUE,dimnames=list(NULL, c( "h01", "h02", "hcov1_age", "hcov1_sex", "hcov2_age", "hcov2_sex","dummy"))),
                      transform.par=c(1,1,0,0,0,0,0), covariance.model=diag(c(0,0,0,0,0,0,1)),covariate.model = t(matrix(c(0,0,0,0,1,0,0,1,1,0,0,1,0,0),ncol=2,byrow=T)),
                      fixed.estim = c(1,1,0,0,0,0,1), omega.init = diag(c(1,1,1,1,1,1,1)))



################################################# Initialisation
# Computational function
source(file.path(workDir,"multi_aux2.R"))
source(file.path(workDir,"multi_initializeMainAlgo.R"))

# Initialisation - debug
saemix.data<-dataJM
saemix.model<-jointTTE
saemix.options<-saemixControl(seed=12345)

# Initialisation - function
saemixObject<-new(Class="SaemixObject",data=saemix.data, model=saemix.model,options=saemix.options)
saemix.options<-saemixObject["options"]
saemix.model<-saemixObject["model"]
saemix.data<-saemixObject["data"]
saemix.data@ocov<-saemix.data@ocov[saemix.data@data[,"mdv"]==0,,drop=FALSE]
saemix.data@data<-saemix.data@data[saemix.data@data[,"mdv"]==0,]
saemix.data@ntot.obs<-dim(saemix.data@data)[1]
#  showall(saemixObject)

# Initialising random generator
set.seed(saemix.options$seed)

# Initialise lists
xinit<-initialiseMainAlgo.multi(saemix.data,saemix.model,saemix.options)
saemix.model<-xinit$saemix.model
Dargs<-xinit$Dargs
Uargs<-xinit$Uargs
varList<-xinit$varList
phiM<-xinit$phiM
mean.phi<-xinit$mean.phi
DYF<-xinit$DYF
opt<-xinit$opt
betas<-betas.ini<-xinit$betas
fixed.psi<-xinit$fixedpsi.ini
var.eta<-varList$diag.omega

# Initialise estimation results
if(length(grep("structural",Dargs$modeltype))>0){
  theta0<-c(fixed.psi,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res])
  parpop<-matrix(data=0,nrow=(saemix.options$nbiter.tot+1),ncol=(Uargs$nb.parameters+length(Uargs$i1.omega2)+length(saemix.model["indx.res"])))
  colnames(parpop)<-c(saemix.model["name.modpar"], saemix.model["name.random"], saemix.model["name.sigma"][saemix.model["indx.res"]])
  allpar<-matrix(data=0,nrow=(saemix.options$nbiter.tot+1), ncol=(Uargs$nb.betas+length(Uargs$i1.omega2)+length(saemix.model["indx.res"])))
  colnames(allpar)<-c(saemix.model["name.fixed"],saemix.model["name.random"], saemix.model["name.sigma"][saemix.model["indx.res"]])
} else{
  theta0<-c(fixed.psi,var.eta[Uargs$i1.omega2])
  parpop<-matrix(data=0,nrow=(saemix.options$nbiter.tot+1),ncol=(Uargs$nb.parameters+length(Uargs$i1.omega2)))
  colnames(parpop)<-c(saemix.model["name.modpar"], saemix.model["name.random"])
  allpar<-matrix(data=0,nrow=(saemix.options$nbiter.tot+1), ncol=(Uargs$nb.betas+length(Uargs$i1.omega2)))
  colnames(allpar)<-c(saemix.model["name.fixed"],saemix.model["name.random"])
}

parpop[1,]<-theta0
allpar[1,]<-xinit$allpar0

# using several Markov chains - only useful if passed back to main routine...
# 	chdat<-new(Class="SaemixRepData",data=saemix.data, nb.chains=saemix.options$nb.chains)
# 	NM<-chdat["NM"]
# 	IdM<-chdat["dataM"]$IdM
# 	yM<-chdat["dataM"]$yM
# 	XM<-chdat["dataM"][,saemix.data["name.predictors"],drop=FALSE]

# List of sufficient statistics - change during call to stochasticApprox
nstruct<-length(grep("structural",saemix.model@modeltype))
suffStat<-vector(mode="list", length=3+nstruct)
names(suffStat)[1:3]<-c("statphi1","statphi2","statphi3")
for(i in 1:length(suffStat)) suffStat[[i]]<-0
phi<-array(data=0,dim=c(Dargs$N, Uargs$nb.parameters, saemix.options$nb.chains))
structural.model<-saemix.model["model"]

################################################# E-step
source(file.path(workDir,"multi_estep.R"))

# E-step
kiter <- 1
xmcmc<-estep.multi(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
varList<-xmcmc$varList
DYF<-xmcmc$DYF
phiM<-xmcmc$phiM

################################################# M-step
source(file.path(workDir,"multi_mstep.R"))

# M-step
xstoch<-mstep.multi(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat)
dim(xstoch$suffStat[[1]])
head(xstoch$suffStat[[1]])
xstoch$suffStat[[2]]

################################################# Burn-in (kiter=1 to 5)
for (kiter in 1:5) {
  xmcmc<-estep.multi(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
  varList<-xmcmc$varList
  DYF<-xmcmc$DYF
  phiM<-xmcmc$phiM
  # no M-step as stepsize==0
  allpar[(kiter+1),]<-allpar[kiter,]
  
  # Actually useless, theta itself never used ?? => used here to print results
  if(length(grep("structural",Dargs$modeltype))>0)
    theta<-c(fixed.psi,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res]) else
      theta<-c(fixed.psi,var.eta[Uargs$i1.omega2])
  print(theta)
}

################################################# SA iterations after burn-in (kiter=6 to 150)
for (kiter in 6:149) {
  xmcmc<-estep.multi(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
  varList<-xmcmc$varList
  DYF<-xmcmc$DYF
  phiM<-xmcmc$phiM
  
  if(opt$stepsize[kiter]>0) {
    ############# Stochastic Approximation
    xstoch<-mstep.multi(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat)
    varList<-xstoch$varList
    mean.phi<-xstoch$mean.phi
    phi<-xstoch$phi
    betas<-xstoch$betas
    suffStat<-xstoch$suffStat
    
    beta.I<-betas[Uargs$indx.betaI]
    fixed.psi<-transphi(matrix(beta.I,nrow=1),saemix.model["transform.par"])
    betaC<-betas[Uargs$indx.betaC]
    var.eta<-mydiag(varList$omega)
    l1<-betas.ini
    l1[Uargs$indx.betaI]<-fixed.psi
    l1[Uargs$indx.betaC]<-betaC
    
    if(length(grep("structural",Dargs$modeltype))>0)
      allpar[(kiter+1),]<-c(l1,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res]) else
        allpar[(kiter+1),]<-c(l1,var.eta[Uargs$i1.omega2])
    
  } else  #end of loop on if(opt$stepsize[kiter]>0)
    allpar[(kiter+1),]<-allpar[kiter,]
  print(allpar[(kiter+1),])
}

################################################# Stopping SA (kiter=150 to 300)
kiter=150
if(opt$flag.fmin && kiter==saemix.options$nbiter.sa) {
  Uargs$COV1<-Uargs$COV[,Uargs$ind.fix11]
  ind.prov<-!(varList$ind.eta %in% Uargs$i0.omega2)
  varList$domega2<-varList$domega2[ind.prov,ind.prov,drop=FALSE] # keep in domega2 only indices of parameters with IIV
  varList$ind0.eta<-Uargs$i0.omega2
  varList$ind.eta<-1:(Uargs$nb.parameters)  	
  if(length(varList$ind0.eta)>0) varList$ind.eta<-varList$ind.eta[!(varList$ind.eta %in% varList$ind0.eta)] # update ind.eta, now only parameters with IIV
  Uargs$nb.etas<-length(varList$ind.eta)
  suffStat$statphi1<-0
  suffStat$statphi2<-0
  suffStat$statphi3<-0
}

for (kiter in 150:300) {
  xmcmc<-estep.multi(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
  varList<-xmcmc$varList
  DYF<-xmcmc$DYF
  phiM<-xmcmc$phiM
  
  if(opt$stepsize[kiter]>0) {
    ############# Stochastic Approximation
    xstoch<-mstep.multi(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat)
    varList<-xstoch$varList
    mean.phi<-xstoch$mean.phi
    phi<-xstoch$phi
    betas<-xstoch$betas
    suffStat<-xstoch$suffStat
    
    beta.I<-betas[Uargs$indx.betaI]
    fixed.psi<-transphi(matrix(beta.I,nrow=1),saemix.model["transform.par"])
    betaC<-betas[Uargs$indx.betaC]
    var.eta<-mydiag(varList$omega)
    l1<-betas.ini
    l1[Uargs$indx.betaI]<-fixed.psi
    l1[Uargs$indx.betaC]<-betaC
    
    if(length(grep("structural",Dargs$modeltype))>0)
      allpar[(kiter+1),]<-c(l1,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res]) else
        allpar[(kiter+1),]<-c(l1,var.eta[Uargs$i1.omega2])
    
  } else  #end of loop on if(opt$stepsize[kiter]>0)
    allpar[(kiter+1),]<-allpar[kiter,]
  print(allpar[(kiter+1),])
}
################################################# Smoothing phase (kiter=301 to 400)
for (kiter in 301:saemix.options$nbiter.tot) {
  xmcmc<-estep.multi(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
  varList<-xmcmc$varList
  DYF<-xmcmc$DYF
  phiM<-xmcmc$phiM
  
  if(opt$stepsize[kiter]>0) {
    ############# Stochastic Approximation
    xstoch<-mstep.multi(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat)
    varList<-xstoch$varList
    mean.phi<-xstoch$mean.phi
    phi<-xstoch$phi
    betas<-xstoch$betas
    suffStat<-xstoch$suffStat
    
    beta.I<-betas[Uargs$indx.betaI]
    fixed.psi<-transphi(matrix(beta.I,nrow=1),saemix.model["transform.par"])
    betaC<-betas[Uargs$indx.betaC]
    var.eta<-mydiag(varList$omega)
    l1<-betas.ini
    l1[Uargs$indx.betaI]<-fixed.psi
    l1[Uargs$indx.betaC]<-betaC
    
    if(length(grep("structural",Dargs$modeltype))>0)
      allpar[(kiter+1),]<-c(l1,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res]) else
        allpar[(kiter+1),]<-c(l1,var.eta[Uargs$i1.omega2])
    
  } else  #end of loop on if(opt$stepsize[kiter]>0)
    allpar[(kiter+1),]<-allpar[kiter,]
  print(allpar[(kiter+1),])
}

#################################################
source(file.path(workDir,"multi_estep.R"))
source(file.path(workDir,"multi_mstep.R"))
source(file.path(workDir,"multi_main.R"))
saemix.data<-dataJM
saemix.model<-jointTTE
saemix.options<-saemixControl(seed=12345, map=FALSE, fim=FALSE, ll.is=FALSE)

yfit <- saemix.multi(saemix.model, saemix.data, saemix.options)

#################################################
source(file.path(progDir,"func_FIM.R"))
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
# Replace main with the beginning of the algorithm
source("copy_estep_print.R")


