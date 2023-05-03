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
data<- read.csv("C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/lin_mixed.csv", header=TRUE)

dataJM<-saemixData(name.data=data, name.group=c("id"), name.predictors=c("time"),
                   name.response="obs")
# model
modlin<-function(psi,id,xidep) {
  b0 <- psi[id,1] 
  b1 <- psi[id,2]
  
  tim<-xidep[,1] 
  
  ypred <- b0+b1*tim
  return(ypred)
}

# Proportional error model
param<-c(15,0.3)
mod<-saemixModel(model=modlin,description="lmem",modeltype="structural",
                 psi0=matrix(param,ncol=2,byrow=TRUE,dimnames=list(NULL, c("b0","b1"))),
                 transform.par=c(0,0), covariance.model=diag(c(1,1)),fixed.estim = c(1,1),
                 omega.init = diag(c(0.5,0.5)),error.model = "constant")


################################################# Initialisation
# Computational function
source(file.path(workDir,"multi_aux2.R"))
source(file.path(workDir,"multi_initializeMainAlgo.R"))

# Initialisation - debug
saemix.data<-dataJM
saemix.model<-mod
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

### ici on calacule les SE methode 1 de la feuille: on dérive par omega**2

phibis = matrix(phi,nrow = 100, ncol=2)
pssi = transphi(phibis,Dargs$transform.par)
pred = Dargs$structural.model(suffStat$statphi1,Dargs$IdM,Dargs$XM)

mu0 = l1[1]
mu1 = l1[2]
omega0.2 = var.eta[Uargs$i1.omega2][1]
omega1.2 = var.eta[Uargs$i1.omega2][2]
sigma = varList$pres[Uargs$ind.res]

bigmat=matrix(data=0,nrow = 5,ncol = 5)

for (i in unique(Dargs$IdM)){
  yobsi = Dargs$yobs[Dargs$IdM==i]
  ni = length(yobsi)
  
  dphi_mu0 = matrix(c(rep(0,ni),1/omega0.2,0,0,0))
  dphi_mu1 = matrix(c(rep(0,ni),0,1/omega1.2,0,0))
  dphi_omega0 = matrix(c(rep(0,ni),-mu0/(omega0.2**2),0,1/(2*(omega0.2**2)),0))
  dphi_omega1 = matrix(c(rep(0,ni),0,-mu1/(omega1.2**2),0,1/(2*(omega1.2**2))))
  dphi_sigma = matrix(c(rep(2/sigma**3,ni),0,0,0,0))
  
  Si = as.matrix(c(1/2 * (yobsi-pred[Dargs$IdM==i])**2, suffStat$statphi1[i,], suffStat$statphi3[i,]))
  
  delta_mu0 = sum(Si*dphi_mu0) - mu0/omega0.2
  delta_mu1 = sum(Si*dphi_mu1) - mu1/omega1.2
  delta_omega0 = sum(Si*dphi_omega0) -1/(2*omega0.2) + mu0**2/(2*omega0.2**2)
  delta_omega1 = sum(Si*dphi_omega1) -1/(2*omega1.2) + mu1**2/(2*omega1.2**2)
  delta_sigma = sum(Si*dphi_sigma) -1/(sigma) *ni
  
  mat = matrix(c(delta_mu0,delta_mu1,delta_omega0,delta_omega1,delta_sigma),ncol=1)
  bigmat = bigmat + mat%*%t(mat)
}

bigmat = bigmat/length(unique(Dargs$IdM))

sol = solve(bigmat)
sqrt(diag(sol))/10


#### test gradient numerique

phi = matrix(c(rep(-1/sigma**2,ni),mu0/omega0.2,mu1/omega1.2,-1/(2*omega0.2),-1/(2*omega1.2)))
dphi0 = (matrix(c(rep(-1/sigma**2,ni),(mu0+0.01)/omega0.2,mu1/omega1.2,-1/(2*omega0.2),-1/(2*omega1.2)))-phi)/0.01

psi = ni*log(sqrt(2*pi*sigma**2))+log(sqrt(2*pi*omega0.2))+log(sqrt(2*pi*omega1.2))+mu0**2/(2*omega0.2)+mu1**2/(2*omega1.2)
dpsi0 = (ni*log(sqrt(2*pi*sigma**2))+log(sqrt(2*pi*omega0.2))+log(sqrt(2*pi*omega1.2))+(mu0+0.01)**2/(2*omega0.2)+mu1**2/(2*omega1.2)-(ni*log(sqrt(2*pi*sigma**2))+log(sqrt(2*pi*omega0.2))+log(sqrt(2*pi*omega1.2))+mu0**2/(2*omega0.2)+mu1**2/(2*omega1.2)))/0.01
dpsi0
### ici on calcule les SE: methode 2, on dérive par rapport à omega 

phibis = matrix(phi,nrow = 100, ncol=2)
pssi = transphi(phibis,Dargs$transform.par)
pred = Dargs$structural.model(suffStat$statphi1,Dargs$IdM,Dargs$XM)

mu0 = l1[1]
mu1 = l1[2]
omega0.2 = var.eta[Uargs$i1.omega2][1]
omega1.2 = var.eta[Uargs$i1.omega2][2]
sigma = varList$pres[Uargs$ind.res]

bigmat=matrix(data=0,nrow = 5,ncol = 5)

for (i in unique(Dargs$IdM)){
  yobsi = Dargs$yobs[Dargs$IdM==i]
  ni = length(yobsi)
  
  dphi_mu0 = matrix(c(rep(0,ni),1/omega0.2,0,0,0))
  dphi_mu1 = matrix(c(rep(0,ni),0,1/omega1.2,0,0))
  dphi_omega0 = matrix(c(rep(0,ni),-2*mu0/(sqrt(omega0.2)**3),0,1/(sqrt(omega0.2)**3),0))
  dphi_omega1 = matrix(c(rep(0,ni),0,-2*mu1/(sqrt(omega1.2)**3),0,1/(sqrt(omega1.2))**3))
  dphi_sigma = matrix(c(rep(2/(sigma**3),ni),0,0,0,0))
  
  Si = as.matrix(c(1/2 * (yobsi-pred[Dargs$IdM==i])**2,  suffStat$statphi1[i,], suffStat$statphi3[i,]))
  
  delta_mu0 = sum(Si*dphi_mu0) - mu0/omega0.2
  delta_mu1 = sum(Si*dphi_mu1) - mu1/omega1.2
  delta_omega0 = sum(Si*dphi_omega0) -1/(sqrt(omega0.2)) + mu0**2/(sqrt(omega0.2)**3)
  delta_omega1 = sum(Si*dphi_omega1) -1/(sqrt(omega1.2)) + mu1**2/(sqrt(omega1.2)**3)
  delta_sigma = sum(Si*dphi_sigma) -1/(sigma) *ni
  
  mat = matrix(c(delta_mu0,delta_mu1,delta_omega0,delta_omega1,delta_sigma),ncol=1)
  bigmat = bigmat + mat%*%t(mat)
}

bigmat = bigmat/length(unique(Dargs$IdM))

sol = solve(bigmat)
sqrt(sol[1,1])
sqrt(sol[2,2])
sqrt(sol[3,3])
sqrt(sol[4,4])
sqrt(sol[5,5])


#################################################
# Ne pas lancer, script automatique # 
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