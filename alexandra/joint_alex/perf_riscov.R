data = read.table("W:/RISCOV/simu/back7/datas/data1.txt", header = T)
data = data[-c(which(data$ytype %in% c(4,5,6,7))),]

souslongi = data[data$ytype==1 | data$ytype==2 | data$ytype==3,]
soussurv = data[data$ytype==8 | data$ytype==9,]

surv = data.frame(id = unique(data$id), time = 0, obs = 0, ytype=4)
for (i in surv$id){
  if (soussurv$obs[soussurv$id==i & soussurv$ytype==8][2]!=0){
    surv$obs[surv$id==i] = 1
    surv$time[surv$id==i] = soussurv$time[soussurv$ytype==8 & soussurv$id==i][2]
  }
  else if (soussurv$obs[soussurv$id==i & soussurv$ytype==9][2]!=0){
    surv$obs[surv$id==i] = 2
    surv$time[surv$id==i] = soussurv$time[soussurv$ytype==9 & soussurv$id==i][2]
  }
  else {
    surv$time[surv$id==i] = 30
  }
}

data_joint = rbind(souslongi,surv)


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

dataJM<-saemixData(name.data=data_joint, name.group=c("id"), name.predictors=c("time","obs"),
                   name.response="obs",name.ytype = "ytype")

#model JM lin + 2 competing risks
JMmodel<-function(psi,id,xidep) {
  ytype<-xidep$ytype  
  
  b01 <- psi[id,1]
  b11 <- psi[id,2]
  b02 <- psi[id,3]
  b12 <- psi[id,4]
  b03 <- psi[id,5]
  b13 <- psi[id,6]
  b23 <- psi[id,7]
  a3 <- psi[id,8]
  p1 <- psi[id,9] 
  g1 <- psi[id,10]
  alpha1 <- psi[id,11]
  alpha2 <- psi[id,12]
  alpha3 <- psi[id,13]
  
  T2 = xidep[ytype==2,1]
  T3 = xidep[ytype==3,1]
  T<-xidep[ytype==4,1] # vector of times partie survie
  ev = xidep$obs[ytype==4]
  Nj <- length(T)
  cens<-which(ev==0)  # indices of censored observations
  ind1 <- which(ev==1) # indices of event 1
  ind2 <- which(ev==2) # indices of event 2
  
  b01b <- unique(b01)
  b11b <- unique(b11)
  b02b <- unique(b02)
  b12b <- unique(b12)
  b03b = unique(b03)
  b13b = unique(b13)
  b23b = unique(b23)
  a3b = unique(a3)
  p1b <- unique(p1)
  g1b <- unique(g1)
  alpha1b = unique(alpha1)
  alpha2b = unique(alpha2)
  alpha3b = unique(alpha3)
  
  f=function(x) seq(0,x,length.out=100)
  tab = mapply(f,T)
  tab = t(tab)
  pas = tab[,2]-tab[,1]
  
  f2=function(x) seq(0,x,length.out=1000)
  tab2 = replicate(Nj,f2(1000))
  tab2 = t(tab2)
  pas2 = tab2[,2]-tab2[,1]
  
  haz1 = p1b*g1b*exp(-g1b*tab)/(1-p1b*(1-exp(-g1b*tab)))*exp(alpha1b*(b01b+b11b*tab-7.42)+alpha2b*(b02b+b12b*tab-2.89)+alpha3b*(b03b+a3b*(exp(b13b*tab)-exp(b23b*tab))-5.78))
  H1 = apply(haz1,1,sum)*pas
  hazt1 = haz1[,100]
  
  haz1b = p1b*g1b*exp(-g1b*tab2)/(1-p1b*(1-exp(-g1b*tab2)))*exp(alpha1b*(b01b+b11b*tab2-7.42)+alpha2b*(b02b+b12b*tab2-2.89)+alpha3b*(b03b+a3b*(exp(b13b*tab2)-exp(b23b*tab2))-5.78))
  H1b = apply(haz1b,1,sum)*pas2
  F1 = 1-exp(-H1b)
  
  hazt2 = 1/10*(1-F1)*exp(-T/10)/(1-(1-F1)*(1-exp(-T/10)))
  H2 = -log(1-(1-F1)*(1-exp(-T/10)))
  
  logpdf <- rep(0,Nj)
  logpdf[cens] <- log(exp(-H1[cens])+exp(-H2[cens])-1)
  logpdf[ind1] <- -H1[ind1] + log(hazt1[ind1]) 
  logpdf[ind2] <- -H2[ind2] + log(hazt2[ind2]) 
  
  ypred = rep(NA,length(xidep[,1]))
  
  ypred[ytype==1] = b01[ytype==1]+b11[ytype==1]*xidep[ytype==1,1]
  ypred[ytype==2] = b02[ytype==2]+b12[ytype==2]*xidep[ytype==2,1]
  ypred[ytype==3] = b03[ytype==3]+a3[ytype==3]*(exp(b13[ytype==3]*xidep[ytype==3,1])-exp(b23[ytype==3]*xidep[ytype==3,1]))
  
  ypred[ytype==4] = logpdf
  
  return(ypred)
}


# joint TTE  
param<-c(7.4,0.003,4.2,-0.16,4.6,-0.15,-0.16,5.3,0.05,0.1,-11,0.14,0.6) # b01, b11, b02, b12, b03, b13, b23, a3, p1, g1, alpha1, alpha2, alpha3 
jointTTE<-saemixModel(model=JMmodel,description="JM lin+competing risks",modeltype=c("structural","structural","structural","likelihood"),
                      psi0=matrix(param,ncol=13,byrow=TRUE,dimnames=list(NULL, c( "b01", "b11", "b02", "b12", "b03", "b13", "b23", "a3", "p1", "g1", "alpha1", "alpha2", "alpha3"))),
                      transform.par=c(0,0,0,0,0,0,0,1,1,1,0,0,0), covariance.model=diag(c(1,1,1,1,1,1,1,1,0,0,0,0,0)),
                      fixed.estim = c(1,1,1,1,1,1,1,1,1,1,1,1,1), omega.init = diag(c(5,0.1,0.0016,0.0001,0.8,0.022,0.26,0.9,0.1,0.1,1,1,1)),
                      error.model = c("constant","constant","proportional","likelihood"))


## clean run 
source(file.path(workDir,"multi_aux2.R"))
source(file.path(workDir,"multi_initializeMainAlgo.R"))
source(file.path(workDir,"multi_estep.R"))
source(file.path(workDir,"multi_mstep.R"))
source(file.path(workDir,"multi_main.R"))
source(file.path(workDir,"multi_map.R"))

saemix.data<-dataJM
saemix.model<-jointTTE
saemix.options<-saemixControl(seed=12345, map=FALSE, fim=FALSE, ll.is=FALSE, displayProgress = T, nbiter.saemix=c(400,150))
yfit <- saemix.multi(saemix.model, saemix.data, saemix.options)

yfit <- map.saemix(yfit)

source(file.path(workDir,"func_fim_JMCR.R"))

fim_fit = try(fim.saemix(yfit))

d = try(data.frame(par=c(fim_fit@results@name.fixed,fim_fit@results@name.random,fim_fit@results@name.sigma),est = c(fim_fit@results@fixed.effects,diag(fim_fit@results@omega)[1:2],fim_fit@results@respar[1:2]),
                   se = c(fim_fit@results@se.fixed,fim_fit@results@se.omega,fim_fit@results@se.respar)))


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



######################### RESULTS ########################################
est = data.frame(sim = 1:100,b01=NA,b11=NA,b02=NA,b12=NA,b03=NA,b13=NA,b23=NA,a3=NA,p1=NA,g1=NA,alpha11=NA,alpha12=NA,alpha13=NA,err_a1=NA,err_a2=NA,err_b3=NA,
                 omega_b01=NA,omega_b11=NA,omega_b02=NA,omega_b12=NA,omega_b03=NA,omega_b13=NA,omega_b23=NA,omega_a3=NA)

for (i in 1:100){
  t = read.table(paste0("W:/saemix/dev/simu_perf/riscov/res_sim",i,".txt"),sep=' ',header=T)
  est$b01[i] = t$est[t$par=='b01']
  est$b11[i] = t$est[t$par=='b11']
  est$b02[i] = t$est[t$par=='b02']
  est$b12[i] = t$est[t$par=='b12']
  est$b03[i] = t$est[t$par=='b03']
  est$b13[i] = t$est[t$par=='b13']
  est$b23[i] = t$est[t$par=='b23']
  est$a3[i] = t$est[t$par=='a3']
  est$p1[i] = t$est[t$par=='p1']
  est$g1[i] = t$est[t$par=='g1']
  est$alpha11[i] = t$est[t$par=='alpha1']
  est$alpha12[i] = t$est[t$par=='alpha2']
  est$alpha13[i] = t$est[t$par=='alpha3']
  est$err_a1[i] = t$est[t$par=='a1']
  est$err_a2[i] = t$est[t$par=='a2']
  est$err_b3[i] = t$est[t$par=='b3']
  est$omega_b01[i] = t$est[t$par=='omega2.b01']
  est$omega_b11[i] = t$est[t$par=='omega2.b11']
  est$omega_b02[i] = t$est[t$par=='omega2.b02']
  est$omega_b12[i] = t$est[t$par=='omega2.b12']
  est$omega_b03[i] = t$est[t$par=='omega2.b03']
  est$omega_b13[i] = t$est[t$par=='omega2.b13']
  est$omega_b23[i] = t$est[t$par=='omega2.b23']
  est$omega_a3[i] = t$est[t$par=='omega2.a3']
}

e.b01 = (est$b01-b01)/b01
e.b11 = (est$b11-b11)/b11
e.b02 = (est$b02-b02)/b02
e.b12 = (est$b12-b12)/b12
e.b03 = (est$b03-b03)/b03
e.b13 = (est$b13-b13)/b13
e.b23 = (est$b23-b23)/b23
e.a3 = (est$a3-a3)/a3
e.p1 = (est$p1-p1)/p1
e.g1 = (est$g1-g1)/g1
e.alpha11 = (est$alpha11-alpha1)/alpha1
e.alpha12 = (est$alpha12-alpha2)/alpha2
e.alpha13 = (est$alpha13-alpha3)/alpha3
e.err_a1 = (est$err_a1-sigma_a1)/sigma_a1
e.err_a2 = (est$err_a2-sigma_a2)/sigma_a2
e.err_b3 = (est$err_b3-sigma_b3)/sigma_b3

e.omega_b01 = (est$omega_b01-omega_b01)/omega_b01
e.omega_b11 = (est$omega_b11-omega_b11)/omega_b11
e.omega_b02 = (est$omega_b02-omega_b02)/omega_b02
e.omega_b12 = (est$omega_b12-omega_b12)/omega_b12
e.omega_b03 = (est$omega_b03-omega_b03)/omega_b03
e.omega_b13 = (est$omega_b13-omega_b13)/omega_b13
e.omega_b23 = (est$omega_b23-omega_b23)/omega_b23
e.omega_a3 = (est$omega_a3-omega_a3)/omega_a3


########### Biais relatifs ########
(mean(est$b01)-b01)/b01
(mean(est$b11)-b11)/b11
(mean(est$b02)-b02)/b02
(mean(est$b12)-b12)/b12
(mean(est$b03)-b03)/b03
(mean(est$b13)-b13)/b13
(mean(est$b23)-b23)/b23
(mean(est$a3)-a3)/a3
(mean(est$p1)-p1)/p1
(mean(est$g1)-g1)/g1
(mean(est$alpha11)-alpha1)/alpha1
(mean(est$alpha12)-alpha2)/alpha2
(mean(est$alpha13)-alpha3)/alpha3

(mean(est$err_a1)-sigma_a1)/sigma_a1
(mean(est$err_a2)-sigma_a2)/sigma_a2
(mean(est$err_b3)-sigma_b3)/sigma_b3

(mean(est$omega_b01)-omega_b01)/omega_b01
(mean(est$omega_b11)-omega_b11)/omega_b11
(mean(est$omega_b02)-omega_b02)/omega_b02
(mean(est$omega_b12)-omega_b12)/omega_b12
(mean(est$omega_b03)-omega_b03)/omega_b03
(mean(est$omega_b13)-omega_b13)/omega_b13
(mean(est$omega_b23)-omega_b23)/omega_b23
(mean(est$omega_a3)-omega_a3)/omega_a3

m_p1 = mean(est$p1)-p1
s_p1 = sd(est$p1-p1)
m_p1/s_p1

m_g1 = mean(est$g1)-g1
s_g1 = sd(est$g1-g1)
m_g1/s_g1

m_alpha1 = mean(est$alpha11)-alpha1
s_alpha1 = sd(est$alpha11-alpha1)
m_alpha1/s_alpha1

m_alpha2 = mean(est$alpha12)-alpha2
s_alpha2 = sd(est$alpha12-alpha2)
m_alpha2/s_alpha2

m_alpha3 = mean(est$alpha13)-alpha3
s_alpha3 = sd(est$alpha13-alpha3)
m_alpha3/s_alpha3

m_b13 = mean(est$b13)-b13
s_b13 = sd(est$b13-b13)
m_b13/s_b13

m_b23 = mean(est$b23)-b23
s_b23 = sd(est$b23-b23)
m_b23/s_b23

m_a3 = mean(est$a3)-a3
s_a3 = sd(est$a3-a3)
m_a3/s_a3



########### RRMSE ######

sqrt((1/100)*sum(((est$b01-b01)/b01)^2))*100
sqrt(100*sum((est$b01-b01)/b01)**2)
sqrt(mean((est$b01-b01)^2))/b01*100
sqrt(mean((est$b11-b11)^2))/b11
sqrt(mean((est$b02-b02)^2))/b02
sqrt(mean((est$b12-b12)^2))/b12
sqrt(mean((est$b03-b03)^2))/b03
sqrt(mean((est$b13-b13)^2))/b13
sqrt(mean((est$b23-b23)^2))/b23
sqrt(mean((est$a3-a3)^2))/a3

sqrt(mean((est$err_a1-sigma_a1)^2))/sigma_a1
sqrt(mean((est$err_a2-sigma_a2)^2))/sigma_a2
sqrt(mean((est$err_b3-sigma_b3)^2))/sigma_b3

sqrt(mean((est$omega_b01-omega_b01)^2))/omega_b01
sqrt(mean(((est$omega_b01-omega_b01)/omega_b01)^2))
sqrt(mean((est$omega_b11-omega_b11)^2))/omega_b11
sqrt(mean((est$omega_b02-omega_b02)^2))/omega_b02
sqrt(mean((est$omega_b12-omega_b12)^2))/omega_b12
sqrt(mean((est$omega_b03-omega_b03)^2))/omega_b03
sqrt(mean((est$omega_b13-omega_b13)^2))/omega_b13
sqrt(mean((est$omega_b23-omega_b23)^2))/omega_b23
sqrt(mean((est$omega_a3-omega_a3)^2))/omega_a3

sqrt(mean((est$p1-p1)^2))/p1
sqrt(mean((est$g1-g1)^2))/g1
sqrt(mean((est$alpha11-alpha1)^2))/alpha1
sqrt(mean((est$alpha12-alpha2)^2))/alpha2
sqrt(mean((est$alpha13-alpha3)^2))/alpha3


########### violin plots ###########
liste1=c(rep("μ01",nrow(est)),rep("μ11",nrow(est)),rep("μ21",nrow(est)),rep("μa1",nrow(est)),
         rep("omega_01",nrow(est)),rep("omega_11",nrow(est)),
         rep("omega_21",nrow(est)),rep("omega_a1",nrow(est)),
         rep("sigma_b1",nrow(est)))
values1 = c(e.b03,e.b13,e.b23,e.a3,e.omega_b03,e.omega_b13,e.omega_b23,e.omega_a3,
            e.err_b3)
values1= values1*100


tb1 = data.frame(param=liste1,val=values1)
tb1$param=factor(tb1$param, levels = c("μ01", "μ11", "μ21", "μa1", "omega_01", "omega_11", "omega_21","omega_a1","sigma_b1"))


vioplot = ggplot(data=tb1, aes(x = param, y = val)) +
  geom_violin(fill="#99CCFF",width=0.6,scale = "width") + geom_boxplot(width=0.07,outlier.size = 0.5,size=0.2) + 
  geom_hline(yintercept = 0,col='red',linetype='dashed',lwd=0.2)+
  theme_classic()+ylab("Relative error (%)")+xlab("Longitudinal parameters for biomarker 1")+
  scale_x_discrete(labels=c("μ01", "μ11", "μ21", "μa1", "omega_01"=expression(omega*"01"), "omega_11"=expression(omega*"11"), "omega_21"=expression(omega*"21"),"omega_a1"=expression(omega*"a1"),"sigma_b1"=expression(sigma*"b1")))+
  scale_y_continuous(limits = c(-100,200))+
  theme(axis.text.x = element_text(angle=45,hjust=1,size = 8),axis.title = element_text(size=10),axis.text.y = element_text(size=8))
vioplot


liste2=c(rep("μ02",nrow(est)),rep("μ12",nrow(est)),
         rep("omega_02",nrow(est)),rep("omega_12",nrow(est)),
         rep("sigma_a2",nrow(est)))
values2 = c(e.b01,e.b11,e.omega_b01,e.omega_b11,e.err_a1)
values2= values2*100

tb2 = data.frame(param=liste2,val=values2)
tb2$param=factor(tb2$param, levels = c("μ02", "μ12", "omega_02", "omega_12","sigma_a2"))


vioplot2 = ggplot(data=tb2, aes(x = param, y = val)) +
  geom_violin(fill="#99CCFF",width=0.6,scale = "width") + geom_boxplot(width=0.07,outlier.size = 0.5,size=0.2) + 
  geom_hline(yintercept = 0,col='red',linetype='dashed', lwd=0.2)+
  theme_classic()+ylab("Relative error (%)")+xlab("Longitudinal parameters for biomarker 2")+
  scale_x_discrete(labels=c("μ02", "μ12","omega_02"=expression(omega*"02"), "omega_12"=expression(omega*"12"),"sigma_a2"=expression(sigma*"a2")))+
  scale_y_continuous(limits = c(-100,200))+
  theme(axis.text.x = element_text(angle=45,hjust=1,size = 8),axis.title = element_text(size=10),axis.text.y = element_text(size=8))
vioplot2


liste3=c(rep("μ03",nrow(est)),rep("μ13",nrow(est)),
         rep("omega_03",nrow(est)),rep("omega_13",nrow(est)),
         rep("sigma_a3",nrow(est)))
values3 = c(e.b02,e.b12,e.omega_b02,e.omega_b12,e.err_a2)
values3 = values3*100

tb3 = data.frame(param=liste3,val=values3)
tb3$param=factor(tb3$param, levels = c("μ03", "μ13", "omega_03", "omega_13","sigma_a3"))


vioplot3 = ggplot(data=tb3, aes(x = param, y = val)) +
  geom_violin(fill="#99CCFF",width=0.6,scale = "width") + geom_boxplot(width=0.07,outlier.size = 0.5,size=0.2) +
  geom_hline(yintercept = 0,col='red',linetype='dashed',lwd=0.2) +
  theme_classic()+ylab("Relative error (%)")+xlab("Longitudinal parameters for biomarker 3")+
  scale_y_continuous(limits = c(-100,200))+
  scale_x_discrete(labels=c("μ03", "μ13","omega_03"=expression(omega*"03"), "omega_13"=expression(omega*"13"),"sigma_a3"=expression(sigma*"a3")))+
  theme(axis.text.x = element_text(angle=45,hjust=1,size = 8),axis.title = element_text(size=10),axis.text.y = element_text(size=8))
vioplot3


liste4=c(rep("p1",nrow(est)),rep("g1",nrow(est)),
         rep("alpha11",nrow(est)),rep("alpha12",nrow(est)),
         rep("alpha13",nrow(est)))
values4 = c(e.p1,e.g1,e.alpha11,e.alpha12,e.alpha13)
values4 = values4*100

tb4 = data.frame(param=liste4,val=values4)
tb4$param=factor(tb4$param, levels = c("p1", "g1", "alpha11", "alpha12","alpha13"))


vioplot4 = ggplot(data=tb4, aes(x = param, y = val)) +
  geom_violin(fill="#99CCFF",width=0.6,scale = "width") + geom_boxplot(width=0.07,outlier.size = 0.5,size=0.2) +
  geom_hline(yintercept = 0,col='red',linetype='dashed',lwd=0.2) +
  theme_classic()+ylab("Relative error (%)")+xlab("Survival parameters")+
  scale_y_continuous(limits = c(-100,200))+
  scale_x_discrete(labels=c("p1", "g1","alpha11"=expression(alpha*"11"), "alpha12"=expression(alpha*"12"),"alpha13"=expression(alpha*"13")))+
  theme(axis.text.x = element_text(angle=45,hjust=1,size = 8),axis.title = element_text(size=10),axis.text.y = element_text(size=8))
vioplot4


plot_grid(vioplot,vioplot2,vioplot3,vioplot4,nrow=4)

ggsave2("simu_riscov.png",path = "C:/Users/AlexandraLAVALLEY/Documents/These/projet 3/figures",width = 15, height = 20, units = "cm",dpi = 600)
