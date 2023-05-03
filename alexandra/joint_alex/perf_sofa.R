### 100 simu pour voir estimatin des param sur la base de SOFA 

data = read.table("W:/SOFA/simu_2/datas/data1.txt", header = T)
colnames(data) = c("id", "time", "obs", "ytype", "cov")
souslongi = data[data$ytype==1,]
soussurv = data[data$ytype==2 | data$ytype==3,]

surv = data.frame(id = unique(data$id), time = 0, obs = 0, ytype=2, cov2 = 0, cov3 = 0)
for (i in surv$id){
  if (soussurv$obs[soussurv$id==i & soussurv$ytype==2][2]!=0){
    surv$obs[surv$id==i] = 1
    surv$time[surv$id==i] = soussurv$time[soussurv$ytype==2 & soussurv$id==i][2]
  }
  else if (soussurv$obs[soussurv$id==i & soussurv$ytype==3][2]!=0){
    surv$obs[surv$id==i] = 2
    surv$time[surv$id==i] = soussurv$time[soussurv$ytype==3 & soussurv$id==i][2]
  }
  else {
    surv$time[surv$id==i] = 30
  }
  surv$cov2[surv$id==i] = ifelse(soussurv$cov[soussurv$id==i][1]==2,1,0)
  surv$cov3[surv$id==i] = ifelse(soussurv$cov[soussurv$id==i][1]==3,1,0)
  souslongi$cov2[souslongi$id==i] = ifelse(soussurv$cov[soussurv$id==i][1]==2,1,0)
  souslongi$cov3[souslongi$id==i] = ifelse(soussurv$cov[soussurv$id==i][1]==3,1,0)
}

souslongi = souslongi[,-5]
data_joint = rbind(souslongi,surv)

data_joint = data_joint[data_joint$id %in% seq(1,100),]
###############################################################################
###############################################################################

# ou on re simule avec erreur additive pour estimer le modèle correct




###############################################################################
###############################################################################
###############################################################################
###############################################################################


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
data_joint$cov2 = as.factor(data_joint$cov2)
data_joint$cov3 = as.factor(data_joint$cov3)
dataJM<-saemixData(name.data=data_joint, name.group=c("id"), name.predictors=c("time","obs"),
                   name.response="obs",name.ytype = "ytype", name.covariates = c("cov2","cov3"))

#model JM lin + 2 competing risks
JMmodel<-function(psi,id,xidep) {
  ytype<-xidep$ytype  
  
  b0 <- psi[id,1]
  b1 <- psi[id,2]
  b2 <- psi[id,3]
  a <- psi[id,4]
  tlag <- psi[id,5]
  p1 <- psi[id,6] 
  g1 <- psi[id,7]
  alpha1 <- psi[id,8]
  beta <- psi[id,9]
  
  T<-xidep[ytype==2,1] # vector of times partie survie ev 1
  ev = xidep$obs[ytype==2]
  Nj <- length(T)
  cens<-which(ev==0)  # indices of censored observations
  ind1 <- which(ev==1) # indices of event 1
  ind2 <- which(ev==2)
  
  schem = sapply(1:Nj, function(i) sum(floor(T[1:i])+1)+i)
  
  b0b = b0[schem]
  b1b = b1[schem]
  b2b = b2[schem]
  ab = a[schem]
  tlagb = tlag[schem]
  p1b <- p1[schem]
  g1b <- g1[schem]
  alpha1b = alpha1[schem]
  betab = beta[schem]
  
  f=function(x) seq(0,x,length.out=100)
  tab = mapply(f,T)
  tab = t(tab)
  pas = tab[,2]-tab[,1]
  
  f2=function(x) seq(0,x,length.out=1000)
  tab2 = replicate(Nj,f2(1000))
  tab2 = t(tab2)
  pas2 = tab2[,2]-tab2[,1]
  
  haz1 = p1b*g1b*exp(-g1b*tab)/(1-p1b*(1-exp(-g1b*tab)))*exp(pmin(pmax(alpha1b*(b0b+ab*(exp(b1b*(tab-tlagb))-exp(b2b*(tab-tlagb)))),0),24)+betab)
  H1 = apply(haz1,1,sum)*pas
  hazt1 = haz1[,100]
  
  haz1b = p1b*g1b*exp(-g1b*tab2)/(1-p1b*(1-exp(-g1b*tab2)))*exp(pmin(pmax(alpha1b*(b0b+ab*(exp(b1b*(tab2-tlagb))-exp(b2b*(tab2-tlagb)))),0),24)+betab)
  H1b = apply(haz1b,1,sum)*pas2
  F1 = 1-exp(-H1b)
  
  hazt2 = 1/10*(1-F1)*exp(-T/10)/(1-(1-F1)*(1-exp(-T/10)))
  H2 = -log(1-(1-F1)*(1-exp(-T/10)))
  
  logpdf <- rep(0,Nj)
  logpdf[cens] <- log(exp(-H1[cens])+exp(-H2[cens])-1)
  logpdf[ind1] <- -H1[ind1] + log(hazt1[ind1]) 
  logpdf[ind2] <- -H2[ind2] + log(hazt2[ind2]) 
  
  ypred = rep(NA,length(xidep[,1]))
  
  ypred[ytype==1] = pmin(pmax(b0[ytype==1]+a[ytype==1]*(exp(b1[ytype==1]*(xidep[ytype==1,1]-tlag[ytype==1]))-exp(b2[ytype==1]*(xidep[ytype==1,1]-tlag[ytype==1]))),0),24)
  ypred[ytype==2] = logpdf
  
  return(ypred)
}


# joint TTE  
param<-c(2,-0.25,-1,10,-1,0.01,0.05,0.35,0) # b0, b1, b2, a, tlag, p1, g1, alpha1, beta
jointTTE<-saemixModel(model=JMmodel,description="JM lin+competing risks",modeltype=c("structural","likelihood"),
                      psi0=matrix(param,ncol=9,byrow=TRUE,dimnames=list(NULL, c( "b0", "b1", "b2", "a", "tlag","p1", "g1", "alpha1", "beta"))),
                      transform.par=c(0,0,0,1,0,1,1,0,0), covariance.model=diag(c(1,1,1,1,1,0,0,0,0)),
                      fixed.estim = c(1,1,1,1,1,1,1,1,0), omega.init = diag(c(9,0.04,0.25,0.25,0.64,0.1,0.1,0.1,0.1)), error.model = "combined",
                      covariate.model = matrix(c(0,0,0,0,0,0,0,0,1,
                                                 0,0,0,0,0,0,0,0,1),ncol = 9,nrow=2, byrow = T))


## clean run 
source(file.path(workDir,"multi_aux2.R"))
source(file.path(workDir,"multi_initializeMainAlgo.R"))
source(file.path(workDir,"multi_estep.R"))
source(file.path(workDir,"multi_mstep.R"))
source(file.path(workDir,"multi_main.R"))
source(file.path(workDir,"multi_map.R"))

saemix.data<-dataJM
saemix.model<-jointTTE

saemix.options<-saemixControl(seed=12345, map=FALSE, fim=FALSE, ll.is=FALSE)


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


yfit <- saemix.multi(saemix.model, saemix.data, saemix.options)
yfit <- map.saemix(yfit)

source(file.path(workDir,"func_fim_JMCR.R"))

fim_fit = try(fim.saemix(yfit))

d = data.frame(par=c(yfit@results@name.fixed,yfit@results@name.random,yfit@results@name.sigma),est = c(yfit@results@fixed.effects,diag(yfit@results@omega)[1:5],yfit@results@respar[1:2]))


for (i in 1:100){
  code = readLines("W:/saemix/dev/simu_perf/sofa/2000it/code/code.R")
  code[25] = paste0("data = read.table(\"/home/alexandra.lavalley2/SOFA/simu_2/datas/data",i,".txt\", header = T)")
  code[156] = paste0("  saemix.options<-saemixControl(seed=12345, map=FALSE, fim=FALSE, ll.is=FALSE, nb.chains = 1, directory = paste0(\"sofa\",",i,"), displayProgress = T, nbiter.saemix=c(2000,500))")
  code[161] = paste0("  write.table(d,file=\"/home/alexandra.lavalley2/saemix/dev/simu_perf/sofa/2000it/res_sim",i,".txt\",row.names=F)")
  writeLines(code,paste0("W:/saemix/dev/simu_perf/sofa/2000it/code/code",i,".R"))
}

############################# RESULTS 

beta0.pop<-2
beta1.pop<- -0.25
beta2.pop <- -1 
a.pop<- 10
tlag.pop = -1
alpha1<-0.35
alpha2=0

omega_b0<-3
omega_b1<-0.2
omega_b2<-0.5
omega_a<-0.5
omega_tlag = 0.8

err_a<-0.9
err_b = 0.1

p=0.01
g = 0.05

N=1000

est = data.frame(sim = 1:100,b0=NA,b1=NA,b2=NA,a=NA,tlag=NA,p1=NA,g1=NA,alpha1=NA,alpha2=NA,cov1Age2=NA,cov1Age3=NA,cov2Age2=NA,cov2Age3=NA,err_a=NA,err_b=NA,
                 omega_b0=NA,omega_b1=NA,omega_b2=NA,omega_a=NA,omega_tlag=NA)
for (i in 1:nrow(est)){
  t = read.table(paste0("W:/saemix/dev/simu_perf/sofa/1000it/res_sim",i,".txt"),sep=' ',header=T)
  est$b0[i] = t$est[t$par=='b0']
  est$b1[i] = t$est[t$par=='b1']
  est$b2[i] = t$est[t$par=='b2']
  est$a[i] = t$est[t$par=='a']
  est$tlag[i] = t$est[t$par=='tlag']
  est$p1[i] = t$est[t$par=='p1']
  est$g1[i] = t$est[t$par=='g1']
  est$alpha1[i] = t$est[t$par=='alpha1']
  est$cov1Age2[i] = t$est[t$par=='beta_cov2(beta)']
  est$cov1Age3[i] = t$est[t$par=='beta_cov3(beta)']
  est$err_a[i] = t$est[t$par=='a.1']
  est$err_b[i] = t$est[t$par=='b.1']
  est$omega_b0[i] = t$est[t$par=='omega2.b0']
  est$omega_b1[i] = t$est[t$par=='omega2.b1']
  est$omega_b2[i] = t$est[t$par=='omega2.b2']
  est$omega_a[i] = t$est[t$par=='omega2.a']
  est$omega_tlag[i] = t$est[t$par=='omega2.tlag']
}

e.b0 = (est$b0-beta0.pop)/beta0.pop
e.b1 = (est$b1-beta1.pop)/beta1.pop
e.b2 = (est$b2-beta2.pop)/beta2.pop
e.a = (est$a-a.pop)/a.pop
e.tlag = (est$tlag-tlag.pop)/tlag.pop
e.p1 = (est$p1-p)/p
e.g1 = (est$g1-g)/g
e.alpha1 = (est$alpha1-alpha1)/alpha1
e.age2 = (est$cov1Age2-0.346)/0.346
e.age3 = (est$cov1Age3-0.74)/0.74 
e.err_a = (est$err_a-err_a)/err_a
e.err_b = (est$err_b-err_b)/err_b
e.omega_b0 = (est$omega_b0-omega_b0**2)/(omega_b0**2)
e.omega_b1 = (est$omega_b1-omega_b1**2)/(omega_b1**2)
e.omega_b2 = (est$omega_b2-omega_b2**2)/omega_b2**2
e.omega_a = (est$omega_a-omega_a**2)/omega_a**2
e.omega_tlag = (est$omega_tlag-omega_tlag**2)/omega_tlag**2

# biais relatifs
(mean(est$b0)-beta0.pop)/beta0.pop
(mean(est$b1)-beta1.pop)/beta1.pop
(mean(est$b2)-beta2.pop)/beta2.pop
(mean(est$a)-a.pop)/a.pop
(mean(est$tlag)-tlag.pop)/tlag.pop
(mean(est$p1)-p)/p
(mean(est$g1)-g)/g
(mean(est$alpha1)-alpha1)/alpha1
(mean(est$cov1Age2)-0.346)/0.346
(mean(est$cov1Age3)-0.74)/0.74 
(mean(est$err_a)-err_a)/err_a
(mean(est$err_b)-err_b)/err_b
(mean(est$omega_b0)-omega_b0**2)/omega_b0**2
(mean(est$omega_b1)-omega_b1**2)/omega_b1**2
(mean(est$omega_b2)-omega_b2**2)/omega_b2**2
(mean(est$omega_a)-omega_a**2)/omega_a**2
(mean(est$omega_tlag)-omega_tlag**2)/omega_tlag**2

# RMSE
sqrt((1/200)*sum((est$b0-beta0.pop)^2))
sqrt(mean((est$b0-beta0.pop)^2))
sqrt(mean((est$b1-beta1.pop)^2))
sqrt(mean((est$b2-beta2.pop)^2))
sqrt(mean((est$a-a.pop)^2))
sqrt(mean((est$tlag-tlag.pop)^2))
sqrt(mean((est$err_a-err_a)^2))
sqrt(mean((est$err_b-err_b)^2))

sqrt(mean((est$omega_b0-omega_b0)^2))
sqrt(mean((est$omega_b1-omega_b1)^2))
sqrt(mean((est$omega_b2-omega_b2)^2))
sqrt(mean((est$omega_a-omega_a)^2))
sqrt(mean((est$omega_tlag-omega_tlag)^2))

sqrt(mean((est$alpha1-alpha1)^2))
sqrt(mean((est$cov1Age2-0.346)^2))
sqrt(mean((est$cov1Age3-0.74)^2)) 


liste=c()
liste=c(rep("μ0",nrow(est)),rep("μ1",nrow(est)),rep("μ2",nrow(est)),rep("a",nrow(est)),
        rep("tlag",nrow(est)),rep("omega_0",nrow(est)),rep("omega_1",nrow(est)),
        rep("omega_2",nrow(est)),rep("omega_a",nrow(est)),rep("omega_tlag",nrow(est)),
        rep("sigma_a",nrow(est)),rep("sigma_b",nrow(est)),rep("alpha1",nrow(est)),
        rep("β1_age2",nrow(est)),rep("β1_age3",nrow(est)))
values = c(e.b0,e.b1,e.b2,e.a,e.tlag,e.omega_b0,e.omega_b1,e.omega_b2,e.omega_a,e.omega_tlag,
           e.err_a,e.err_b,e.alpha1,e.age2,e.age3)
values= values*100

liste_longi = c(rep("beta0",nrow(est)),rep("beta1",nrow(est)),rep("beta2",nrow(est)),rep("a",nrow(est)),
                rep("tlag",nrow(est)),rep("omega_0",nrow(est)),rep("omega_1",nrow(est)),
                rep("omega_2",nrow(est)),rep("omega_a",nrow(est)),rep("omega_tlag",nrow(est)),
                rep("σ_a",nrow(est)),rep("σ_b",nrow(est)))
values_longi = c(e.b0,e.b1,e.b2,e.a,e.tlag,e.omega_b0,e.omega_b1,e.omega_b2,e.omega_a,e.omega_tlag,
                 e.err_a,e.err_b)
values_longi= values_longi*100

liste_survie=c(rep("alpha1",nrow(est)),rep("gamma1_age2",nrow(est)),rep("gamma1_age3",nrow(est)))
values_survie = c(e.alpha1,e.age2,e.age3)
values_survie= values_survie*100

tb = data.frame(param=liste,val=values)
tb$param=factor(tb$param, levels = c("μ0", "μ1", "μ2", "a", "tlag", "omega_0", "omega_1", "omega_2","omega_a","omega_tlag","sigma_a", "sigma_b","alpha1","β1_age2","β1_age3"))

tb_longi = data.frame(param=liste_longi,val=values_longi)
tb_survie = data.frame(param=liste_survie,val=values_survie)

vioplot = ggplot(data=tb, aes(x = param, y = val)) + geom_hline(yintercept = 0,col='red',linetype='dashed') +
  geom_violin(fill="#99CCFF",width=0.6,scale = "width") + geom_boxplot(width=0.07,outlier.size = 0.5,size=0.2) +
  theme_classic()+ylab("Relative error (%)")+xlab("Parameters (joint model)")+
  scale_x_discrete(labels=c("μ0", "μ1", "μ2", "a", "tlag", "omega_0"=expression(omega*"0"), "omega_1"=expression(omega*"1"), "omega_2"=expression(omega*"2"),"omega_a"=expression(omega*"a"),"omega_tlag"=expression(omega*"tlag"),"sigma_a"=expression(sigma*"a"), "sigma_b"=expression(sigma*"b"),"alpha1"=expression(alpha*"1"),"β1_age2"=expression(beta*"1_age2"),"β1_age3"=expression(beta*"1_age3")))+
  theme(axis.text.x = element_text(angle=45,hjust=1,size = 5),axis.title = element_text(size=10),axis.text.y = element_text(size=8))+
  scale_y_continuous(limits = c(-200,200))
vioplot
ggsave2("simu_sofa.png",path = "C:/Users/AlexandraLAVALLEY/Documents/These/projet 3/figures",width = 15, height = 10, units = "cm",dpi = 600)


vioplot = ggplot(data=tb, aes(x = param, y = val)) + geom_hline(yintercept = 0,col='red',linetype='dashed',size=0.2) +
  geom_violin(fill="#99CCFF",width=0.7,scale = "width",size=0.2) + geom_boxplot(width=0.07,outlier.size = 0.5,size=0.2) +
  theme_classic()+ylab("Relative error (%)")+xlab("Parameters")+
  scale_x_discrete(limits=c("μ0", "μ1", "μ2", "a", "tlag", "omega_0", "omega_1", "omega_2","omega_a","omega_tlag","sigma_a", "sigma_b","alpha1","β1_age[60,75]","β1_age3"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,size = 5),axis.title = element_text(size=10),axis.text.y = element_text(size=8))
vioplot

ggsave2("simubis.png",path = "C:/Users/AlexandraLAVALLEY/Documents/These/article 1/figures",width = 15, height = 10, units = "cm",dpi = 600)


vioplot_longi = ggplot(data=tb_longi, aes(x = param, y = val)) + geom_hline(yintercept = 0,col='red',linetype='dashed') +
  geom_violin(fill="#99CCFF",width=1.5) + geom_boxplot(width=0.06) +
  theme_classic()+ylab("Relative error (%)")+xlab("Parameters")+
  scale_x_discrete(limits=c("b0", "b1", "b2", "a", "tlag", "omega_b0", "omega_b1", "omega_b2","omega_a","omega_tlag","sigma_a", "sigma_b"))+
  theme(axis.text.x = element_text(angle=45,hjust=1))+ylim(-50,50)
vioplot_longi

vioplot_survie = ggplot(data=tb_survie, aes(x = param, y = val)) + geom_hline(yintercept = 0,col='red',linetype='dashed') +
  geom_violin(fill="#99CCFF",width=1) + geom_boxplot(width=0.05) +
  theme_classic()+ylab("Relative error (%)")+xlab("Parameters")+
  scale_x_discrete(limits=c("alpha1","gamma1_age2","gamma1_age3"))+
  theme(axis.text.x = element_text(angle=45,hjust=1))
vioplot_survie

