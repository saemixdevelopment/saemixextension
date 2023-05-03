#################################################
saemixDir <- "C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension"
workDir <- file.path(saemixDir, "alexandra","joint_alex")
setwd(workDir)

#library(Cairo)
library("viridis")
library(rlang)

# Chargement des fonctions originelles de la librairie
progDir<-file.path(saemixDir, "R")
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"func_plots.R")) # for saemix.plot.setoptions

source(file.path(workDir,"multi_aux2.R"))
source(file.path(workDir,"multi_initializeMainAlgo.R"))
source(file.path(workDir,"multi_estep.R"))
source(file.path(workDir,"multi_mstep.R"))
source(file.path(workDir,"multi_main.R"))
source(file.path(workDir,"multi_map.R"))
################################################# Data and model (original files)
# Creating data and model objects

####### NLMEM #####

  data_longi<- read.table(paste0("W:/saemix/dev/comput_se/nonlin/datas/data",it,".txt"), header=TRUE)
  datalongi<-saemixData(name.data=data_longi, name.group=c("id"), name.predictors=c("time"), 
                        name.response="obs")
  # Simulation parameters
  param<-c(10,0.15,0.20,10)
  omega.sim<-c(4, 0.1, 0.05, 0.06)
  sigma.sim <- 0.9
  
  # model
  mod<-function(psi,id,xidep) {
    b0 <- psi[id,1] 
    b1 <- psi[id,2]
    b2 <- psi[id,3]
    a <- psi[id,4]
    time<-xidep[,1] 
    
    ypred <- b0+a*(exp(-b1*(time))-exp(-b2*(time)))
    return(ypred)
  }
  
  
  # Proportional error model
  modlongi<-saemixModel(model=mod,description="longi seul",modeltype="structural",
                        psi0=matrix(param,ncol=4,byrow=TRUE,dimnames=list(NULL, c("b0","b1","b2","a"))),
                        transform.par=c(0,1,1,0), covariance.model=diag(c(1,1,1,1)),fixed.estim = c(1,1,1,1),
                        omega.init = diag(omega.sim),error.model = "constant", error.init = c(sigma.sim,0))
  
  saemix.data<-datalongi
  saemix.model<-modlongi
  saemix.options<-saemixControl(seed=12345, map=F, fim=T, ll.is=FALSE, nb.chains = 3)
  yfit <- saemix.multi(saemix.model, saemix.data, saemix.options)
  tab = data.frame(par=c(yfit@results@name.fixed,yfit@results@name.random,yfit@results@name.sigma[1]),est = c(yfit@results@fixed.effects,diag(yfit@results@omega)[1:4],yfit@results@respar[1]), se = sqrt(diag(yfit@results@fim)))
  

##### JM LMEM 
  data_joint<- read.table(paste0("W:/saemix/dev/comput_se/type1err/jointlin/datas/data",it,".txt"), header=TRUE)
  dataJM<-saemixData(name.data=data_joint, name.group=c("id"), name.predictors=c("time","obs"),
                     name.response="obs",name.ytype = "ytype")
  
  # Simulation parameters
  param<-c(4,0.2,0.025,0.001)
  omega.sim<-c(3, 0.1, 0.00004, 0.00001)
  sigma.sim <- 1
  
  # model
  JMmodel<-function(psi,id,xidep) {
    ytype<-xidep$ytype  # type of response (1: continuous, 2: event)
    b0 <- psi[id,1] 
    b1 <- psi[id,2] 
    h0 <- psi[id,3]
    alpha <- psi[id,4]  ## coeff de lien
    
    ypred <- b0+b1*xidep[,1]  ## pred longi 
    
    T<-xidep[ytype==2,1]# vector of times partie survie
    Nj <- length(T)
    ev = xidep$obs[ytype==2]
    cens<-which(ev==0)  # censoring time=30
    ind <- which(ev==1)
    b0b = b0[ytype==2]
    b1b = b1[ytype==2]
    h0b = h0[ytype==2]
    alphab = alpha[ytype==2]
    
    haz <- h0b*exp(alphab*(b0b+b1b*T))
    H <- (h0b/(alphab*b1b))*exp((b0b+b1b*T)*alphab)-(h0b/(alphab*b1b))*exp(alphab*b0b)
    
    logpdf <- rep(0,Nj)
    logpdf[cens] <- -H[cens] #+ H[cens-1]
    logpdf[ind] <- -H[ind] + log(haz[ind]) #+ H[ind-1]
    
    ypred[ytype==2] = logpdf
    return(ypred)
  }
  
  # Proportional error model
  
  jointTTE<-saemixModel(model=JMmodel,description="JM lin longi one tte",modeltype=c("structural","likelihood"),
                        psi0=matrix(param,ncol=4,byrow=TRUE,dimnames=list(NULL, c("b0","b1","h0","alpha"))),
                        transform.par=c(0,0,1,0), covariance.model=diag(c(1,1,0,0)),
                        fixed.estim = c(1,1,1,1),error.model = "constant",
                        omega.init = diag(omega.sim))
  
  saemix.data<-dataJM
  saemix.model<-jointTTE
  saemix.options<-saemixControl(seed=12345, map=F, fim=T, ll.is=FALSE, nb.chains = 1)
  
  ######### JOINT NLMEM CR
  data_joint<- read.table(paste0("W:/saemix/dev/comput_se/type1err/jointnonlinCR/datas/data",it,".txt"), header=TRUE)
  dataJM<-saemixData(name.data=data_joint, name.group=c("id"), name.predictors=c("time","obs"),
                     name.response="obs",name.ytype = "ytype")
  
  param<-c(4,0.15,0.2,10,0.5,0.10,0.0001,15)
  omega.sim<-c(4, 0.1, 0.1, 0.5, 0.005, 0.001, 0.000001, 0.15)
  sigma.sim <- 1
  
  #model JM nonlin + 2 competing risks
  JMmodel<-function(psi,id,xidep) {
    ytype<-xidep$ytype  
    
    b0 <- psi[id,1]
    b1 <- psi[id,2]
    b2 <- psi[id,3]
    a <- psi[id,4]
    p1 <- psi[id,5] 
    g1 <- psi[id,6]
    alpha1 <- psi[id,7]
    b <- psi[id,8]
    
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
    p1b <- p1[schem]
    g1b <- g1[schem]
    alpha1b = alpha1[schem]
    bb = b[schem]
    
    f=function(x) seq(0,x,length.out=100)
    tab = mapply(f,T)
    tab = t(tab)
    pas = tab[,2]-tab[,1]
    
    f2=function(x) seq(0,x,length.out=1000)
    tab2 = replicate(Nj,f2(1000))
    tab2 = t(tab2)
    pas2 = tab2[,2]-tab2[,1]
    
    haz1 = p1b*g1b*exp(-g1b*tab)/(1-p1b*(1-exp(-g1b*tab)))*exp(alpha1b*(b0b+ab*(exp(-b1b*(tab))-exp(-b2b*(tab)))))
    H1 = apply(haz1,1,sum)*pas
    hazt1 = haz1[,100]
    
    haz1b = p1b*g1b*exp(-g1b*tab2)/(1-p1b*(1-exp(-g1b*tab2)))*exp(alpha1b*(b0b+ab*(exp(-b1b*(tab2))-exp(-b2b*(tab2)))))
    H1b = apply(haz1b,1,sum)*pas2
    F1 = pmin(1-exp(-H1b),0.99999)
    
    hazt2 = 1/bb*(1-F1)*exp(-T/bb)/(1-(1-F1)*(1-exp(-T/bb)))
    H2 = -log(1-(1-F1)*(1-exp(-T/bb)))
    
    logpdf <- rep(0,Nj)
    logpdf[cens] <- log(exp(-H1[cens])+exp(-H2[cens])-1)
    logpdf[ind1] <- -H1[ind1] + log(hazt1[ind1]) 
    logpdf[ind2] <- -H2[ind2] + log(hazt2[ind2]) 
    
    ypred = rep(NA,length(xidep[,1]))
    
    ypred[ytype==1] = b0[ytype==1]+a[ytype==1]*(exp(-b1[ytype==1]*(xidep[ytype==1,1]))-exp(-b2[ytype==1]*(xidep[ytype==1,1])))
    ypred[ytype==2] = logpdf
    
    return(ypred)
  }

  jointTTE<-saemixModel(model=JMmodel,description="JM lin+competing risks",modeltype=c("structural","likelihood"),
                        psi0=matrix(param,ncol=8,byrow=TRUE,dimnames=list(NULL, c( "b0", "b1", "b2", "a", "p1", "g1", "alpha1", "b"))),
                        transform.par=c(0,1,1,1,3,1,0,0), covariance.model=diag(c(1,1,1,1,0,0,0,0)),
                        fixed.estim = c(1,1,1,1,1,1,1,1), omega.init = diag(omega.sim))
  
  saemix.data<-dataJM
  saemix.model<-jointTTE
  saemix.options<-saemixControl(seed=12345, map=F, fim=T, ll.is=FALSE, nb.chains = 3)

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
  deltai = xinit$deltai
  
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
    #print(theta)
  }
  
  ################################################# SA iterations after burn-in (kiter=6 to 150)
  for (kiter in 6:149) {
    xmcmc<-estep.multi(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
    varList<-xmcmc$varList
    DYF<-xmcmc$DYF
    phiM<-xmcmc$phiM
    
    if(opt$stepsize[kiter]>0) {
      ############# Stochastic Approximation
      xstoch<-mstep.multi(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat, deltai)
      varList<-xstoch$varList
      mean.phi<-xstoch$mean.phi
      phi<-xstoch$phi
      betas<-xstoch$betas
      suffStat<-xstoch$suffStat
      deltai = xstoch$deltai
      
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
    if (sum(is.na(deltai))>0) print(kiter)
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
      xstoch<-mstep.multi(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat,deltai)
      varList<-xstoch$varList
      mean.phi<-xstoch$mean.phi
      phi<-xstoch$phi
      betas<-xstoch$betas
      suffStat<-xstoch$suffStat
      deltai = xstoch$deltai
      
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
    if (sum(is.na(deltai))>0) print(kiter)
  }
  ################################################# Smoothing phase (kiter=301 to 400)
  for (kiter in 301:saemix.options$nbiter.tot) {
    xmcmc<-estep.multi(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
    varList<-xmcmc$varList
    DYF<-xmcmc$DYF
    phiM<-xmcmc$phiM
    
    if(opt$stepsize[kiter]>0) {
      ############# Stochastic Approximation
      xstoch<-mstep.multi(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat, deltai)
      varList<-xstoch$varList
      mean.phi<-xstoch$mean.phi
      phi<-xstoch$phi
      betas<-xstoch$betas
      suffStat<-xstoch$suffStat
      deltai = xstoch$deltai
      
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
    if (sum(is.na(deltai))>0) print(kiter)
  }
  
  d = array(data=NA,dim=c(9,9,100))
  for (i in 1:100){
    ddi = deltai[,i] %*% t(deltai[,i])
    d[,,i] = ddi
  }
  
  mat = 0
  for (i in 1:100){
    mat = mat+d[,,i]
  }
  inv=solve(mat)
  sqrt(diag(inv))
  
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
  
  sol = solve(bigmat)  