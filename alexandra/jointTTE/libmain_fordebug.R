# input: model, data, control
#saemix<-function(model,data,control=list()) {
  
convplot.infit<-function(allpar,K1,niter=0) {
  # Convergence plots for all the fixed effects, random effects and residual variability
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1 
  np<-dim(allpar)[2]
  K<-dim(allpar)[1]
  n1<-round(sqrt(np))
  n2<-ceiling(np/n1)
  if(n1>5 | n2>5) {n1<-3;n2<-4}
  if(niter==0) niter<-K
  par(mfrow=c(n1,n2))
  for(j in 1:np) {
    plot(1:niter,allpar[1:niter,j],type="l", xlab="Iteration", ylab=colnames(allpar)[j])
    abline(v=K1)
  }
}

saemixObject<-new(Class="SaemixObject",data=data,model=model,options=control)
#  saemixObject<-new(Class="SaemixObject",data=saemix.data, model=saemix.model,options=saemix.options)
opt.warn<-getOption("warn")
if(!saemixObject["options"]$warnings) options(warn=-1)

saemix.options<-saemixObject["options"]
saemix.model<-saemixObject["model"]
saemix.data<-saemixObject["data"]
saemix.data@ocov<-saemix.data@ocov[saemix.data@data[,"mdv"]==0,,drop=FALSE]
saemix.data@data<-saemix.data@data[saemix.data@data[,"mdv"]==0,]
saemix.data@ntot.obs<-dim(saemix.data@data)[1]
#  showall(saemixObject)

# Initialising random generator
set.seed(saemix.options$seed)

############################################
#  Main Algorithm
############################################

# Initialisation - creating several lists with necessary information extracted (Uargs, Dargs, opt,varList, suffStat)
xinit<-initialiseMainAlgo(saemix.data,saemix.model,saemix.options)
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

if (Dargs$modeltype=="structural"){
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
suffStat<-list(statphi1=0,statphi2=0,statphi3=0,statrese=0)
phi<-array(data=0,dim=c(Dargs$N, Uargs$nb.parameters, saemix.options$nb.chains))

# structural model, check nb of parameters
structural.model<-saemix.model["model"]
#  nb.parameters<-saemix.model["nb.parameters"]

keepAR <- NULL

# Running the algorithm
# hw=waitbar(1,'Estimating the population parameters (SAEM). Wait...');
if(saemix.options$displayProgress) par(ask=FALSE)
if(saemix.options$warnings) cat("Running main SAEM algorithm\n")
if(saemix.options$warnings) print(date())
for (kiter in 1:saemix.options$nbiter.tot) { # Iterative portion of algorithm
#for (kiter in 1:160) { # Iterative portion of algorithm
# for (kiter in 1:50) { # Iterative portion of algorithm
  # Burn-in - resetting sufficient statistics
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
  
  # E-step
  xmcmc<-estep(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
  varList<-xmcmc$varList
  DYF<-xmcmc$DYF
  phiM<-xmcmc$phiM
  if(!is.null(xmcmc$keepAR)) keepAR<-rbind(keepAR, xmcmc$keepAR)
  #  psiM<-transphi(phiM,saemix.model["transform.par"])
  
  # M-step
  if(opt$stepsize[kiter]>0) {
    ############# Stochastic Approximation
    xstoch<-mstep(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat)
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
    
    if(Dargs$modeltype=="structural") {
      allpar[(kiter+1),]<-c(l1,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res])
    } else{
      allpar[(kiter+1),]<-c(l1,var.eta[Uargs$i1.omega2])
    }
    
  } else { #end of loop on if(opt$stepsize[kiter]>0)
    allpar[(kiter+1),]<-allpar[kiter,]
  }
  if(Dargs$modeltype=="structural") {
    theta<-c(fixed.psi,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res])
  } else{
    theta<-c(fixed.psi,var.eta[Uargs$i1.omega2])
  }
  # End of loop on kiter
}

convplot.infit(allpar,saemix.options$nbiter.saemix[1],niter=(kiter-2))

etaM<-xmcmc$etaM # only need etaM here (re-created in estep otherwise)
