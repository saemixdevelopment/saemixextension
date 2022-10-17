
############################### Initialising main algorithm #############################

initialiseMainAlgo.alex<-function(saemix.data,saemix.model,saemix.options) {
  n_mark = 2
  # Function to reformat covariance structure and initialise lists used in the main algorithm
  # Input: data, model and options
  # Output:
  ### saemix.model: added elements betaest, covariate model, indices (); modified/formatted psi0 adjusting to the nb of covariates
  ### Dargs: data elements - passed on to functions, unchanged
  ### Uargs: list of indices and variables (fixed) - passed on to functions, unchanged
  ### varList: variability-related elements - passed to functions and optimised
  ### opt: list of options and settings (fixed) - passed on to functions, unchanged
  ### DYF: used for the acceptance/rejection algorithm
  ### parameters optimised during the fit: phiM, mean.phi, betas, fixedpsi.ini
  ### allpar0: array holding the successive values of population parameters
  
  # Elements of the lists
  #   Dargs<-list(IdM=IdM, XM=XM, yM=yM, NM=NM, N=N, nobs=saemix.data["ntot.obs"],
  #               yobs=saemix.data["data"][,saemix.data["name.response"]],transform.par=saemix.model["transform.par"],
  #               error.model=saemix.model["error.model"],structural.model=structural.model)
  #   Uargs<-list(nchains=saemix.options$nb.chains,nb.parameters=nb.parameters, nb.betas=nb.betas, nb.etas=nb.etas, 
  #      t,indx.betaC=indx.betaC, indx.betaI=indx.betaI, ind.res=ind.res,
  #               indest.omega=         nb.parest=nb.paresindest.omega, i0.omega2=i0.omega2, i1.omega2=i1.omega2,	j.covariate=j.covariate, 
  #               ind.fix10=ind.fix10, ind.fix11=ind.fix11, ind.fix1=ind.fix1, ind.fix0=ind.fix0,
  #               MCOV0=MCOV0, COV=COV, COV0=COV0, COV1=COV1, LCOV=LCOV, COV2=COV2, dstatCOV=dstatCOV, 
  #               Mcovariates=Mcovariates, ind.ioM=ind.ioM)
  #   varList<-list(pres=pres,ind0.eta=ind0.eta,ind.eta=ind.eta,omega=omega, MCOV=MCOV,
  #                 domega2=do.call(cbind,rep(list((sqrt(mydiag(omega.eta)))*saemix.options$rw.ini),nb.etas)),diag.omega=mydiag(omega))
  #   opt<-list(stepsize.rw=saemix.options$stepsize.rw,stepsize=stepsize,
  #             proba.mcmc=saemix.options$proba.mcmc,nbiter.mcmc=saemix.options$nbiter.mcmc,
  #             nbiter.sa=saemix.options$nbiter.sa,nbiter.map=saemix.options$nbiter.map,alpha1.sa=saemix.options$alpha.sa,
  #             alpha0.sa=10^(-3/saemix.options$nbiter.sa),nbiter.saemix=saemix.options$nbiter.saemix,
  #             maxim.maxiter=saemix.options$maxim.maxiter,flag.fmin=flag.fmin)
  
  structural.model<-saemix.model["model"]
  nb.parameters<-saemix.model["nb.parameters"]
  N<-saemix.data["N"]
  
  # Initialising residual error model
  # error models :
  #   constant            y = f + a*e
  #   proportional        y = f + b*f*e
  #   combined            y = f + sqrt(a^2+b^2*f^2)*e
  #   exponential         y = f*exp(a*e)    ( <=>  log(y) = log(f) + a*e )
  # error models are a + bf described by [a b], [1]=constant coefficient, [2]= proportional coefficient
  if(saemix.model["modeltype"]=="structural"){
    pres<-c(saemix.model["error.init"][1],saemix.model["error.init"][2],saemix.model["error.init"][3],saemix.model["error.init"][4])
  }
  
  # ECO TODO: integrate all this section in the object creation ?
  # Initialisation: 
  # create local copies modified of omega.init and covariate.model in saemix.model
  # A la fin: i1.omega2 renomme en indx.omega et ind.res en indx.res
  i0.omega2<-which((1-mydiag(saemix.model["covariance.model"]))>0) # index of parameters without IIV
  indest.omega<-which(saemix.model["covariance.model"]>0)
  #  i1.omega2<-which(mydiag(saemix.model$covariance.model)>0)
  i1.omega2<-saemix.model@indx.omega # index of parameters with IIV
  ind.res<-saemix.model["indx.res"]
  
  # Covariate model & design matrix
  id<-saemix.data["data"][,saemix.data["name.group"]]
  if(length(saemix.data["name.covariates"])==0) tab<-data.frame(id=id) else
    tab<-data.frame(id=id,saemix.data["data"][, saemix.data["name.covariates",drop=FALSE]])
  temp2<-unique(tab)
  temp<-tab[!duplicated(id),,drop=FALSE]
  if(dim(temp)[1]!=dim(temp2)[1]) {
    if(saemix.options$warnings) cat("Some covariates have time-varying values; only the first is taken into account in the current version of the algorithm.\n")
  }
  #temp<-temp[order(temp[,1]),]
  if(length(saemix.data["name.covariates"])>0) {
    Mcovariates<-data.frame(id=rep(1,N),temp[,2:dim(temp)[2]])} else {
      Mcovariates<-data.frame(id=rep(1,N))
    }
  # removing from model unused lines
  j.cov<-which(rowSums(saemix.model["betaest.model"])>0)
  betaest.model<-saemix.model["betaest.model"][j.cov,,drop=FALSE]
  Mcovariates<-Mcovariates[,j.cov,drop=FALSE] # eliminate all the unused covariates
  for(icol in dim(Mcovariates)[2])
    if(is.factor(Mcovariates[,icol])) Mcovariates[,icol]<-as.numeric(Mcovariates[,icol])-1
  
  #  if(length(j.cov)==1) {
  #    betaest.model<-matrix(betaest.model,nrow=1, dimnames=list(c("Fixed"),colnames(saemix.model["betaest.model"])))
  #    Mcovariates<-matrix(Mcovariates)
  #  }
  saemix.model["betaest.model"]<-betaest.model
  temp1<-betaest.model[-c(1),,drop=FALSE]
  #  if(is.null(dim(temp1))) temp1<-matrix(temp1,nrow=1, dimnames=list(rownames(betaest.model)[-c(1)], colnames(betaest.model)))  
  saemix.model["covariate.model"]<-temp1
  
  fixedpsi.ini<-saemix.model["psi0"][1,] # initial fixed effects (original parametrization)
  betaI.ini<-transpsi(matrix(fixedpsi.ini,nrow=1),saemix.model["transform.par"]) #initial fixed effects (Gaussian parametrization)
  fixed.ini<-saemix.model["betaest.model"]*0
  fixed.ini[1,]<-betaI.ini
  nr.psi0<-dim(saemix.model["psi0"])[1]
  nr.cov<-dim(saemix.model["betaest.model"])[1]
  if(nr.psi0>nr.cov) {
    saemix.model["psi0"]<-saemix.model["psi0"][1:nr.cov,,drop=FALSE]
    nr.psi0<-dim(saemix.model["psi0"])[1]
  }
  #t1<-NULL
  if(nr.psi0<nr.cov) {
    #  t1<-t(covariate.model[(nr.psi0+1):nr.cov,])
    psi1<-saemix.model["psi0"][nr.psi0,]
    for(j in (nr.psi0+1):nr.cov)
      saemix.model["psi0"]<-rbind(saemix.model["psi0"],psi1)
    nr.psi0<-dim(saemix.model["psi0"])[1]
  }
  if(nr.psi0>1) fixed.ini[2:nr.psi0,]<-saemix.model["psi0"][2:nr.psi0,]
  
  #covariate.estim<-matrix(c(rep(saemix.model$fixed.estim,nr.psi0),t1),byrow=TRUE, nrow=nr.cov)
  #	covariate.estim<-matrix(rep(saemix.model["fixed.estim"],nr.psi0),byrow=TRUE, ncol=length(saemix.model["fixed.estim"]))
  #if(!is.null(dim(t1))) covariate.estim<-rbind(covariate.estim,t1) 
  #	covariate.estim<-covariate.estim*saemix.model["betaest.model"]
  # 29/05/2020 - changing definition of covariate.estim to 
  covariate.estim<-saemix.model["betaest.model"]
  covariate.estim[1,]<-saemix.model["fixed.estim"]
  
  betas.ini<-fixed.ini[which(saemix.model["betaest.model"]>0)]
  betas.ini<-matrix(betas.ini,ncol=1)
  
  nb.betas<-sum(saemix.model["betaest.model"])
  ind.covariate<-which(saemix.model["betaest.model"]==1)
  #matrix(which(covariate.model==1),nrow=1)
  
  # the covariates
  LCOV<-MCOV<-matrix(data=0,nrow=nb.betas,ncol=nb.parameters)
  j1<-1
  COV<-matrix(nrow=dim(Mcovariates)[1],ncol=0)
  pfix<-matrix(data=0,nrow=1,ncol=nb.parameters)
  mean.phi<-matrix(data=0,nrow=N,ncol=nb.parameters)
  for(j in 1:nb.parameters) {
    jcov<-which(saemix.model["betaest.model"][,j]==1)
    lambdaj<-fixed.ini[jcov,j]
    aj<-as.matrix(Mcovariates[,jcov])
    COV<-cbind(COV,aj)
    nlj<-length(lambdaj)
    j2<-j1+nlj-1
    LCOV[j1:j2,j]<-matrix(data=1,nrow=nlj,ncol=1)
    j1<-j2+1
    if(length(jcov)<=1) mean.phi[,j]<-aj*lambdaj else mean.phi[,j]<-aj%*%lambdaj
    pfix[j]<-length(lambdaj)
  }
  
  
  if(nb.parameters>1){
    indx.betaI<-cumsum(c(0,pfix[1:(nb.parameters-1)]))+1	
  } else{
    indx.betaI<-1
  }
  
  
  
  idx<-1:nb.betas
  indx.betaC<-idx[is.na(match(idx,indx.betaI))]
  saemix.model["indx.fix"]<-indx.betaI
  saemix.model["indx.cov"]<-indx.betaC
  
  COV2<-t(COV)%*%COV
  j.covariate<-which(LCOV==1)
  MCOV[j.covariate]<-betas.ini
  betas<-betas.ini
  
  ind.fix1<-which(covariate.estim[ind.covariate]==1)
  ind.fix0<-which(covariate.estim[ind.covariate]==0)
  COV1<-COV[,ind.fix1]
  #if(length(ind.fix0)==1) dstatCOV<-matrix(COV[,ind.fix0],ncol=1)%*%MCOV[ind.fix0,] else 
  dstatCOV<-COV[,ind.fix0,drop=FALSE]%*%MCOV[ind.fix0,]
  
  covariate.estim1<-covariate.estim
  covariate.estim1[,i0.omega2]<-0
  ind.fix11<-which(covariate.estim1[ind.covariate]==1)
  covariate.estim0<-covariate.estim
  covariate.estim0[,i1.omega2]<-0
  ind.fix10<-which(covariate.estim0[ind.covariate]==1)
  MCOV0<-MCOV[ind.fix10,i0.omega2,drop=FALSE]
  #if(is.null(dim(MCOV0)) & length(MCOV0)>0) MCOV0<-matrix(MCOV0,ncol=1)
  COV0<-COV[,ind.fix10]
  j0.covariate<-which(LCOV[ind.fix10,i0.omega2]==1)
  flag.fmin<-as.integer(sum(covariate.estim0[1,])>0)
  
  # using several Markov chains
  chdat<-new(Class="SaemixRepData",data=saemix.data, nb.chains=saemix.options$nb.chains)
  NM<-chdat["NM"]
  IdM<-chdat["dataM"]$IdM
  yM<-chdat["dataM"]$yM
  XM<-chdat["dataM"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
  io<-matrix(data=0,nrow=N,ncol=max(saemix.data["nind.obs"]))
  for(i in 1:N)
    io[i,1:saemix.data["nind.obs"][i]]<-1
  ioM<-do.call(rbind,rep(list(io),saemix.options$nb.chains))
  ind.ioM <- which(t(ioM)!=0)
  DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])
  
  # Initialisation of phiM
  if(length(i0.omega2)>0) {
    xmat<-covariate.estim[,i0.omega2]
    if(is.null(dim(xmat))) xmat<-matrix(xmat,ncol=length(i0.omega2))
    i0.temp<-which(colSums(xmat)==0)
    ind0.eta<-i0.omega2[i0.temp] # ind0.eta: index of parameters without IIV
  } else ind0.eta<-c()
  if(length(ind0.eta)>0) { # ind.eta: index of parameters with IIV
    idx<-1:nb.parameters
    ind.eta<-idx[-c(ind0.eta)]
  } else ind.eta<-1:nb.parameters
  nb.etas<-length(ind.eta)
  
  itest.phi<-1:NM
  ltest.phi<-length(itest.phi)
  phiM<-matrix(data=0,nrow=NM,ncol=nb.parameters)
  etaM<-matrix(data=0,nrow=NM,ncol=nb.etas)
  
  mean.phiM<-do.call(rbind,rep(list(mean.phi),saemix.options$nb.chains))
  kt<-0
  omega<-saemix.model["omega.init"]
  chol.omega<-try(chol(omega[ind.eta,ind.eta]),silent=TRUE)
  if(inherits(chol.omega,"try-error")) {
    #	cat("ind.eta=",ind.eta,"\n")
    #	print(saemix.model["omega.init"])
    #	print(omega[ind.eta,ind.eta])
    chol.omega<-saemix.model["omega.init"][ind.eta,ind.eta]<-omega[ind.eta, ind.eta]<-mydiag(nrow=length(ind.eta),ncol=length(ind.eta))
    if(saemix.options$warnings) cat("Problem inverting covariance matrix, setting initial Omega to diagonal.\n")
  }
  
  # Find a valid set of parameters wrt to the structural.model.
  # Any parameter set that does not generate NaN, inf or imaginary numbers
  # will satisfy this criteria.
  phiMc<-mean.phiM
  while (ltest.phi>0) {
    kt<-kt+1
    if (kt==100) 
      stop("stats:fit.saemix:FailedInitialParameterGuess\nFailed to find a valid initial parameter guess\n")
    end   
    etaMc<-0.5*matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
    phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
    etaM[itest.phi,]<-etaMc[itest.phi,]
    phiM[itest.phi,]<-phiMc[itest.phi,]
    psiM<-transphi(phiM,saemix.model["transform.par"])
    fpred<-structural.model(psiM, IdM, XM)
    inan<-(is.na(fpred)+is.infinite(fpred)+(Im(fpred)!=0))
    itest.phi<-unique(IdM[inan])
    ltest.phi<-length(itest.phi)
  }
  
  if(saemix.model["modeltype"]=="structural"){
    var.eta<-mydiag(saemix.model["omega.init"])
    theta0<-c(fixedpsi.ini,var.eta[i1.omega2],pres[saemix.model["indx.res"]])
    l1<-betas.ini
    l1[indx.betaI]<-transphi(matrix(l1[indx.betaI],nrow=1),saemix.model["transform.par"])
    allpar0<-c(l1,var.eta[i1.omega2],pres[saemix.model["indx.res"]]) # Eco modifié
  } else {
    var.eta<-mydiag(saemix.model["omega.init"])
    theta0<-c(fixedpsi.ini,var.eta[i1.omega2])
    l1<-betas.ini
    l1[indx.betaI]<-transphi(matrix(l1[indx.betaI],nrow=1),saemix.model["transform.par"])
    allpar0<-c(l1,var.eta[i1.omega2])
  }
  # Data - passed on to functions, unchanged
  Dargs<-list(IdM=IdM, XM=XM, yM=yM, NM=NM, N=N, nobs=saemix.data["ntot.obs"],
              yobs=saemix.data["data"][,saemix.data["name.response"]],transform.par=saemix.model["transform.par"],
              error.model=saemix.model["error.model"],structural.model=structural.model , etype.exp=which(saemix.model["error.model"] == "exponential"),modeltype=saemix.model["modeltype"])
  
  # List of indices and variables (fixed) - passed on to functions, unchanged
  nb.parest<-sum(covariate.estim)+ sum(saemix.model["covariance.model"][upper.tri(saemix.model["covariance.model"], diag=TRUE)])+n_mark+sum(as.integer(saemix.model["error.model"]=="combined"))
  
  Uargs<-list(nchains=saemix.options$nb.chains,nb.parameters=nb.parameters, nb.betas=nb.betas, nb.etas=nb.etas, 
              nb.parest=nb.parest,indx.betaC=indx.betaC, indx.betaI=indx.betaI, ind.res=ind.res,indest.omega=indest.omega,
              i0.omega2=i0.omega2, i1.omega2=i1.omega2,j.covariate=j.covariate, j0.covariate=j0.covariate,
              ind.fix10=ind.fix10, ind.fix11=ind.fix11, ind.fix1=ind.fix1, ind.fix0=ind.fix0,
              MCOV0=MCOV0, COV=COV, COV0=COV0, COV1=COV1, LCOV=LCOV, COV2=COV2, dstatCOV=dstatCOV, 
              Mcovariates=Mcovariates, ind.ioM=ind.ioM)
  # Variability-related elements
  omega.eta<-omega[ind.eta,ind.eta] # IIV matrix for estimated parameters
  varList<-list()
  if(saemix.model["modeltype"]=="structural") varList$pres<-pres
#  varList$indind0.eta<-ind0.eta # Eco modifié Alex
  varList$ind0.eta<-ind0.eta 
  varList$ind.eta<-ind.eta
  varList$omega<-omega
  varList$MCOV=MCOV
  varList$domega2<-do.call(cbind,rep(list((sqrt(mydiag(omega.eta)))*saemix.options$rw.ini),nb.etas))
  varList$diag.omega<-mydiag(omega)
  
  # List of options and settings (fixed) - passed on to functions, unchanged
  stepsize<-rep(1,saemix.options$nbiter.tot)
  stepsize[(saemix.options$nbiter.saemix[1]+1):saemix.options$nbiter.tot]<-1/
    (1:saemix.options$nbiter.saemix[2])
  stepsize[1:saemix.options$nbiter.burn]<-0
  opt<-list(stepsize.rw=saemix.options$stepsize.rw,stepsize=stepsize,
            proba.mcmc=saemix.options$proba.mcmc,nbiter.mcmc=saemix.options$nbiter.mcmc,
            nbiter.sa=saemix.options$nbiter.sa,nbiter.map=saemix.options$nbiter.map,alpha1.sa=saemix.options$alpha.sa,
            alpha0.sa=10^(-3/saemix.options$nbiter.sa),nbiter.saemix=saemix.options$nbiter.saemix,
            maxim.maxiter=saemix.options$maxim.maxiter,flag.fmin=flag.fmin)
  saemix.model["name.sigma"] = c("a.1","b.1","a.2","b.2")   # rajout alex pour deux modèles longi 
  return(list(saemix.model=saemix.model, Dargs=Dargs, Uargs=Uargs, varList=varList, opt=opt, DYF=DYF, phiM=phiM, mean.phi=mean.phi,betas=betas, fixedpsi.ini=fixedpsi.ini, allpar0=allpar0))
}

############################### Simulation - MCMC kernels (E-step) #############################

estep.alex<-function(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM) {
  # E-step - simulate unknown parameters
  # Input: kiter, Uargs, mean.phi (unchanged)
  # Output: varList, DYF, phiM (changed)
  # E-step - simulate unknown parameters
  # Input: kiter, Uargs, mean.phi (unchanged)
  # Output: varList, DYF, phiM (changed)
  
  # Function to perform MCMC simulation
  nb.etas<-length(varList$ind.eta)
  domega<-cutoff(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]),.Machine$double.eps)
  omega.eta<-varList$omega[varList$ind.eta,varList$ind.eta,drop=FALSE]
  omega.eta<-omega.eta-mydiag(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]))+mydiag(domega)
  chol.omega<-try(chol(omega.eta))
  somega<-solve(omega.eta)
  
  # "/" dans Matlab = division matricielle, selon la doc "roughly" B*INV(A) (et *= produit matriciel...)
  
  VK<-rep(c(1:nb.etas),2)
  mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))
  phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]
  
  U.y<-compute.LLy.alex(phiM,Uargs,Dargs,DYF,varList$pres, kiter=kiter)
  
  #  ind1<-c()
  etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]
  U.eta<-0.5*rowSums(etaM*(etaM%*%somega))
  if((kiter%%5)==1) print(apply(etaM,2, var))
  
  indchange1<-c()
  phiMc<-phiM
  for(u in 1:opt$nbiter.mcmc[1]) { # 1er noyau
    etaMc<-matrix(rnorm(Dargs$NM*nb.etas),ncol=nb.etas)%*%chol.omega
    phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
    Uc.y<-compute.LLy.alex(phiMc,Uargs,Dargs,DYF,varList$pres, kiter=kiter)
    Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
    deltau<-Uc.y-U.y+Uc.eta-U.eta
    ind<-which(deltau<(-1)*log(runif(Dargs$NM)))
    indchange1 <- c(indchange1, ind)
    #    ind1<-c(ind1,ind)
    etaM[ind,]<-etaMc[ind,]
    U.eta[ind]<-Uc.eta[ind]
    U.y[ind]<-Uc.y[ind]
  }
#  U.eta<-0.5*rowSums(etaM*(etaM%*%somega))
  if((kiter%%5)==1) print(apply(etaM,2, var))
  
  # Second stage
  
  indchange2<-c()
  if(opt$nbiter.mcmc[2]>0) {
    nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
    nrs2<-1
    for (u in 1:opt$nbiter.mcmc[2]) {
      for(vk2 in 1:nb.etas) {
        etaMc<-etaM
        #				cat('vk2=',vk2,' nrs2=',nrs2,"\n")
        etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(Dargs$NM*nrs2), ncol=nrs2)%*%mydiag(varList$domega2[vk2,nrs2],nrow=1) # 2e noyau ? ou 1er noyau+permutation?
        phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
        Uc.y<-compute.LLy.alex(phiMc,Uargs,Dargs,DYF,varList$pres, kiter=kiter)
        Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
        deltu<-Uc.y-U.y+Uc.eta-U.eta
        ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
        indchange2 <- c(indchange2, ind)
        #        ind1<-c(ind1,ind)
        etaM[ind,]<-etaMc[ind,]
        U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
        U.eta[ind]<-Uc.eta[ind]
        nbc2[vk2]<-nbc2[vk2]+length(ind)
        nt2[vk2]<-nt2[vk2]+Dargs$NM
      }
    }
    varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))
    if((kiter%%5)==1) print(apply(etaM,2, var))
  }
  
  indchange3<-c()
  if(opt$nbiter.mcmc[3]>0) {
    nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
    nrs2<-kiter%%(nb.etas-1)+2
    if(is.nan(nrs2)) nrs2<-1 # to deal with case nb.etas=1
    for (u in 1:opt$nbiter.mcmc[3]) {
      if(nrs2<nb.etas) {
        vk<-c(0,sample(c(1:(nb.etas-1)),nrs2-1))
        nb.iter2<-nb.etas
      } else {
        vk<-0:(nb.etas-1)
        #        if(nb.etas==1) vk<-c(0)
        nb.iter2<-1
      }
      for(k2 in 1:nb.iter2) {
        vk2<-VK[k2+vk]
        etaMc<-etaM
        etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(Dargs$NM*nrs2), ncol=nrs2)%*%mydiag(varList$domega2[vk2,nrs2])
        phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
        Uc.y<-compute.LLy.alex(phiMc,Uargs,Dargs,DYF,varList$pres, kiter=kiter)
        Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
        deltu<-Uc.y-U.y+Uc.eta-U.eta
        ind<-which(deltu<(-log(runif(Dargs$NM))))
        indchange3 <- c(indchange3, ind)
        #        ind1<-c(ind1,ind)
        etaM[ind,]<-etaMc[ind,]
        #        if(kiter<20 | (kiter>150 & kiter<170)) {
        #        	cat("kiter=",kiter,length(ind),"  varList$ind.eta=",varList$ind.eta,"  nrs2=",nrs2,"\n")
        #        	print(head(etaMc))
        #        }
        U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
        U.eta[ind]<-Uc.eta[ind]
        nbc2[vk2]<-nbc2[vk2]+length(ind)
        nt2[vk2]<-nt2[vk2]+Dargs$NM
      }
    }
    varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))
    if((kiter%%5)==1) print(apply(etaM,2, var))
  }
  
  
  if(opt$nbiter.mcmc[4]>0 & kiter<opt$nbiter.map) {
    etaMc<-etaM
    propc <- U.eta
    prop <- U.eta
    phi.map<-mean.phi
    i1.omega2<-varList$ind.eta
    iomega.phi1<-solve(omega.eta[i1.omega2,i1.omega2])
    
    # Setup for MAP calculation (MAP is identical no matter the chain)
    id<-Dargs$IdM[1:Dargs$nobs]
    xind<-Dargs$XM[1:Dargs$nobs,]
    yobs<-Dargs$yM[1:Dargs$nobs]
    id.list<-unique(id)
    
    if(Dargs$modeltype=="structural"){
      for(i in 1:length(id.list)) {
        isuj<-id.list[i]
        xi<-xind[id==isuj,,drop=FALSE]
        yi<-yobs[id==isuj]
        idi<-rep(1,length(yi))
        mean.phi1<-mean.phiM[i,i1.omega2]
        phii<-phiM[i,]
        phi1<-phii[i1.omega2]
        phi1.opti<-optim(par=phi1, fn=conditional.distribution_c, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=Dargs$transform.par, model=Dargs$structural.model, pres=varList$pres, err=Dargs$error.model)
        phi.map[i,i1.omega2]<-phi1.opti$par
      }
      
      # Repeat the map nchains time
      phi.map <- phi.map[rep(seq_len(nrow(phi.map)),Uargs$nchains ), ]
      
      map.psi<-transphi(phi.map,Dargs$transform.par)
      map.psi<-data.frame(id=id.list,map.psi)
      map.phi<-data.frame(id=id.list,phi.map)
      psi_map <- as.matrix(map.psi[,-c(1)])
      phi_map <- as.matrix(map.phi[,-c(1)])
      eta_map <- phi_map - mean.phiM
      
      fpred1<-Dargs$structural.model(psi_map, Dargs$IdM, Dargs$XM)
      gradf <- matrix(0L, nrow = length(fpred1), ncol = nb.etas) 
      
      ## Compute gradient of structural model (gradf)
      for (j in 1:nb.etas) {
        psi_map2 <- psi_map
        psi_map2[,j] <- psi_map[,j]+psi_map[,j]/1000
        fpred1<-Dargs$structural.model(psi_map, Dargs$IdM, Dargs$XM)
        fpred2<-Dargs$structural.model(psi_map2, Dargs$IdM, Dargs$XM)
        for (i in 1:(Dargs$NM)){
          r = which(Dargs$IdM==i)
          gradf[r,j] <- (fpred2[r] - fpred1[r])/(psi_map[i,j]/1000)
        }
      }
      
      
      ## Compute gradient of mapping (psi to phi) function (gradh)
      gradh <- list(omega.eta,omega.eta)
      for (i in 1:Dargs$NM){
        gradh[[i]] <- gradh[[1]]
      }
      for (j in 1:nb.etas) {
        phi_map2 <- phi_map
        phi_map2[,j] <- phi_map[,j]+phi_map[,j]/1000
        psi_map2 <- transphi(phi_map2,Dargs$transform.par) 
        for (i in 1:(Dargs$NM)){
          gradh[[i]][,j] <- (psi_map2[i,] - psi_map[i,])/(phi_map[i,]/1000)
        }
      }
      
      ## Calculation of the covariance matrix of the proposal
      Gamma <- chol.Gamma <- inv.chol.Gamma <- inv.Gamma <- list(omega.eta,omega.eta)
      for (i in 1:(Dargs$NM)){
        r = which(Dargs$IdM==i)
        temp <- gradf[r,]%*%gradh[[i]]
        Gamma[[i]] <- solve(t(temp)%*%temp/(varList$pres[1])^2+solve(omega.eta)) # Eco: why pres[1] ?
        chol.Gamma[[i]] <- chol(Gamma[[i]])
        inv.chol.Gamma[[i]] <- solve(chol.Gamma[[i]])
        inv.Gamma[[i]] <- solve(Gamma[[i]])
      }
      
    } else {
      for(i in 1:length(id.list)) {
        isuj<-id.list[i]
        xi<-xind[id==isuj,,drop=FALSE]
        yi<-yobs[id==isuj]
        idi<-rep(1,length(yi))
        mean.phi1<-mean.phiM[i,i1.omega2]
        phii<-phiM[i,]
        phi1<-phii[i1.omega2]
        phi1.opti<-optim(par=phi1, fn=conditional.distribution_d, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=Dargs$transform.par, model=Dargs$structural.model)
        phi.map[i,i1.omega2]<-phi1.opti$par
      }
      #rep the map nchains time
      phi.map <- phi.map[rep(seq_len(nrow(phi.map)),Uargs$nchains ), ] 
      map.psi<-transphi(phi.map,Dargs$transform.par)
      map.psi<-data.frame(id=id.list,map.psi)
      map.phi<-data.frame(id=id.list,phi.map)
      
      psi_map <- as.matrix(map.psi[,-c(1)])
      phi_map <- as.matrix(map.phi[,-c(1)])
      eta_map <- phi_map[,varList$ind.eta] - mean.phiM[,varList$ind.eta]
      
      #gradient at the map estimation
      gradp <- matrix(0L, nrow = Dargs$NM, ncol = nb.etas) 
      
      for (j in 1:nb.etas) {
        phi_map2 <- phi_map
        phi_map2[,j] <- phi_map[,j]+phi_map[,j]/100;
        psi_map2 <- transphi(phi_map2,Dargs$transform.par) 
        fpred1<-Dargs$structural.model(psi_map, Dargs$IdM, Dargs$XM)
        DYF[Uargs$ind.ioM]<- fpred1
        l1<-colSums(DYF)
        fpred2<-Dargs$structural.model(psi_map2, Dargs$IdM, Dargs$XM)
        DYF[Uargs$ind.ioM]<- fpred2
        l2<-colSums(DYF)
        
        for (i in 1:(Dargs$NM)){
          gradp[i,j] <- (l2[i] - l1[i])/(phi_map[i,j]/100)
        }
      }
      
      #calculation of the covariance matrix of the proposal
      fpred<-Dargs$structural.model(psi_map, Dargs$IdM, Dargs$XM)
      DYF[Uargs$ind.ioM]<- fpred
      denom <- colSums(DYF)
      
      Gamma <- chol.Gamma <- inv.Gamma <- list(omega.eta,omega.eta)
      z <- matrix(0L, nrow = length(fpred), ncol = 1) 
      for (i in 1:(Dargs$NM)){
        Gamma[[i]] <- solve(gradp[i,]%*%t(gradp[i,])/denom[i]^2+solve(omega.eta))
        chol.Gamma[[i]] <- chol(Gamma[[i]])
        inv.Gamma[[i]] <- solve(Gamma[[i]])
      }
      
    }
    
    etaM <- eta_map
    phiM<-etaM+mean.phiM
    U.eta<-0.5*rowSums(etaM*(etaM%*%somega))
    U.y<-compute.LLy.alex(phiM,Uargs,Dargs,DYF,varList$pres, kiter=kiter)
    
    for (u in 1:opt$nbiter.mcmc[4]) {
      
      #generate candidate eta with new proposal
      for (i in 1:(Dargs$NM)){
        Mi <- rnorm(nb.etas)%*%chol.Gamma[[i]]
        etaMc[i,varList$ind.eta]<- eta_map[i,varList$ind.eta] + Mi
      }
      
      phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc[,varList$ind.eta]
      Uc.y<-compute.LLy.alex(phiM,Uargs,Dargs,DYF,varList$pres, kiter=kiter)
      Uc.eta<-0.5*rowSums(etaMc[,varList$ind.eta]*(etaMc[,varList$ind.eta]%*%somega))
      
      for (i in 1:(Dargs$NM)){
        propc[i] <- 0.5*rowSums((etaMc[i,varList$ind.eta]-eta_map[i,varList$ind.eta])*(etaMc[i,varList$ind.eta]-eta_map[i,varList$ind.eta])%*%inv.Gamma[[i]])
        prop[i] <- 0.5*rowSums((etaM[i,varList$ind.eta]-eta_map[i,varList$ind.eta])*(etaM[i,varList$ind.eta]-eta_map[i,varList$ind.eta])%*%inv.Gamma[[i]])
      }
      
      deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
      ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
      etaM[ind,varList$ind.eta]<-etaMc[ind,varList$ind.eta]
      U.y[ind]<-Uc.y[ind]
      U.eta[ind]<-Uc.eta[ind]
      
    }
  }
  if((kiter%%5)==1) cat("nmod=",length(unique(indchange1)),"/",length(unique(indchange2)),"/",length(unique(indchange3))," --", exp(mean.phiM[1,]),"\n")
  
  phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
  return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM))
  
}



################## Stochastic approximation - compute sufficient statistics (M-step) #####################

mstep.alex<-function(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat) {
  # M-step - stochastic approximation
  # Input: kiter, Uargs, structural.model, DYF, phiM (unchanged)
  # Output: varList, phi, betas, suffStat (changed)
  #					mean.phi (created)
  
  # Update variances - TODO - check if here or elsewhere
  nb.etas<-length(varList$ind.eta)
  domega<-cutoff(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]),.Machine$double.eps)
  omega.eta<-varList$omega[varList$ind.eta,varList$ind.eta,drop=FALSE]
  omega.eta<-omega.eta-mydiag(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]))+mydiag(domega)
  #  print(varList$omega.eta)
  chol.omega<-try(chol(omega.eta))
  d1.omega<-Uargs$LCOV[,varList$ind.eta]%*%solve(omega.eta)
  d2.omega<-d1.omega%*%t(Uargs$LCOV[,varList$ind.eta])
  comega<-Uargs$COV2*d2.omega
  
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  #	if(Dargs$error.model=="exponential")
  #		fpred<-log(cutoff(fpred))
  ff<-matrix(fpred,nrow=Dargs$nobs,ncol=Uargs$nchains)
  for(k in 1:Uargs$nchains) phi[,,k]<-phiM[((k-1)*Dargs$N+1):(k*Dargs$N),]
  # overall speed similar
  #    phi<-aperm(array(phiM,c(N,nchains,3)),c(1,3,2))
  stat1<-apply(phi[,varList$ind.eta,,drop=FALSE],c(1,2),sum) # sum on columns ind.eta of phi, across 3rd dimension
  stat2<-matrix(data=0,nrow=nb.etas,ncol=nb.etas)
  stat3<-apply(phi**2,c(1,2),sum) #  sum on phi**2, across 3rd dimension
  statr1<-0
  statr2<-0
  for(k in 1:Uargs$nchains) {
    phik<-phi[,varList$ind.eta,k]
    stat2<-stat2+t(phik)%*%phik
    fk<-ff[,k]
    if(length(Dargs$error.model)==2) {    # à généraliser !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ici 2 réponses 
      if(!is.na(match(Dargs$error.model[1],c("constant","exponential")))){
        resk1<-sum((Dargs$yobs[Dargs$XM$ytype==1]-fk[Dargs$XM$ytype==1])**2) } else {
          if(Dargs$error.model[1]=="proportional") {
          resk1<-sum((Dargs$yobs[Dargs$XM$ytype==1]-fk[Dargs$XM$ytype==1])**2/cutoff(fk[Dargs$XM$ytype==1]**2,.Machine$double.eps))
          } else resk<-0
        }
      if(!is.na(match(Dargs$error.model[2],c("constant","exponential")))){
        resk2<-sum((Dargs$yobs[Dargs$XM$ytype==2]-fk[Dargs$XM$ytype==2])**2)} else {
          if(Dargs$error.model[2]=="proportional") {
            #		        idx.okpred<-which(fk>.Machine$double.eps)
            #		        vec<-(Dargs$yobs-fk)**2/cutoff(fk**2,.Machine$double.eps)
            #		        resk<-sum(vec[idx.okpred])
            resk2<-sum((Dargs$yobs[Dargs$XM$ytype==2]-fk[Dargs$XM$ytype==2])**2/cutoff(fk[Dargs$XM$ytype==2]**2,.Machine$double.eps))
          } else resk<-0
        }
    } else resk<-0
    statr1<-statr1+resk1
    statr2<-statr2+resk2
  }
  # Update sufficient statistics
  suffStat$statphi1<-suffStat$statphi1+opt$stepsize[kiter]*(stat1/Uargs$nchains-suffStat$statphi1)
  suffStat$statphi2<-suffStat$statphi2+opt$stepsize[kiter]*(stat2/Uargs$nchains-suffStat$statphi2)
  suffStat$statphi3<-suffStat$statphi3+opt$stepsize[kiter]*(stat3/Uargs$nchains-suffStat$statphi3)
  suffStat$statrese1<-suffStat$statrese1+opt$stepsize[kiter]*(statr1/Uargs$nchains-suffStat$statrese1)
  suffStat$statrese2<-suffStat$statrese2+opt$stepsize[kiter]*(statr2/Uargs$nchains-suffStat$statrese2)
  
  ############# Maximisation
  ##### fixed effects
  
  if (opt$flag.fmin && kiter>=opt$nbiter.sa) {
    temp<-d1.omega[Uargs$ind.fix11,]*(t(Uargs$COV1)%*%(suffStat$statphi1-Uargs$dstatCOV[,varList$ind.eta]))
    betas[Uargs$ind.fix11]<-solve(comega[Uargs$ind.fix11,Uargs$ind.fix11],rowSums(temp)) 
    # ECO TODO: utiliser optimise dans le cas de la dimension 1
    #		if(length(Uargs$ind.fix10)>1) 
    beta0<-optim(par=betas[Uargs$ind.fix10],fn=compute.Uy.alex,phiM=phiM,pres=varList$pres,args=Uargs,Dargs=Dargs,DYF=DYF,control=list(maxit=opt$maxim.maxiter))$par # else
    #		beta0<-optimize(f=compute.Uy, interval=c(0.01,100)*betas[Uargs$ind.fix10],phiM=phiM,pres=varList$pres,args=Uargs,Dargs=Dargs,DYF=DYF)
    #		if(kiter==opt$nbiter.sa) {
    #		  cat("ind.fix10=",Uargs$ind.fix10,"ind.fix11=",Uargs$ind.fix11,"ind.fix1=",Uargs$ind.fix1,"ind.fix0=",Uargs$ind.fix0,"\n")
    #  		cat(betas,"\n")
    #		}
    betas[Uargs$ind.fix10]<-betas[Uargs$ind.fix10]+opt$stepsize[kiter]*(beta0-betas[Uargs$ind.fix10])
  } else {
    temp<-d1.omega[Uargs$ind.fix1,]*(t(Uargs$COV1)%*%(suffStat$statphi1-Uargs$dstatCOV[,varList$ind.eta]))
    betas[Uargs$ind.fix1]<-solve(comega[Uargs$ind.fix1,Uargs$ind.fix1],rowSums(temp)) 
  }
  
  varList$MCOV[Uargs$j.covariate]<-betas
  mean.phi<-Uargs$COV %*% varList$MCOV
  e1.phi<-mean.phi[,varList$ind.eta,drop=FALSE]
  
  # Covariance of the random effects
  omega.full<-matrix(data=0,nrow=Uargs$nb.parameters,ncol=Uargs$nb.parameters)
  omega.full[varList$ind.eta,varList$ind.eta]<-suffStat$statphi2/Dargs$N + t(e1.phi)%*%e1.phi/Dargs$N - t(suffStat$statphi1)%*%e1.phi/Dargs$N - t(e1.phi)%*%suffStat$statphi1/Dargs$N
  ## rajout alex : 
  #omega.full[1,1]=saemix.model["omega.init"][1,1]
  varList$omega[Uargs$indest.omega]<-omega.full[Uargs$indest.omega]
  
  # Simulated annealing (applied to the diagonal elements of omega)
  if (kiter<=opt$nbiter.sa) {
    diag.omega.full<-mydiag(omega.full)
    vec1<-diag.omega.full[Uargs$i1.omega2]
    vec2<-varList$diag.omega[Uargs$i1.omega2]*opt$alpha1.sa
    idx<-as.integer(vec1<vec2)
    varList$diag.omega[Uargs$i1.omega2]<-idx*vec2+(1-idx)*vec1
    varList$diag.omega[Uargs$i0.omega2]<-varList$diag.omega[Uargs$i0.omega2]*opt$alpha0.sa
  } else {
    varList$diag.omega<-mydiag(varList$omega)
  }
  varList$omega<-varList$omega-mydiag(mydiag(varList$omega))+mydiag(varList$diag.omega)
  
  # Residual error
  n_mark=2
  # Modified to add SA to constant and exponential residual error models (Edouard Ollier 10/11/2016)
  if(Dargs$modeltype=="structural") {
    if(length(Uargs$ind.res)/n_mark==1) { # 
      if (Dargs$error.model[1] %in% c("constant","exponential")) {  # a généraliser  
        sig21<-suffStat$statrese1/length(Dargs$XM$ytype[Dargs$XM$ytype==1])
        sig22<-suffStat$statrese2/length(Dargs$XM$ytype[Dargs$XM$ytype==2])
        if (kiter<=opt$nbiter.sa) {
          varList$pres[1]<-max(varList$pres[1]*opt$alpha1.sa,sqrt(sig21))
          varList$pres[3]<-max(varList$pres[3]*opt$alpha1.sa,sqrt(sig22)) # idem a généraliser 
        } else {
          varList$pres[1]<-sqrt(sig21)
          varList$pres[3]<-sqrt(sig22)
        }
      }
      if (Dargs$error.model[1]=="proportional") {    ## a faire cette partie, ici pas du tout adaptée à 2longi
        sig21<-suffStat$statrese1/length(Dargs$XM$ytype[Dargs$XM$ytype==1]) 
        sig22<-suffStat$statrese2/length(Dargs$XM$ytype[Dargs$XM$ytype==2]) 
        if (kiter<=opt$nbiter.sa) {
          varList$pres[2]<-max(varList$pres[2]*opt$alpha1.sa,sqrt(sig21))
          varList$pres[4]<-max(varList$pres[4]*opt$alpha1.sa,sqrt(sig22))
        } else {
          varList$pres[2]<-sqrt(sig21)
          varList$pres[4]<-sqrt(sig22)
        }
      }
      
    } else {
      #	if (Dargs$error.model=="combined") {
      # ECO TODO: check and secure (when fpred<0 => NaN, & what happens if bres<0 ???)
      ABres<-optim(par=varList$pres,fn=ssq,y=Dargs$yM,f=fpred,etype=Dargs$XM$ytype)$par
      if (kiter<=opt$nbiter.sa) {
        for(i in 1:length(varList$pres)) varList$pres[i]<-max(varList$pres[i]*opt$alpha1.sa,ABres[i])
      }  else {
        if (kiter<=opt$nbiter.saemix[1]) {
          for(i in 1:length(varList$pres)) varList$pres[i]<-ABres[i]
        } else {
          for(i in 1:length(varList$pres)) varList$pres[i]<-varList$pres[i]+opt$stepsize[kiter]*(ABres[i]-varList$pres[i])
        }
      }
    }
  }
  return(list(varList=varList,mean.phi=mean.phi,phi=phi,betas=betas,suffStat=suffStat))
}

################## Likelihood ##################
compute.LLy.alex = function(phiM, args, Dargs, DYF, pres, kiter) {
  psiM <- transphi(phiM, Dargs$transform.par)
  fpred <- Dargs$structural.model(psiM, Dargs$IdM, Dargs$XM)
  for (ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype == ityp] <- log(cutoff(fpred[Dargs$XM$ytype == ityp]))
  
  gpred <- error(fpred, pres, Dargs$XM$ytype)
  
  DYF[args$ind.ioM[Dargs$XM$ytype==1]] <- 0.5 * ((Dargs$yM[Dargs$XM$ytype==1] - fpred[Dargs$XM$ytype==1])/gpred[Dargs$XM$ytype==1])^2 + log(gpred[Dargs$XM$ytype==1])
  DYF[args$ind.ioM[Dargs$XM$ytype==2]] <- 0.5 * ((Dargs$yM[Dargs$XM$ytype==2] - fpred[Dargs$XM$ytype==2])/gpred[Dargs$XM$ytype==2])^2 + log(gpred[Dargs$XM$ytype==2])
  
  #DYF[Uargs$ind.ioM] <- 0.5 * ((Dargs$yM - fpred)/gpred)^2 + log(gpred)
  # ca doit faire la même chose ça  # oui (si les réponses sont toutes Gaussiennes)
  # DYF1<-DYF
  # DYF1[Uargs$ind.ioM] <- 0.5 * ((Dargs$yM - fpred)/gpred)^2 + log(gpred)
  # summary(c(DYF-DYF1))
  
  U <- colSums(DYF)
  return(U)
}

compute.Uy.alex = function(b0, phiM, pres, args, Dargs, DYF, kiter) {
  args$MCOV0[args$j0.covariate] <- b0
  phi0 <- args$COV0 %*% args$MCOV0
  phiM[, args$i0.omega2] <- do.call(rbind, rep(list(phi0), 
                                               args$nchains))
  psiM <- transphi(phiM, Dargs$transform.par)
  #ytype = Dargs[["XM"]][["ytype"]]
  
  fpred <- Dargs$structural.model(psiM, Dargs$IdM, Dargs$XM)
  for (ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype == ityp] <- log(cutoff(fpred[Dargs$XM$ytype == ityp]))
  
  gpred <- error(fpred, pres, Dargs$XM$ytype)
  DYF[args$ind.ioM] <- 0.5 * ((Dargs$yM - fpred)/gpred)^2 + log(gpred)

  #DYF[args$ind.ioM[ytype==2]] <- -fpred[ytype==2]
  #cat(sum(DYF[args$ind.ioM[ytype==1]])," ", sum(DYF[args$ind.ioM[ytype==2]]), "\n")
  
  U <- sum(DYF)
  return(U)
}
