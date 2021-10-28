############################### Initialising main algorithm #############################
initialiseMainAlgo<-function(saemix.data,saemix.model,saemix.options) {
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
#               nb.parest=nb.parest,indx.betaC=indx.betaC, indx.betaI=indx.betaI, ind.res=ind.res,
#               indest.omega=indest.omega, i0.omega2=i0.omega2, i1.omega2=i1.omega2,	j.covariate=j.covariate, 
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
		pres<-c(saemix.model["error.init"][1],saemix.model["error.init"][2])
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
	allpar0<-c(l1,var.eta[i1.omega2],pres[ind.res])
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
	nb.parest<-sum(covariate.estim)+ sum(saemix.model["covariance.model"][upper.tri(saemix.model["covariance.model"], diag=TRUE)])+1+as.integer(saemix.model["error.model"]=="combined")
	
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
	varList$indind0.eta<-ind0.eta
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
	
	return(list(saemix.model=saemix.model, Dargs=Dargs, Uargs=Uargs, varList=varList, opt=opt, DYF=DYF, phiM=phiM, mean.phi=mean.phi,betas=betas, fixedpsi.ini=fixedpsi.ini, allpar0=allpar0))
}
