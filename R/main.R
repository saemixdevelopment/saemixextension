#' Stochastic Approximation Expectation Maximization (SAEM) algorithm
#' 
#' SAEM algorithm perform parameter estimation for nonlinear mixed effects
#' models without any approximation of the model (linearization, quadrature
#' approximation, . . . )
#' 
#' The SAEM algorithm is a stochastic approximation version of the standard EM
#' algorithm proposed by Kuhn and Lavielle (see reference). Details of the
#' algorithm can be found in the pdf file accompanying the package.
#' 
#' @name saemix
#' @aliases initialiseMainAlgo mstep estep
#' 
#' @param model an object of class SaemixModel, created by a call to the
#' function \code{\link{saemixModel}}
#' @param data an object of class SaemixData, created by a call to the function
#' \code{\link{saemixData}}
#' @param control a list of options, see \code{\link{saemixControl}}
#' @return An object of class SaemixObject containing the results of the fit of
#' the data by the non-linear mixed effect model. A summary of the results is
#' printed out to the terminal, and, provided the appropriate options have not
#' been changed, numerical and graphical outputs are saved in a directory.
#' 
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{SaemixObject}}, \code{\link{saemixControl}},
#' \code{\link{plot.saemix}}
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' @author Audrey Lavenu
#' @author Marc Lavielle
#' 
#' @keywords models
#' @examples
#' 
#' data(theo.saemix)
#' 
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
#'    name.group=c("Id"),name.predictors=c("Dose","Time"),
#'    name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'    units=list(x="hr",y="mg/L", covariates=c("kg","-")), name.X="Time")
#' 
#' model1cpt<-function(psi,id,xidep) { 
#' 	  dose<-xidep[,1]
#' 	  tim<-xidep[,2]  
#' 	  ka<-psi[id,1]
#' 	  V<-psi[id,2]
#' 	  CL<-psi[id,3]
#' 	  k<-CL/V
#' 	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#' 	  return(ypred)
#' }
#' 
#' saemix.model<-saemixModel(model=model1cpt,modeltype="structural",
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' \donttest{
#' # Not run (strict time constraints for CRAN)
#'  saemix.fit<-saemix(saemix.model,saemix.data,list(seed=632545,directory="newtheo",
#'  save=FALSE,save.graphs=FALSE, print=FALSE))
#'  
#' # Prints a summary of the results
#' print(saemix.fit)
#' 
#' # Outputs the estimates of individual parameters
#' psi(saemix.fit)
#' 
#' # Shows some diagnostic plots to evaluate the fit
#' plot(saemix.fit)
#' }
#' @export saemix

saemix<-function(model,data,control=list()) {

# Convergence plots during fit (special function, not user-level)
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
  if(!is(model,"SaemixModel")) {
    message("Please provide a valid model object (see the help page for SaemixModel)\n")
    return("Missing model")
  }
  if(!is(data,"SaemixData")) {
    message("Please provide a valid data object (see the help page for SaemixData)\n")
    return("Missing data")
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
  OLDRAND<-TRUE
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
  
# Running the algorithm
# hw=waitbar(1,'Estimating the population parameters (SAEM). Wait...');
  if(saemix.options$displayProgress) par(ask=FALSE)
  if(saemix.options$warnings) cat("Running main SAEM algorithm\n")
  if(saemix.options$warnings) print(date())
  for (kiter in 1:saemix.options$nbiter.tot) { # Iterative portion of algorithm
    
# SAEM convergence plots
	if(kiter%%saemix.options$nbdisplay==0) {
	  if(saemix.options$warnings) cat(".")
    if(saemix.options$displayProgress)    
      try(convplot.infit(allpar,saemix.options$nbiter.saemix[1],niter=(kiter-2)))
  }
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
  
  etaM<-xmcmc$etaM # only need etaM here (re-created in estep otherwise)
  if(saemix.options$warnings) cat("\n    Minimisation finished\n")
  if(saemix.options$warnings) print(date())
  
  ############# After end of iterations
  fixed.effects<-0*betas
  fixed.effects[Uargs$indx.betaI]<-fixed.psi
  fixed.effects[Uargs$indx.betaC]<-betaC
  varList$omega[Uargs$i0.omega2,]<-0
  varList$omega[,Uargs$i0.omega2]<-0
  
  ##### Compute the individual parameters (MAP)
  phi[,Uargs$i0.omega2,1:saemix.options$nb.chains]<-mean.phi[,Uargs$i0.omega2]
  phi.samp<-phi
  phi<-apply(phi,c(1,2),mean)
  
  ##### Conditional means and variances used for the estimation of the log-likelihood via Importance Sampling
  cond.mean.phi<-phi
  sphi1<-phi
  sphi1[,varList$ind.eta]<-suffStat$statphi1
  cond.mean.phi[,Uargs$i1.omega2]<-sphi1[,Uargs$i1.omega2]
  cond.var.phi<-array(data=0,dim=dim(phi))
  cond.var.phi[,Uargs$i1.omega2]<-suffStat$statphi3[,Uargs$i1.omega2]-cond.mean.phi[,Uargs$i1.omega2]**2
  cond.mean.psi<-transphi(cond.mean.phi,saemixObject["model"]["transform.par"])
  
  cond.mean.eta<-matrix(0,nrow=dim(etaM)[1],ncol=Uargs$nb.parameters)
  cond.mean.eta[,varList$ind.eta]<-etaM
  cond.mean.eta<-array(t(cond.mean.eta),dim=c(Uargs$nb.parameters, Dargs$N, saemix.options$nb.chains))
  cond.mean.eta<-t(apply(cond.mean.eta,c(1,2),mean))
  
# Updating objects
  saemix.model["Mcovariates"]<-Uargs$Mcovariates
  saemix.model["indx.res"]<-Uargs$ind.res
  saemix.model["indx.fix"]<-Uargs$indx.betaI
  saemix.model["indx.cov"]<-Uargs$indx.betaC
  saemix.model["indx.omega"]<-Uargs$i1.omega2

# Filling in result object
  if(Dargs$modeltype=="structural") {
    saemix.res<-new(Class="SaemixRes",status="fitted",modeltype=Dargs$modeltype,name.fixed=saemix.model["name.fixed"], name.random=saemix.model["name.random"],name.sigma=saemix.model["name.sigma"], fixed.effects=c(fixed.effects),fixed.psi=c(fixed.psi),betas=betas,betaC=betaC, omega=varList$omega,respar=varList$pres,cond.mean.phi=cond.mean.phi,cond.var.phi=cond.var.phi, mean.phi=mean.phi, phi=phi,phi.samp=phi.samp,parpop=parpop,allpar=allpar,MCOV=varList$MCOV)
  } else{
    saemix.res<-new(Class="SaemixRes",status="fitted",modeltype=Dargs$modeltype,name.fixed=saemix.model["name.fixed"], name.random=saemix.model["name.random"],name.sigma=saemix.model["name.sigma"], fixed.effects=c(fixed.effects),fixed.psi=c(fixed.psi),betas=betas,betaC=betaC, omega=varList$omega,cond.mean.phi=cond.mean.phi,cond.var.phi=cond.var.phi, mean.phi=mean.phi, phi=phi,phi.samp=phi.samp,parpop=parpop,allpar=allpar,MCOV=varList$MCOV)
  }
  saemix.res["indx.res"]<-Uargs$ind.res
  saemix.res["indx.fix"]<-Uargs$indx.betaI
  saemix.res["indx.cov"]<-Uargs$indx.betaC
  saemix.res["indx.omega"]<-Uargs$i1.omega2
  saemix.res["npar.est"]<-Uargs$nb.parest
  saemix.res["nbeta.random"]<- sum(saemix.model["betaest.model"]%*%diag(saemix.model["fixed.estim"])%*%as.matrix(diag(saemix.model["covariance.model"])))
  saemix.res["nbeta.fixed"]<-  sum(saemix.model["betaest.model"]%*%diag(saemix.model["fixed.estim"])%*%as.matrix(-1*diag(saemix.model["covariance.model"])+1))
  saemix.res["cond.mean.psi"]<-cond.mean.psi
  saemix.res["cond.mean.eta"]<-cond.mean.eta

# Updating elements of saemixObject
  saemixObject["model"]<-saemix.model
  saemixObject["results"]<-saemix.res
  saemixObject["options"]<-saemix.options
#  saemixObject["rep.data"]<-chdat # Utile ? maybe remove rep.data

# ECO TODO check
# a la fin: mais verifier, pe pb de distribution ??? ie allpar sur l'echelle des betas et pas parpop ? a verifier
# saemix.res["allpar"]<-allpar
# saemix.res["parpop"]<-allpar[,-c(indx.betaC)]
#### Final computations
# Compute the MAP estimates of the PSI_i's 
  if(saemix.options$map) {
    x<-try(saemixObject<-map.saemix(saemixObject))
    if(inherits(x,"try-error") & saemixObject@options$warnings) message("Problem estimating the MAP parameters\n") else 
      if(!saemix.options$save.graphs) saemixObject<-saemix.predict(saemixObject, type=c("ipred","ppred")) # if no graphs, compute predictions all the same
  }
  
# Compute the Fisher Information Matrix & update saemix.res
  if(saemix.options$fim) {
    x<-try(saemixObject<-fim.saemix(saemixObject))
    if(inherits(x,"try-error") & saemixObject@options$warnings) message("Problem estimating the FIM\n")
  }
  
# Estimate the log-likelihood via importance Sampling/Gaussian quadrature
  if(saemix.options$ll.is) {
    x<-try(saemixObject<-llis.saemix(saemixObject))
    if(inherits(x,"try-error") & saemixObject@options$warnings) message("Problem estimating the likelihood by IS\n")
  }
  if(saemix.options$ll.gq) {
    x<-try(saemixObject<-llgq.saemix(saemixObject))
    if(inherits(x,"try-error") & saemixObject@options$warnings) message("Problem estimating the likelihood by GQ\n")
  }
  
#### Pretty printing the results (TODO finish in particular cov2cor)
  if(saemix.options$print) print(saemixObject,digits=2)

#### Save the results to a file
  if(saemix.options$save | saemix.options$save.graphs) {
# create directory to save the results
     if(saemix.options$directory!="") xsave<-dir.create(saemix.options$directory)
     if(!xsave) {
# Check that we're not trying to create a directory with the same name as a file
       if(!file_test("-d",saemix.options$directory)) {
         if(saemix.options$warnings) message("Unable to create directory",saemix.options$directory)
         saemix.options$directory<-"newdir"
         dir.create(saemix.options$directory)         
         xsave<-file_test("-d",saemix.options$directory)
         if(!xsave) {
           saemix.options$directory<-""
           xsave<-TRUE
           if(saemix.options$warnings) message(", saving in current directory.\n")
         } else {if(saemix.options$warnings) message(", saving results in newdir instead.\n")}
       } else {
       xsave<-TRUE
       if(saemix.options$warnings) message("Overwriting files in directory",saemix.options$directory,"\n")
       }
     }
   }
  if(saemix.options$save) {
    namres<-ifelse(saemix.options$directory=="","pop_parameters.txt", file.path(saemix.options$directory,"pop_parameters.txt"))
    xtry<-try(sink(namres))
    if(!inherits(xtry,"try-error")) {
    print(saemixObject)
    sink()
    namres<-ifelse(saemix.options$directory=="","indiv_parameters.txt", file.path(saemix.options$directory,"indiv_parameters.txt"))
    if(length(saemixObject["results"]["map.psi"])>0) {
      tab<-cbind(Id=unique(saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]),saemixObject["results"]["map.psi"])
      colnames(tab)[1]<-saemixObject["data"]["name.group"]
      write.table(tab,namres,quote=FALSE, row.names=FALSE)
    }
    } else {
       message("Unable to save results, check writing permissions and/or path to directory.\n")
     }
  }
# ECO TODO finish, adding all
  if(saemix.options$save.graphs) {
    saemixObject<-saemix.predict(saemixObject)
    if(saemix.options$directory=="") namgr<-"diagnostic_graphs.ps" else
      namgr<-file.path(saemix.options$directory,"diagnostic_graphs.ps")
    xtry<-try(postscript(namgr,horizontal=TRUE))
    if(!inherits(xtry,"try-error")) {
      oldpar <- par(no.readonly = TRUE)    # code line i
      on.exit(par(oldpar))            # code line i + 1 
      par(mfrow=c(1,1))
      try(plot(saemixObject,plot.type="data"))
      
      try(plot(saemixObject,plot.type="convergence"))
      
      if(length(saemixObject["results"]["ll.is"])>0) {
        par(mfrow=c(1,1))
        try(plot(saemixObject, plot.type="likelihood"))
      }
      
      try(plot(saemixObject,plot.type="observations.vs.predictions"))
      
      try(plot(saemixObject,plot.type="random.effects"))
      
      try(plot(saemixObject,plot.type="correlations"))

# Note: can replace all this by:
#    default.saemix.plots(saemixObject)

    dev.off()
    
    if(saemix.options$directory=="") namgr<-"individual_fits.ps" else
      namgr<-file.path(saemix.options$directory,"individual_fits.ps")
    postscript(namgr,horizontal=FALSE)
    try(plot(saemixObject,plot.type="individual.fit"))
    dev.off()
    } else {
       message("Unable to save results, check writing permissions and/or path to directory.\n")
     }
  }

  options(warn=opt.warn)
  return(saemixObject)
}
