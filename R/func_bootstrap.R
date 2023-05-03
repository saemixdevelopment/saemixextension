#########################################################################################
# Bootstrap saemix

#' Bootstrap for saemix fits
#' 
#' This function provides bootstrap estimates for a saemixObject run. Different bootstrap approaches have been 
#' implemented (see below for details), with the default method being the conditional non-parametric bootstrap ()
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' @param method the name of the bootstrap algorithm to use (one of: "case", "residual", "parametric" or "conditional") (defaults to "conditional")
#' @param nboot number of bootstrap samples
#' @param nsamp number of samples from the conditional distribution (for method="conditional")
#' @param saemix.options list of options to run the saemix algorithm. Defaults to the options in the object, but suppressing the estimation of individual parameters (map=FALSE) and likelihood (ll.is=FALSE), and the options to print out intermediate results and graphs (displayProgress=FALSE,save.graphs=FALSE,print=FALSE)
#' 
#' @details Different bootstrap algorithms have been proposed for non-linear mixed effect models, 
#' to account for the hierarchical nature of these models (see review in Thai et al. 2013, 2014)
#' In saemix, we implemented the following bootstrap algorithms, which can be selected using the "method" argument:
#' * case refers to case bootstrap, where resampling is performed at the level of the individual and the bootstrap sample consists in sampling with replacement from the observed individuals
#' * residual refers to non-parametric residual bootstrap, where for each individual the residuals for the random effects and for the residual error are resampled from the individual estimates after the fit
#' * parametric refers to parametric residual bootstrap, where these residuals are sampled from their theoretical distributions. In
#' saemix, random effects are drawn from a normal distribution using the estimated variance in saemixObject, and transformed 
#' to the distribution specified by the transform.par component of the model, and the residual error is drawn from a normal distribution
#' using the estimated residual variance.
#' * conditional refers to the conditional non-parametric bootstrap, where the random effects are resampled from the individual 
#'   conditional distributions as in (Comets et al. 2021) and the residual errors resampled from the corresponding residuals
#'   
#' **Important note:** for discrete data models, all residual-based bootstraps (residual, parametric and conditional) need 
#' a simulate.function slot to be included in the model object, as the algorithm will need to generate predictions with the 
#' resampled individual parameters in order to generate bootstrap samples. See \code{\link{SaemixModel}} and the examples
#' provided as notebooks for details.
#'   
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>
#' 
#' @references Thai H, Mentré F, Holford NH, Veyrat-Follet C, Comets E. 
#' A comparison of bootstrap approaches for estimating uncertainty of parameters in 
#' linear mixed-effects models. Pharmaceutical Statistics, 2013 ;12:129–40.
#' 
#' Thai H, Mentré F, Holford NH, Veyrat-Follet C, Comets E. 
#' Evaluation of bootstrap methods for estimating uncertainty of parameters in 
#' nonlinear mixed-effects models : a simulation study in population pharmacokinetics.
#' Journal of Pharmacokinetics and Pharmacodynamics, 2014; 41:15–33.
#' 
#' Comets E, Rodrigues C, Jullien V, Moreno U. Conditional non-parametric bootstrap for 
#' non-linear mixed effect models. Pharmaceutical Research, 2021; 38, 1057–66.
#' 
#' @examples 
#' 
#' # Bootstrap for the theophylline data 
#' data(theo.saemix)
#' 
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
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
#' saemix.model<-saemixModel(model=model1cpt,
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' saemix.options<-list(algorithm=c(1,0,0),seed=632545,save=FALSE,save.graphs=FALSE, 
#' displayProgress=FALSE)
#' 
#' # Not run (strict time constraints for CRAN)
#' \donttest{
#' saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' # Only 10 bootstrap samples here for speed, please increase to at least 200
#' theo.case <- saemix.bootstrap(saemix.fit, method="case", nboot=10)
#' theo.cond <- saemix.bootstrap(saemix.fit, nboot=10)
#' }
#' 
#' # Bootstrap for the toenail data
#' data(toenail.saemix)
#' saemix.data<-saemixData(name.data=toenail.saemix,name.group=c("id"), name.predictors=c("time","y"), 
#'  name.response="y", name.covariates=c("treatment"),name.X=c("time"))
#'  
#' binary.model<-function(psi,id,xidep) {
#'   tim<-xidep[,1]
#'   y<-xidep[,2]
#'   inter<-psi[id,1]
#'   slope<-psi[id,2]
#'   logit<-inter+slope*tim
#'   pevent<-exp(logit)/(1+exp(logit))
#'   pobs = (y==0)*(1-pevent)+(y==1)*pevent
#'   logpdf <- log(pobs)
#'   return(logpdf)
#' }
#' simulBinary<-function(psi,id,xidep) {
#'     tim<-xidep[,1]
#'     y<-xidep[,2]
#'     inter<-psi[id,1]
#'     slope<-psi[id,2]
#'     logit<-inter+slope*tim
#'     pevent<-1/(1+exp(-logit))
#'     ysim<-rbinom(length(tim),size=1, prob=pevent)
#'     return(ysim)
#'     }
#' 
#' saemix.model<-saemixModel(model=binary.model,description="Binary model",
#'      modeltype="likelihood", simulate.function=simulBinary,
#'      psi0=matrix(c(-5,-.1,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("inter","slope"))),
#'      transform.par=c(0,0))
#'      
#' \donttest{
#' saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, 
#'    nb.chains=10, fim=FALSE)
#' binary.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' # Only 10 bootstrap samples here for speed, please increase to at least 200
#' toenail.case <- saemix.bootstrap(binary.fit, method="case", nboot=10)
#' toenail.cond <- saemix.bootstrap(binary.fit, nboot=10)
#' }
#' 
#' @export saemix.bootstrap


saemix.bootstrap<-function(saemixObject, method="conditional", nboot=200, nsamp=100, saemix.options=NULL) {
  if(method!="case") {
    if(saemixObject@model@modeltype!="structural" & (is.null(body(saemixObject@model@simulate.function)) | length(formals(saemixObject@model@simulate.function))!=3)) {
      if(saemixObject@options$warnings) message("A simulation function needs to be provided for non Gaussian models to obtain bootstrap distributions by other methods than the Case bootstrap. This function needs to have the same structure as the model function and return simulated values based on the same model. \nPlease provide a simulation function the simulate.function slot of the model or use method='case' for Case bootstrap. \nExiting bootstrap.")
      return(NULL)
    }
  }
  if(method=="residual" | method=="conditional") {
    saemixObject<-saemix.predict(saemixObject) # estimate individual parameters and compute residuals (currently iwres are needed also for conditional but need to modify this in a further extension to cNP ECO TODO)
  }
  if(method=="conditional") {
    ndone <- dim(saemixObject@results@phi.samp)
    if(!is.null(ndone)) ndone<-ndone[3] else ndone<-0
    if(ndone<nsamp) {
      if(saemixObject@options$warnings) message("Not enough samples in the object, sampling from the conditional distribution\n")
      saemixObject<-conddist.saemix(saemixObject, nsamp=nsamp) # estimate conditional distributions and sample residuals
    }
    eta.sampc<-centerDist.NPcond(saemixObject, nsamp=nsamp) # Center eta samples from the conditional distribution, to avoid doing this repeatedly
  }
  if(is.null(saemix.options)) {
    #      saemix.options<-list(directory="current",fix.seed=FALSE,map=FALSE,ll.is=FALSE,displayProgress=FALSE,save.graphs=FALSE,print=FALSE)
    saemix.options<-saemixObject["options"]
    saemix.options$directory<-"current"
    saemix.options$fix.seed<-FALSE
    saemix.options$map<-FALSE   # Only parameter estimates are required for bootstrap
    saemix.options$fim<-FALSE
    saemix.options$displayProgress<-FALSE 
    saemix.options$save.graphs<-FALSE
    saemix.options$save<-FALSE
    saemix.options$ll.is<-FALSE
    saemix.options$print<-FALSE
  }
  verbose <- saemix.options$warnings
  if(saemixObject@model@modeltype=="structural") idx.eps<-saemixObject@model@indx.res else idx.eps<-integer(0)
  idx.iiv<-saemixObject@model@indx.omega
  idx.rho<-which(saemixObject@model@covariance.model[lower.tri(saemixObject@model@covariance.model)]==1)
  bootstrap.distribution<-failed.runs<-data.frame()
  nelements <- length(saemixObject@results@fixed.effects)+length(idx.iiv)+length(idx.rho)+length(idx.eps)
  # Starting point: estimates from the fit 
  model.boot<-saemixObject["model"]
  model.boot@psi0 <- model.boot["betaest.model"]
  model.boot@psi0[model.boot["betaest.model"]==1]<-saemixObject@results@fixed.effects
  for(iboot in 1:nboot) {
    if(method=="case")  
      data.boot <- dataGen.case(saemixObject)
    if(method=="residual")
      data.boot <- dataGen.NP(saemixObject,conditional=FALSE)
    if(method=="conditional")
      data.boot <- dataGen.NP(saemixObject, nsamp=nsamp,eta.sampc=eta.sampc, conditional=TRUE)
    if(method=="parametric")
      data.boot <- dataGen.Par(saemixObject)
    fit.boot<-try(saemix(model.boot, data.boot, saemix.options))
    if(is(fit.boot,"try-error")) {
      l1<-c(iboot,rep(NA,nelements))
      failed.runs <- rbind(failed.runs, c(iboot, fit.boot))
    } else {
      res<-fit.boot@results
      l1<-c(iboot,res@fixed.effects, diag(res@omega)[idx.iiv])
      if(length(idx.rho)>0) l1<-c(l1,res@omega[lower.tri(res@omega)][idx.rho])
      if(length(idx.eps)>0) l1<-c(l1, res@respar[idx.eps])
      if(length(res@ll.lin)>0) l1<-c(l1, res@ll.lin)
      
    }
#    l1<-c(iboot,res@fixed.effects, diag(res@omega)[idx.iiv],res@omega[lower.tri(res@omega)][idx.rho],res@respar[idx.eps], res@se.fixed, res@se.omega[idx.iiv],res@se.cov[lower.tri(res@se.cov)][idx.rho], res@se.respar[idx.eps],res@ll.lin)
    bootstrap.distribution<-rbind(bootstrap.distribution,l1) 
  }
  # Names
  nampar<-colnames(saemixObject@model@covariance.model)
  namcol<-c(saemixObject@results@name.fixed, saemixObject@results@name.random)
  if(length(idx.rho)>0) {
    for(i in 1:(length(nampar)-1)) {
        for(j in (i+1):length(nampar)) {
        if(saemixObject@model@covariance.model[i,j]==1) {
            namcol<-c(namcol,paste("cov.",nampar[i],nampar[j],sep=""))
        }
        }
      }
  }
  if(length(idx.eps)>0) namcol<-c(namcol,saemixObject@model@name.sigma[idx.eps])
    if(length(res@ll.lin)>0) namcol<-c(namcol,"LL.lin")
#  namcol<-c(saemixObject@results@name.fixed,saemixObject@results@name.random,namcol, saemixObject@results@name.sigma[saemixObject@results@indx.res])
#  colnames(bootstrap.distribution)<-c("Replicate",namcol,paste("SE",namcol,sep="."),"LL.lin")
  colnames(bootstrap.distribution)<-c("Replicate",namcol)
  if(verbose && dim(failed.runs)[1]>0) {
    cat(dim(failed.runs)[1],"failed:\n")
    print(head(failed.runs))
  }
  return(bootstrap.distribution)
}


#' Bootstrap datasets
#' 
#' These functions create bootstrapped datasets using the bootstrap methods described in \code{\link{saemix.bootstrap}}
#' 
#' @rdname boostrap.data
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' 
#' @aliases dataGen.case dataGen.NP dataGen.Par
#' @aliases sampDist.Par sampDist.NP sampDist.NPcond
#'  
#' @export dataGen.case


# Generating bootstrap datasets
dataGen.case<- function(saemixObject) {
  idx.boot<-sample(1:saemixObject@data@N, size=saemixObject@data@N, replace=T)
  smx.data<-saemixObject@data
  data.boot<-NULL
  for(id in idx.boot)
    data.boot<-rbind(data.boot, saemixObject@data@data[saemixObject@data@data$index==id,])
  smx.data@ntot.obs<-dim(data.boot)[1]
  smx.data@nind.obs <-saemixObject@data@nind.obs[idx.boot]
  smx.data@data<-data.boot
  orig.idx<-1:saemixObject@data@N
  smx.data@data$index <- smx.data@data[smx.data@name.group] <- rep(1:smx.data@N,times=smx.data@nind.obs)
  # As yorig and ocov not used in fit, not included in dataset
  return(smx.data)
}

#' @rdname boostrap.data
#' @param nsamp number of samples from the conditional distribution (for method="conditional")
#' @param eta.sampc if available, samples from the conditional distribution (otherwise, they are obtained within the function)
#' @param conditional if TRUE, sample from the conditional distribution, if FALSE, sample within the 
#' empirical Bayes estimates (EBE) as in the traditional non-parametric residual bootstrap
#' @export dataGen.NP

dataGen.NP<-function(saemixObject,nsamp=0,eta.sampc=NULL,conditional=FALSE) {
  if(conditional) {
    if(nsamp==0) nsamp<-dim(saemixObject@results@phi.samp)[3]
    if(is.null(eta.sampc)) x<-sampDist.NPcond(saemixObject,nsamp=nsamp) else x<-sampDist.NPcond(saemixObject,nsamp=nsamp,eta.sampc)
  } else x<-sampDist.NP(saemixObject)
  id<-saemixObject@data@data$index
  xdep<-saemixObject@data@data[,c(saemixObject@data["name.predictors"]),drop=FALSE]
  smx.data<-saemixObject@data
  if(saemixObject@model@modeltype=="structural") {
    fpred<-saemixObject@model@model(x$psi.boot,id,xdep)
    if(saemixObject@model@error.model=="exponential") fpred<-log(cutoff(fpred))
    gpred<-error.typ(fpred,saemixObject@results@respar)
    smx.data@data[,saemixObject@data["name.response"]] <- fpred+gpred*x$eps.boot # Bootstrapped data
  } else {
    smx.data@data[,saemixObject@data["name.response"]] <-saemixObject@model@simulate.function(x$psi.boot,id,xdep)
  }
  return(smx.data)
}

#' @rdname boostrap.data
#' 
#' @importFrom MASS mvrnorm
#' 
#' @export dataGen.Par

dataGen.Par<-function(saemixObject) { # Probably could reintegrate in function .NP as only the first line changes... TODO
  x<-sampDist.Par(saemixObject)
  id<-saemixObject@data@data$index
  xdep<-saemixObject@data@data[,c(saemixObject@data["name.predictors"]),drop=FALSE]
  smx.data<-saemixObject@data
  if(saemixObject@model@modeltype=="structural") {
    fpred<-saemixObject@model@model(x$psi.boot,id,xdep)
    if(saemixObject@model@error.model=="exponential") fpred<-log(cutoff(fpred))
    gpred<-error.typ(fpred,saemixObject@results@respar)
    smx.data@data[,saemixObject@data["name.response"]] <- fpred+gpred*x$eps.boot # Bootstrapped data
  } else {
    smx.data@data[,saemixObject@data["name.response"]] <-saemixObject@model@simulate.function(x$psi.boot,id,xdep)
  }
  return(smx.data)
}


sampDist.Par <- function(saemixObject) {
  omega.est<-saemixObject@results@omega[saemixObject@model@indx.omega,saemixObject@model@indx.omega, drop=FALSE] # Estimated var-cov matrix
  #sigma.est<-saemixObject@results@ # Estimated residual error
  phicond<-saemixObject@results@cond.mean.phi
  etacond<-saemixObject@results@cond.mean.eta # Estimated eta_i
  Cimu<-phicond-etacond
  eta.sim<-mvrnorm(n = dim(etacond)[1], mu=rep(0,dim(omega.est)[1]), Sigma=omega.est)
  etacond[,saemixObject@model@indx.omega]<-eta.sim
  phi.boot <- Cimu + etacond
  psi.boot <- transphi(phi.boot,saemixObject@model["transform.par"])
  if(saemixObject@model@modeltype=="structural") 
    eps.boot <- rnorm(saemixObject@data@ntot.obs)  else 
      eps.boot<-NULL
  return(list(psi.boot=psi.boot, eta.boot=etacond, eps.boot=eps.boot))
}

# Sampling distributions for non-parametric bootstraps
sampDist.NP <- function(saemixObject) {
  omega.est<-saemixObject@results@omega[saemixObject@model@indx.omega,saemixObject@model@indx.omega, drop=FALSE] # Estimated var-cov matrix
  #sigma.est<-saemixObject@results@ # Estimated residual error
  if(saemixObject@model@modeltype=="structural") {
      eps.est<-saemixObject@results@iwres # Estimated epsilon_ij
      epsc<-center.eps(eps.est)
      epscorr<-epsc/sd(epsc)
      eps.boot <- sample(epscorr,size=saemixObject@data@ntot.obs,replace=T)
  } else eps.boot<-NULL
  phicond<-saemixObject@results@cond.mean.phi
  etacond<-saemixObject@results@cond.mean.eta # Estimated eta_i
  Cimu<-phicond-etacond
  etacorr<-normalise.eta.svd(etacond[,saemixObject@model@indx.omega, drop=FALSE], omega.est)
  etasamp<-etacond
  etasamp[,saemixObject@model@indx.omega]<-etacorr
  idx.boot<-sample(1:saemixObject@data@N, size=saemixObject@data@N, replace=T)
  eta.boot<-etasamp[idx.boot,]
  phi.boot <- Cimu + eta.boot
  psi.boot <- transphi(phi.boot,saemixObject@model["transform.par"])
  return(list(psi.boot=psi.boot, eta.boot=eta.boot, eps.boot=eps.boot))
}


# Added eta.sampc=NULL argument for the simulations
# in the simulations we can create the data.frame only once and pass it to reduce computation time
sampDist.NPcond <- function(saemixObject,nsamp, eta.sampc=NULL, population=TRUE) {
  # eta.sampc: etas to sample from, centered 
  # population: centering and resampling method, either at the level of the population or at the level of the individual [only for etas for the moment]
  if(saemixObject@model@modeltype=="structural") {
    eps.est<-saemixObject@results@iwres # Estimated epsilon_ij
    epsc<-center.eps(eps.est)
    epscorr<-center.eps(eps.est)/sd(epsc) # Correcting empirical residuals
    eps.boot <- sample(epscorr,size=saemixObject@data@ntot.obs,replace=T)
  } else eps.boot<-NULL
  omega.est<-saemixObject@results@omega[saemixObject@model@indx.omega,saemixObject@model@indx.omega, drop=FALSE] # Estimated var-cov matrix
  #sigma.est<-saemixObject@results@ # Estimated residual error
  phicond<-saemixObject@results@cond.mean.phi
  etacond<-saemixObject@results@cond.mean.eta # Estimated eta_i
  Cimu<-phicond-etacond
  if(is.null(eta.sampc)) {
    eta.sampc<-centerDist.NPcond(saemixObject,nsamp,population=population)
    # x<-centerDist.NPcond(saemixObject,nsamp,population=population)
    # eta.sampc<-x$eta.sampc
    # epsc<-x$epsc
  }
  if(population) {
    idx.boot<-sample(1:dim(eta.sampc)[1], size=saemixObject@data@N, replace=T)
  } else { # resample also at the individual level
    idx.boot<-sample(1:nsamp, size=saemixObject@data@N, replace=T)
    id1<-1:saemixObject@data@N
    idx.boot<-idx.boot*(saemixObject@data@N-1)+id1 # sample nb idx.boot for subject id1
  }
  eta.boot<-eta.sampc[idx.boot,]
  phi.boot <- Cimu + eta.boot
  psi.boot <- transphi(as.matrix(phi.boot),saemixObject@model["transform.par"])
  return(list(psi.boot=psi.boot, eta.boot=eta.boot, eps.boot=eps.boot))
}

#########################################################################################
# Centering and normalising eta & eps

#' @rdname saemix.internal
#' @aliases center.eps center.eta normalise.eta normalise.eta.svd centerDist.NPcond

center.eps<-function(x) {
  x<-x-mean(x)
}
center.eta<-function(x) {
  for(i in 1:dim(x)[2]) x[,i]<-center.eps(x[,i])
  return(x)
}
normalise.eta<-function(eta.est, omega.est) {
  etacondc<-center.eta(eta.est)
  letacond<-t(chol(cov(etacondc)))
  letaest<-t(chol(omega.est))
  Acorr<-t(letaest %*% solve(letacond))
  etacorr<-etacondc %*% Acorr
  return(etacorr)
}

###### TODO - add checks for singular matrices
normalise.eta.svd<-function(eta.est, omega.est) {
  etacondc<-center.eta(eta.est)
  Semp<-1/(dim(etacondc)[1]-1) * t(etacondc) %*% etacondc
  Sempd<-svd(Semp)
  Sestd <- svd(omega.est)
  demp<-1/sqrt(Sempd$d)
  dest<-sqrt(Sestd$d)
  Acorr <- Sempd$v %*% diag(demp,nrow=length(demp)) %*%  diag(dest,nrow=length(dest)) %*% t(Sestd$v)
  etacorr<-etacondc %*% Acorr
  return(etacorr)
}

centerDist.NPcond<-function(saemixObject, nsamp, population=TRUE) {
  phicond<-saemixObject@results@cond.mean.phi
  etacond<-saemixObject@results@cond.mean.eta # Estimated eta_i
  Cimu<-phicond-etacond  
  eta.samp<-data.frame() # Conditional samples from eta distribution
  for(i in 1:nsamp) {
    eta.samp<-rbind(eta.samp,saemixObject@results@phi.samp[,,i]-Cimu)
  }
  # Centering residuals at the level of the population or at the level of the individual
  if(population) {
    eta.sampc<-center.eta(eta.samp)
  } else {
    id1<-rep(1:dim(phicond)[1],nsamp)
    eta.sampc<-eta.samp
    for(isuj in 1:dim(phicond)[1]) {
      idx<-which(id1==isuj)
      eta.sampc[idx,]<-center.eta(eta.samp[idx,])
    } 
  }
  return(eta.sampc)
}
