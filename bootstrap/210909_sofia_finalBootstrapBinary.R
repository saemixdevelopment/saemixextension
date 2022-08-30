#########################################################################################
# Bootstrap saemix

#' Bootstrap for saemix fits
#' 
#' This function 
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' @param method the name of the bootstrap algorithm to use (here : "epsilon_conditional")
#' @param nboot number of bootstrap samples
#' @param nsamp number of samples from the conditional distribution (for method="conditional")
#' @param saemix.options list of options to run the saemix algorithm. Defaults to the options in the object, but suppressing the estimation of individual parameters (map=FALSE) and likelihood (ll.is=FALSE), and the options to print out intermediate results and graphs (displayProgress=FALSE,save.graphs=FALSE,print=FALSE)
#' 
#' @details 
#' 
#' @md#' 
#' @author 
#' 
#' @references 
#'  
#' 
#' 
#' @export saemix.bootstrap


# Bootstrap function
saemix.bootstrap<-function(saemixObject, method="epsilon_conditional", nboot=200, nsamp=100, saemix.options=NULL) {
    saemixObject<-conddist.saemix(saemixObject, nsamp=nsamp) # estimated conditional distribution and sample residuals
    eta.sampc<-centerDist.NPcond(saemixObject, nsamp=nsamp) # Center eta samples from the conditional distribution, to avoid doing this repeatedly
  if(is.null(saemix.options)) {
    saemix.options<-saemixObject["options"]
    saemix.options$directory<-"current"
    saemix.options$fix.seed<-FALSE
    saemix.options$map<-FALSE   # Only parameter estimates are required for bootstrap
    saemix.options$displayProgress<-FALSE 
    saemix.options$save.graphs<-FALSE
    saemix.options$ll.is<-FALSE
    saemix.options$print<-FALSE
  }
  idx.eps<-saemixObject@model@indx.res
  idx.iiv<-saemixObject@model@indx.omega
  idx.rho<-which(saemixObject@model@covariance.model[lower.tri(saemixObject@model@covariance.model)]==1)
  bootstrap.distribution<-data.frame()
  for(iboot in 1:nboot) {
    data.boot<-dataGen.NPeps(saemixObject,nsamp=nsamp,eta.sampc=eta.sampc)
    fit.boot<-saemix(saemixObject["model"], data.boot, saemix.options)
    res<-fit.boot@results
    l1<-c(iboot,res@fixed.effects, diag(res@omega)[idx.iiv],res@omega[lower.tri(res@omega)][idx.rho],res@respar[idx.eps], res@se.fixed, res@se.omega[idx.iiv],res@se.cov[lower.tri(res@se.cov)][idx.rho], res@se.respar[idx.eps],res@ll.lin)
    bootstrap.distribution<-rbind(bootstrap.distribution,l1) 
  }
  # Names
  nampar<-colnames(saemixObject@model@covariance.model)
  namcol<-c()
  for(i in 1:(length(nampar)-1)) {
    for(j in (i+1):length(nampar)) {
      if(saemixObject@model@covariance.model[i,j]==1) {
        namcol<-c(namcol,paste("cov.",nampar[i],nampar[j],sep=""))
      }
    }
  }
  namcol<-c(saemixObject@results@name.fixed,saemixObject@results@name.random,namcol, saemixObject@results@name.sigma[saemixObject@results@indx.res])
  colnames(bootstrap.distribution)<-c("Replicate",namcol,paste("SE",namcol,sep="."),"LL.lin")
  return(bootstrap.distribution)
}

#' Bootstrap datasets
#' 
#' These functions create bootstrapped datasets using the bootstrap methods described in \code{\link{bootstrap.saemix}}
#' 
#' @rdname boostrap.data
#' @aliases dataGen.NPeps
#' @aliases sampDist.epsCond, sampDist.NPcond 
#'  
#' @export dataGen.NPeps

dataGen.NPeps<-function(saemixObject,nsamp=nsamp,eta.sampc=eta.sampc) {
  eps<-sampDist.epsCond(saemixObject,nsamp=nsamp, eta.sampc=eta.sampc)
  x<-sampDist.NPcond(saemixObject,nsamp=nsamp, eta.sampc=eta.sampc)
  id<-saemixObject@data@data$index
  xdep<-saemixObject@data@data[,c(saemixObject@data["name.predictors"]),drop=FALSE]
  fpred<-saemixObject@model@model(x$psi.boot,id,xdep)
  if(saemixObject@model@error.model=="exponential") fpred<-log(cutoff(fpred))
  gpred<-error.typ(fpred,saemixObject@results@respar)
  smx.data<-saemixObject@data
  smx.data@data[,saemixObject@data["name.response"]] <- fpred+gpred*eps$eps.boot # Bootstrapped data
  return(smx.data)
}

sampDist.epsCond <-function(saemixObject,nsamp,eta.sampc=eta.sampc) {
  eps<-NULL
  psi<-phi.tot<-NULL
  f<-g<-NULL
  omega.est<-saemixObject@results@omega[saemixObject@model@indx.omega,saemixObject@model@indx.omega] # Estimated var-cov matrix
  for (isamp in 1:nsamp) {
    #calculate fpred,gpred  
    phi.samp<-as.data.frame(saemixObject@results@phi.samp[,,isamp])
    phi.tot<-rbind(phi.tot,phi.samp)
    psi.samp<-transphi(as.matrix(phi.samp),saemixObject@model["transform.par"])
    psi<-as.data.frame(rbind(psi,psi.samp))
    id1<-saemix.data["data"][,"index"]
    xidep1<-as.data.frame(saemixObject@data@data[,saemixObject@data@name.predictors])
    fpred<-saemixObject@model@model(psi=psi.samp, id=id1, xidep=xidep1)
    f<-c(f,fpred)
    gpred<-error.typ(fpred,ab=saemixObject@results@respar)
    g<-c(g,gpred)
    
    #calculate epsilon
    response<-saemixObject@data@data[,saemixObject@data@name.response]
    epsilon<-(response-fpred)/gpred
    eps<-c(eps,epsilon)
  }
  
  eps<-center.eps(eps)
  eps.boot <- sample(eps,saemixObject@data@ntot.obs,replace = TRUE)
  
  return(list(eps.boot=eps.boot))
}

sampDist.NPcond <- function(saemixObject,nsamp,eta.sampc=NULL, population=TRUE) {
  phicond<-saemixObject@results@cond.mean.phi
  etacond<-saemixObject@results@cond.mean.eta # Estimated eta_i
  Cimu<-phicond-etacond
  if(is.null(eta.sampc)) {
    eta.sampc<-centerDist.NPcond(saemixObject,nsamp,population=population)
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
  return(list(psi.boot=psi.boot, eta.boot=eta.boot))
}

#########################################################################################
# Centering  eta & eps

#' @rdname saemix.internal
#' @aliases center.eps center.eta centerDist.NPcond

center.eps<-function(x) {
  x<-x-mean(x)
}
center.eta<-function(x) {
  for(i in 1:dim(x)[2]) x[,i]<-center.eps(x[,i])
  return(x)
}

centerDist.NPcond<-function(saemixObject, nsamp, population=FALSE) {
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
