###########################  Individual MAP estimates 	#############################

#' Estimates of the individual parameters (conditional mode)
#' 
#' Compute the estimates of the individual parameters PSI_i (conditional mode -
#' Maximum A Posteriori)
#' 
#' The MCMC procedure is used to estimate the conditional mode (or Maximum A
#' Posteriori) m(phi_i |yi ; hattheta) = Argmax_phi_i p(phi_i |yi ; hattheta)
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' @return \item{saemixObject:}{returns the object with the estimates of the
#' MAP parameters (see example for usage)}
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixObject}},\code{\link{saemix}}
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
#' @keywords models
#' @examples
#'  
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
#' saemix.options<-list(algorithm=c(1,0,0),seed=632545,
#'   save=FALSE,save.graphs=FALSE)
#' 
#' # Not run (strict time constraints for CRAN)
#' # saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' 
#' # Estimating the individual parameters using the result of saemix 
#' # & returning the result in the same object
#' # saemix.fit<-map.saemix(saemix.fit)
#' 
#' 
#' @export map.saemix

map.saemix<-function(saemixObject) {
  # Compute the MAP estimates of the individual parameters PSI_i
  i1.omega2<-saemixObject["model"]["indx.omega"]
  iomega.phi1<-solve(saemixObject["results"]["omega"][i1.omega2,i1.omega2])
  id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
  xind<-saemixObject["data"]["data"][,c(saemixObject["data"]["name.predictors"],saemixObject["data"]["name.cens"],saemixObject["data"]["name.mdv"],saemixObject["data"]["name.ytype"]),drop=FALSE]
  yobs<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  id.list<-unique(id)
  phi.map<-saemixObject["results"]["phi"]
  if(saemixObject["options"]$warnings) cat("Estimating the individual parameters, please wait a few moments...\n")
  
  if(length(grep("structural",saemixObject["model"]["modeltype"]))>0) pres<-saemixObject["results"]["respar"] else pres<-c()
  modtype<-saemixObject@model@error.model
  modtype[saemixObject@model@modeltype=="likelihood"]<-"discrete"
  
  for(i in 1:saemixObject["data"]["N"]) {
    if(saemixObject['options']$warnings) cat(".")
    isuj<-id.list[i]
    xi<-xind[id==isuj,,drop=FALSE]
    #    if(is.null(dim(xi))) xi<-matrix(xi,ncol=1)
    yi<-yobs[id==isuj]
    idi<-rep(1,length(yi))
    mean.phi1<-saemixObject["results"]["mean.phi"][i,i1.omega2]
    phii<-saemixObject["results"]["phi"][i,]
    phi1<-phii[i1.omega2]
    # conditionalLL.multi(phi1=phi1, phii=phii,idi=idi,xi=xi,yi=yi, mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=pres, modtype=modtype)
    
    phi1.opti<-optim(par=phi1, fn=conditionalLL.multi, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=pres, modtype=modtype)
    phi.map[i,i1.omega2]<-phi1.opti$par
  }
  if(saemixObject['options']$warnings)cat("\n")
  map.psi<-transphi(phi.map,saemixObject["model"]["transform.par"])
  map.psi<-data.frame(map.psi)
  map.phi<-data.frame(phi.map)
  colnames(map.psi)<-c(saemixObject["model"]["name.modpar"])
  saemixObject["results"]["map.psi"]<-map.psi
  saemixObject["results"]["map.phi"]<-map.phi
  saemixObject<-compute.eta.map(saemixObject) # compute phi, eta and shrinkage from psi
  return(saemixObject)
}

compute.eta.map<-function(saemixObject) {
  # Compute individual estimates of the MAP random effects from the MAP estimates of the parameters
  # returns the parameters (psi), newly computed if needs be, the corresponding random effects, and the associated shrinkage
  if(length(saemixObject["results"]["map.psi"])==0) {
    saemixObject<-map.saemix(saemixObject)
  }
  psi<-saemixObject["results"]["map.psi"] # [,-c(1)]
  phi<-transpsi(as.matrix(psi),saemixObject["model"]["transform.par"])
  
  # Computing COV again here (no need to include it in results)  
  #  COV<-matrix(nrow=dim(saemix.model["Mcovariates"])[1],ncol=0)
  COV<-matrix(nrow=saemixObject["data"]["N"],ncol=0)
  for(j in 1:saemixObject["model"]["nb.parameters"]) {
    jcov<-which(saemixObject["model"]["betaest.model"][,j]==1)
    aj<-as.matrix(saemixObject["model"]["Mcovariates"][,jcov])
    COV<-cbind(COV,aj)
  }
  eta<-phi-COV%*%saemixObject["results"]["MCOV"] 
  shrinkage<-100*(1-apply(eta,2,var)/mydiag(saemixObject["results"]["omega"]))
  names(shrinkage)<-paste("Sh.",saemixObject["model"]["name.modpar"],".%",sep="")
  colnames(eta)<-paste("eta.",saemixObject["model"]["name.modpar"],sep="")
  #  eta<-cbind(id=saemixObject["results"]["map.psi"][,1],eta)
  
  saemixObject["results"]["map.eta"]<-eta
  saemixObject["results"]["map.shrinkage"]<-shrinkage
  
  return(saemixObject)
}


###########################	Sum of residuals for continuous and discrete models	#############################
#' @rdname saemix.internal
#' 
#' @aliases conditionalLL.multi conditional.distribution 
#' 
#' @keywords internal

# conditionalLL.multi(phi1=phi1, phii=phii,idi=idi,xi=xi,yi=yi, mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=pres, modtype=modtype)


conditionalLL.multi<-function(phi1,phii,idi,xi,yi,mphi,idx,iomega,trpar,model,pres,modtype) {
  phii[idx]<-phi1
  psii<-transphi(matrix(phii,nrow=1),trpar)
  if(is.null(dim(psii))) psii<-matrix(psii,nrow=1)
  fi<-model(psii,idi,xi)
  ind.exp<-which(modtype=="exponential")
  for(ityp in ind.exp) 
    fi[xi$ytype==ityp]<-log(cutoff(fi[xi$ytype==ityp]))
  ind.discrete<-which(modtype=="discrete")
  ind.cont<-which(modtype!="discrete")
  # residuals from model predictions (LL or f)
  Uy<-0
  gi<-error(fi,pres,xi$ytype)      #    Note: this will probably fail when mixing TTE and longitudinal, need better way to track error models
  for(ityp in ind.cont) {
    Uy<-Uy+sum(0.5*((yi[xi$ytype==ityp]-fi[xi$ytype==ityp])/gi[xi$ytype==ityp])**2+log(gi[xi$ytype==ityp]))
  }
  for(ityp in ind.discrete)
    Uy<-Uy+sum(-fi[xi$ytype==ityp])
  
  # residuals from random effects
  dphi<-phi1-mphi
  Uphi<-0.5*sum(dphi*(dphi%*%iomega))
  return(Uy+Uphi)
}

#######################	Residuals - WRES and npde ########################

compute.sres<-function(saemixObject) {
  # Compute standardised residuals (WRES, npd and npde) using simulations
  # saemix.options$nb.sim simulated datasets used to compute npd, npde, and VPC
  # saemix.options$nb.simpred simulated datasets used to compute ypred and WRES ? for the moment saemix.options$nb.sim used for both
  nsim<-saemixObject["options"]$nb.sim
  if(length(saemixObject["sim.data"]["N"])==0 || saemixObject["sim.data"]["nsim"]!=nsim) {
    cat("Simulating data using nsim =",nsim,"simulated datasets\n")
    saemixObject<-simulate(saemixObject,nsim)
  }
  # ECO TODO: maybe here be more clever and use simulations if available (adding some if not enough, truncating if too much ?)  
  ysim<-saemixObject["sim.data"]["datasim"]$ysim
  idsim<-saemixObject["sim.data"]["datasim"]$idsim
  idy<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
  yobsall<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  ypredall<-pd<-npde<-wres<-c()
  #  pde<-c()
  #  if(saemixObject@model@modeltype=="structural" & saemixObject@options$warnings) cat("Computing WRES and npde ")
  if(sum(saemixObject@model@modeltype=="structural")>0) cat("Computing WRES and npde ") else {
    if(saemixObject@options$warnings) cat("Computing pd ")
  }
  isuj1<-1
  for(isuj in unique(saemixObject["data"]["data"][,saemixObject["data"]["name.group"]])) {
    if(isuj1%%10==1 & sum(saemixObject@model@modeltype=="structural")>0) cat(".")
    isuj1<-isuj1+1
    ysimi<-matrix(ysim[idsim==isuj],ncol=nsim)
    #    ysimi.pred<-ysimi[,1:saemixObject["options"]$nb.simpred]
    yobs<-yobsall[idy==isuj]
    tcomp<-apply(cbind(ysimi,yobs),2,"<",yobs)
    if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
    pdsuj<-rowMeans(tcomp)
    #      pdsuj[pdsuj==0]<-1/nsim
    #      pdsuj[pdsuj==1]<-1-1/nsim
    pdsuj<-pdsuj+0.5/(nsim+1)
    # pdsuj=0 pour yobs<min(tabobs$yobs)
    # pdsuj=1 : jamais
    # dc dist va de 0 ? 1-1/(nsim+1)
    # dc dist+0.5/(nsim+1) va de 0.5/(nsim+1) ? 1-0.5/(nsim+1)
    # est-ce que ?a recentre ma distribution ? ECO TODO CHECK
    pd<-c(pd,pdsuj)
    ypred<-rowMeans(ysimi)
    ypredall<-c(ypredall,ypred)
    if(sum(saemixObject@model@modeltype=="structural")>0) {
      xerr<-0
      if(length(yobs)==1) {
        npde<-c(npde,qnorm(pdsuj))
        wres<-c(wres,(yobs-ypred)/sd(t(ysimi)))
      } else {
        # Decorrelation
        vi<-cov(t(ysimi))
        xmat<-try(chol(vi))
        if(is.numeric(xmat)) {
          sqvi<-try(solve(xmat))
          if(!is.numeric(sqvi)) 
            xerr<-2
        } else xerr<-1
        if(xerr==0) {
          #decorrelation of the simulations
          decsim<-t(sqvi)%*%(ysimi-ypred)
          decobs<-t(sqvi)%*%(yobs-ypred)
          wres<-c(wres,decobs)
          #    ydsim<-c(decsim)
          #    ydobs<-decobs
          #Computing the pde
          tcomp<-apply(cbind(decsim,decobs),2,"<",decobs)
          if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
          pdesuj<-rowMeans(tcomp)
          pdesuj<-pdesuj+0.5/(nsim+1)
          #      pde<-c(pde,pdesuj)
          npde<-c(npde,qnorm(pdesuj))
        } else {
          npde<-c(npde,rep(NA,length(yobs)))
          wres<-c(wres,rep(NA,length(yobs)))
        }
      }
    }
  }
  cat("\n")
  if(sum(saemixObject@model@modeltype=="structural")>0) {
    saemixObject["results"]["npde"]<-npde
    saemixObject["results"]["wres"]<-wres
  }
  saemixObject["results"]["ypred"]<-ypredall # = [ E_i(f(theta_i)) ]
  saemixObject["results"]["pd"]<-pd
  if(length(saemixObject["results"]["predictions"])==0) {
    if(sum(saemixObject@model@modeltype=="structural")>0) saemixObject["results"]["predictions"]<-data.frame(ypred=ypredall,wres=wres,pd=pd,npde=npde) else
      saemixObject["results"]["predictions"]<-data.frame(ypred=ypredall,pd=pd)
  } else {
    saemixObject["results"]["predictions"]$ypred<-ypredall
    saemixObject["results"]["predictions"]$pd<-pd
    if(sum(saemixObject@model@modeltype=="structural")>0) {
      saemixObject["results"]["predictions"]$wres<-wres
      saemixObject["results"]["predictions"]$npde<-npde
    }
  }
  #  return(list(ypred=ypredall,pd=pd,npd=npd,wres=wres,sim.data=ysim, sim.pred=x$sim.pred))
  return(saemixObject)
}

