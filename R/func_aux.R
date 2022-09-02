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
    if(saemixObject["model"]["modeltype"]=="structural"){
      suppressWarnings(phi1.opti<-optim(par=phi1, fn=conditional.distribution_c, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"]))
    } else {
      suppressWarnings(phi1.opti<-optim(par=phi1, fn=conditional.distribution_d, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"]))
    }
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
  cat("Computing WRES and npde ")
  isuj1<-1
  for(isuj in unique(saemixObject["data"]["data"][,saemixObject["data"]["name.group"]])) {
    if(isuj1%%10==1) cat(".")
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
  cat("\n")
  saemixObject["results"]["npde"]<-npde
  saemixObject["results"]["wres"]<-wres
  saemixObject["results"]["ypred"]<-ypredall # = [ E_i(f(theta_i)) ]
  saemixObject["results"]["pd"]<-pd
  if(length(saemixObject["results"]["predictions"])==0) saemixObject["results"]["predictions"]<-data.frame(ypred=ypredall,wres=wres,pd=pd,npde=npde) else {
  	saemixObject["results"]["predictions"]$ypred<-ypredall
  	saemixObject["results"]["predictions"]$wres<-wres
  	saemixObject["results"]["predictions"]$npde<-npde
  	saemixObject["results"]["predictions"]$pd<-pd
  }
#  return(list(ypred=ypredall,pd=pd,npd=npd,wres=wres,sim.data=ysim, sim.pred=x$sim.pred))
  return(saemixObject)
}

###########################	Computational fcts	#############################

#' @rdname saemix.internal
#' 
#' @aliases cutoff cutoff.eps cutoff.max cutoff.res
#' @aliases normcdf norminv
#' @aliases error error.typ ssq
#' @aliases transpsi transphi derivphi dtransphi
#' @aliases compute.Uy compute.LLy conditional.distribution trnd.mlx gammarnd.mlx tpdf.mlx
#' @aliases conditional.distribution_c conditional.distribution_d
#' 
#' @keywords internal

# Moved these functions to SaemixOutcome.R (define error models) in Rext
cutoff<-function(x,seuil=.Machine$double.xmin) {x[x<seuil]<-seuil; return(x)}
cutoff.max<-function(x) max(x,.Machine$double.xmin)
cutoff.eps<-function(x) max(x,.Machine$double.eps)
cutoff.res<-function(x,ares,bres) max(ares+bres*abs(x),.Machine$double.xmin)

# Inverse of the normal cumulative distribution fct: using erfcinv from ?pnorm
norminv<-function(x,mu=0,sigma=1)  mu-sigma*qnorm(x,lower.tail=FALSE)

# Truncated gaussian distribution (verifie par rapport a definition de erf/matlab)
normcdf<-function(x,mu=0,sigma=1)
  cutoff(pnorm(-(x-mu)/sigma,lower.tail=FALSE),1e-30)

error<-function(f,ab,etype) { # etype: error model
  g<-f
  for(ityp in sort(unique(etype))) {
    g[etype==ityp]<-error.typ(f[etype==ityp],ab[((ityp-1)*2+1):(ityp*2)])
  }
  return(g)
}
error.typ<-function(f,ab) {
#  g<-cutoff(ab[1]+ab[2]*abs(f))
  g<-cutoff(sqrt(ab[1]^2+ab[2]^2*f^2))  # Johannes 02/21
  return(g)
}

# ssq<-function(ab,y,f) { # Sum of squares
# 	g<-abs(ab[1]+ab[2]*f)
# 	e<-sum(((y-f)/g)**2+2*log(g))
# 	return(e)
# }
# ssq<-function(ab,y,f,ytype) { # Sum of squares
#   g<-f
#   for(ityp in sort(unique(ytype))) {
#     g[ytype==ityp]<-(ab[((ityp-1)*2+1)]+f[ytype==ityp]*ab[(ityp*2)])**2
#   }
#   e<-sum(((y-f)**2/g)+log(g))
#   return(e)
# }

ssq<-function(ab,y,f,etype) { # Sum of squares; need to put ab first as these parameters are optimised by optim
  g<-(error(f,ab,etype))
  e<-sum(((y-f)**2/g**2)+2*log(g))
  return(e)
}

transpsi<-function(psi,tr) {
  phi<-psi
#  if(is.null(dim(psi))) phi<-as.matrix(t(phi),nrow=1)
# ECO TODO: pourquoi ce test ??? Dans le cas ou psi est un vecteur ?
  i1<-which(tr==1) # log-normal
  phi[,i1]<-log(phi[,i1])
  i2<-which(tr==2) # probit
  phi[,i2]<-norminv(phi[,i2])
  i3<-which(tr==3) # logit
  phi[,i3]<-log(phi[,i3]/(1-phi[,i3]))
  if(is.null(dim(psi))) phi<-c(phi)
  return(phi)
}

transphi<-function(phi,tr) {
  psi<-phi
#  if(is.null(dim(psi))) psi<-as.matrix(t(psi),nrow=1)
  i1<-which(tr==1) # log-normal
  psi[,i1]<-exp(psi[,i1])
  i2<-which(tr==2) # probit
  psi[,i2]<-normcdf(psi[,i2])
  i3<-which(tr==3) # logit
  psi[,i3]<-1/(1+exp(-psi[,i3]))
  if(is.null(dim(phi))) psi<-c(psi)
  return(psi)
}

derivphi<-function(phi,tr) {
# Fonction calculant ??? (only used to plot parameter distributions/histograms of individual parameters)
  psi<-phi # identite
  i1<-which(tr==1) # log-normal
  psi[,i1]<-1/exp(phi[,i1])
  i2<-which(tr==2) # probit
  psi[,i2]<-1/(sqrt(2*pi))*exp(-(phi[,i2]**2)/2)
  i3<-which(tr==3) # logit
  psi[,i3]<-2+exp(phi[,i3])+exp(-phi[,i3])
  if(is.null(dim(phi))) psi<-c(psi)
  return(psi)
}

dtransphi<-function(phi,tr) {
# Fonction computing the derivative of h, used to compute the Fisher information matrix
  psi<-phi
  if(is.null(dim(phi))) {
     dpsi<-as.matrix(t(rep(1,length(phi))))
     psi<-as.matrix(t(phi),nrow=1)
  } else 
    dpsi<-matrix(1,dim(phi)[1],dim(phi)[2])
  i1<-which(tr==1) # log-normal
  dpsi[,i1]<-exp(psi[,i1])
  i2<-which(tr==2) # probit
  dpsi[,i2]<-1/dnorm(qnorm(dpsi[,i2]))   # derivee de la fonction probit, dqnorm <- function(p) 1/dnorm(qnorm(p))
  i3<-which(tr==3) # logit
  dpsi[,i3]<-1/(2+exp(-psi[,i3])+exp(psi[,i3]))
  if(is.null(dim(phi))) dpsi<-c(dpsi)
  return(dpsi)
}


compute.Uy<-function(b0,phiM,pres,args,Dargs,DYF) {
# Attention, DYF variable locale non modifiee en dehors
  args$MCOV0[args$j0.covariate]<-b0
  phi0<-args$COV0 %*% args$MCOV0
  phiM[,args$i0.omega2]<-do.call(rbind,rep(list(phi0),args$nchains))
  psiM<-transphi(phiM,Dargs$transform.par)
  if (Dargs$modeltype=="structural"){
  fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  gpred<-error(fpred,pres,Dargs$XM$ytype)
  DYF[args$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
  } else {
    fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
    for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
    DYF[args$ind.ioM]<- -fpred
  }
  U<-sum(DYF)
  return(U)
}

compute.LLy<-function(phiM,args,Dargs,DYF,pres) {
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  if (Dargs$modeltype=="structural"){
    gpred<-error(fpred,pres,Dargs$XM$ytype)
    DYF[args$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
  } else {
    DYF[args$ind.ioM]<- -fpred
  }
  U<-colSums(DYF)
  return(U)
}



conditional.distribution_c<-function(phi1,phii,idi,xi,yi,mphi,idx,iomega,trpar,model,pres,err) {
  phii[idx]<-phi1
  psii<-transphi(matrix(phii,nrow=1),trpar)
  if(is.null(dim(psii))) psii<-matrix(psii,nrow=1)
  fi<-model(psii,idi,xi)
#  if(err=="exponential") # Reverted this bit to the previous version to avoid a compiler error, not sure why it was changed...
#    fi<-log(cutoff(fi))
#  if (!(is.null(pres)) && pres[1] == pres) { # package compile throws an error when comparing a vector of length 2 (pres) to a vector of length 1
  # gi <- cutoff(pres[1])
  # } 
  # else{
  #   gi<-error(fi,pres) #    cutoff((pres[1]+pres[2]*abs(fi)))
  # }
  ind.exp<-which(err=="exponential")
  for(ityp in ind.exp) 
    fi[xi$ytype==ityp]<-log(cutoff(fi[xi$ytype==ityp]))
  gi<-error(fi,pres,xi$ytype)      #    cutoff((pres[1]+pres[2]*abs(fi)))
  Uy<-sum(0.5*((yi-fi)/gi)**2+log(gi))
  dphi<-phi1-mphi
  Uphi<-0.5*sum(dphi*(dphi%*%iomega))
  return(Uy+Uphi)
}

conditional.distribution_d<-function(phi1,phii,idi,xi,yi,mphi,idx,iomega,trpar,model) {
  phii[idx]<-phi1
  psii<-transphi(matrix(phii,nrow=1),trpar)
  if(is.null(dim(psii))) psii<-matrix(psii,nrow=1)
  fi<-model(psii,idi,xi)
  Uy <- sum(-fi)
  dphi<- phi1-mphi
  Uphi<- 0.5*sum(dphi*(dphi%*%iomega))
  return(Uy+Uphi)
}


trnd.mlx<-function(v,n,m) {
  r<-rnorm(n*m)*sqrt(v/2/gammarnd.mlx(v/2,n,m))
  return(r=matrix(r,nrow=n,ncol=m))
}

gammarnd.mlx<-function(a,n,m) {
  nm<-n*m
  y0 <- log(a)-1/sqrt(a)
  c <- a - exp(y0)
  b <- ceiling(nm*(1.7 + 0.6*(a<2)))
  y <- log(runif(b))*sign(runif(b)-0.5)/c + log(a)
  f <- a*y-exp(y) - (a*y0 - exp(y0))
  g <- c*(abs((y0-log(a))) - abs(y-log(a)))
  reject <- ((log(runif(b)) + g) > f)
  y<-y[!reject]
  if(length(y)>=nm) x<-exp(y[1:nm]) else 
    x<-c(exp(y),gammarnd.mlx(a,(nm-length(y)),1))
#  x<-matrix(x,nrow=n,ncol=m) # not useful ?
  return(x)
}

tpdf.mlx<-function(x,v) {
# TPDF_MLX  Probability density function for Student's T distribution

    term<-exp(lgamma((v + 1) / 2) - lgamma(v/2))
    return(term/(sqrt(v*pi)*(1+(x**2)/v)**((v+1)/2)))
}

###########################	Functions for npde	#############################

kurtosis<-function (x) 
{
#from Snedecor and Cochran, p 80
    x<-x[!is.na(x)]
    m4<-sum((x - mean(x))^4)
    m2<-sum((x - mean(x))^2)
    kurt<-m4*length(x)/(m2**2)-3
    return(kurtosis=kurt)
}
skewness<-function (x) 
{
#from Snedecor and Cochran, p 79
    x<-x[!is.na(x)]
    m3<-sum((x - mean(x))^3)
    m2<-sum((x - mean(x))^2)
    skew<-m3/(m2*sqrt(m2/length(x)))
    return(skewness=skew)
}
   


#' Tests for normalised prediction distribution errors
#' 
#' Performs tests for the normalised prediction distribution errors returned by
#' \code{npde}
#' 
#' Given a vector of normalised prediction distribution errors (npde), this
#' function compares the npde to the standardised normal distribution N(0,1)
#' using a Wilcoxon test of the mean, a Fisher test of the variance, and a
#' Shapiro-Wilks test for normality. A global test is also reported.
#' 
#' The helper functions \code{kurtosis} and \code{skewness} are called to
#' compute the kurtosis and skewness of the distribution of the npde.
#' 
#' @aliases testnpde kurtosis skewness
#' @param npde the vector of prediction distribution errors
#' @return a list containing 4 components: \item{Wilcoxon test of
#' mean=0}{compares the mean of the npde to 0 using a Wilcoxon test}
#' \item{variance test }{compares the variance of the npde to 1 using a Fisher
#' test} \item{SW test of normality}{compares the npde to the normal
#' distribution using a Shapiro-Wilks test} \item{global test }{an adjusted
#' p-value corresponding to the minimum of the 3 previous p-values multiplied
#' by the number of tests (3), or 1 if this p-value is larger than 1.}
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>
#' @seealso \code{\link{saemix}}, \code{\link{saemix.plot.npde}}
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F. Mentr\'e.
#' Metrics for external model evaluation with an application to the population
#' pharmacokinetics of gliclazide. \emph{Pharmaceutical Research}, 23:2036--49,
#' 2006.
#' @keywords models
#' @export testnpde
testnpde<-function(npde) 
{
    cat("---------------------------------------------\n")
    cat("Distribution of npde:\n")
    sev<-var(npde)*sqrt(2/(length(npde)-1))
    sem<-sd(npde)/sqrt(length(npde))
    cat("           mean=",format(mean(npde),digits=4),"  (SE=",format(sem,digits=2),")\n")
    cat("       variance=",format(var(npde),digits=4),"  (SE=",format(sev,digits=2),")\n")
    cat("       skewness=",format(skewness(npde),digits=4),"\n")
    cat("       kurtosis=",format(kurtosis(npde),digits=4),"\n")
    cat("---------------------------------------------\n\n")
    myres<-rep(0,4)
    y<-wilcox.test(npde)
    myres[1]<-y$p.val
    y<-shapiro.test(npde)
    myres[3]<-y$p.val

    # test de variance pour 1 ?chantillon
    # chi=s2*(n-1)/sigma0 et test de H0={s=sigma0} vs chi2 ? n-1 df
    semp<-sd(npde)
    n1<-length(npde)
    chi<-(semp**2)*(n1-1)
    y<-2*min(pchisq(chi,n1-1),1-pchisq(chi,n1-1))
    myres[2]<-y
    xcal<-3*min(myres[1:3])
    myres[4]<-min(1,xcal)
    names(myres)<-c("  Wilcoxon signed rank test ","  Fisher variance test      ",
    "  SW test of normality      ","Global adjusted p-value     ")
    cat("Statistical tests\n")
    for(i in 1:4) {
      cat(names(myres)[i],": ")
      #if (myres[i]<1) 
      cat(format(myres[i],digits=3)) 
      #else cat(myres[i])
      if(as.numeric(myres[i])<0.1 & as.numeric(myres[i])>=0.05) cat(" .")
      if(as.numeric(myres[i])<0.05) cat(" *")
      if(as.numeric(myres[i])<0.01) cat("*")
      if(as.numeric(myres[i])<0.001) cat("*")
      cat("\n")
    }
      cat("---\n")
      cat("Signif. codes: '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 \n")
    cat("---------------------------------------------\n")
    return(myres)
}
