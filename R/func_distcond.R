#######################	Conditional means - estimates of PSI_i ########################



#' Estimate conditional mean and variance of individual parameters using the
#' MCMC algorithm
#' 
#' When the parameters of the model have been estimated, we can estimate the
#' individual parameters (psi_i).
#' 
#' Let hattheta be the estimated value of theta computed with the SAEM
#' algorithm and let p(phi_i |y_i; hattheta) be the conditional distribution of
#' phi_i for 1<=i<=N . We use the MCMC procedure used in the SAEM algorithm to
#' estimate these conditional distributions. We empirically estimate the
#' conditional mean E(phi_i |y_i; hattheta) and the conditional standard
#' deviation sd(phi_i |y_i; hattheta).
#' 
#' See PDF documentation for details of the computation. Briefly, the MCMC
#' algorithm is used to obtain samples from the individual conditional
#' distributions of the parameters. The algorithm is initialised for each
#' subject to the conditional estimate of the individual parameters obtained at
#' the end of the SAEMIX fit. A convergence criterion is used to ensure
#' convergence of the mean and variance of the conditional distributions. When
#' nsamp>1, several chains of the MCMC algorithm are run in parallel to obtain
#' samples from the conditional distributions, and the convergence criterion
#' must be achieved for all chains. When nsamp>1, the estimate of the
#' conditional mean is obtained by averaging over the different samples.
#' 
#' The shrinkage for any given parameter for the conditional estimate is
#' obtained as
#' 
#' Sh=1-var(eta_i)/omega(eta)
#' 
#' where var(eta_i) is the empirical variance of the estimates of the
#' individual random effects, and omega(eta) is the estimated variance.
#' 
#' The function adds or modifies the following elements in the results:
#' \describe{ \item{cond.mean.phi}{Conditional mean of the individual
#' distribution of the parameters (obtained as the mean of the samples)}
#' \item{cond.var.phi}{Conditional variance of the individual distribution of
#' the parameters (obtained as the mean of the estimated variance of the
#' samples)} \item{cond.shrinkage}{Estimate of the shrinkage for the
#' conditional estimates} \item{cond.mean.eta}{Conditional mean of the
#' individual distribution of the parameters (obtained as the mean of the
#' samples)} \item{phi.samp}{An array with 3 dimensions, giving nsamp samples
#' from the conditional distributions of the individual parameters}
#' \item{phi.samp.var}{The estimated individual variances for the sampled
#' parameters phi.samp} } A warning is output if the maximum number of
#' iterations is reached without convergence (the maximum number of iterations
#' is saemix.options$nbiter.saemix[2]).
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' @param nsamp Number of samples to be drawn in the conditional distribution
#' for each subject. Defaults to 1
#' @param max.iter Maximum number of iterations for the computation of the
#' conditional estimates. Defaults to twice the total number of iterations
#' (sum(saemixObject["options"]$nbiter.saemix)*2)
#' @param \dots optional arguments passed to the plots. Plots will appear if
#' the option displayProgress in the saemixObject object is TRUE
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{SaemixObject}},\code{\link{saemixControl}},\code{\link{saemix}}
#' @references Comets  E, Lavenu A, Lavielle M. Parameter estimation in nonlinear mixed effect models using saemix, an R implementation of the SAEM algorithm. Journal of Statistical Software 80, 3 (2017), 1-41.
#' 
#' Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.
#' 
#' Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm.
#' 20th meeting of the Population Approach Group in Europe, Athens, Greece
#' (2011), Abstr 2173.
#' @keywords model
#' @examples
#'  
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
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, 
#'   c("ka","V","CL"))),transform.par=c(1,1,1), 
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), 
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="constant")
#' 
#' saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
#' 
#' # Not run (strict time constraints for CRAN)
#' # saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' # saemix.fit<-conddist.saemix(saemix.fit,nsamp=3)
#' 
#' # First sample from the conditional distribution 
#' # (a N (nb of subject) by nb.etas (nb of parameters) matrix)
#' # saemix.fit["results"]["phi.samp"][,,1]
#' 
#' # Second sample
#' # saemix.fit["results"]["phi.samp"][,,2]
#' 
#' 
#' @export conddist.saemix
conddist.saemix<-function(saemixObject,nsamp=1,max.iter=NULL,...) {
  # Estimate conditional means and estimates for the individual parameters PSI_i using the MCMC algorithm
  # nsamp= number of MCMC samples
  # kmax= max nb of iterations
  # returns an array 
  saemix.data<-saemixObject["data"]
  saemix.model<-saemixObject["model"]
  N<-saemix.data["N"]
  nb.parameters<-saemixObject["model"]["nb.parameters"]
  if(is.null(max.iter)) kmax<-sum(saemixObject["options"]$nbiter.saemix)*2 else kmax<-max.iter
  # using several Markov chains
  chdat<-new(Class="SaemixRepData",data=saemixObject["data"], nb.chains=nsamp)
  NM<-chdat["NM"]
  IdM<-chdat["dataM"]$IdM
  yM<-chdat["dataM"]$yM
  XM<-chdat["dataM"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
  io<-matrix(data=0,nrow=N,ncol=max(saemixObject["data"]["nind.obs"]))
  for(i in 1:N)
    io[i,1:saemixObject["data"]["nind.obs"][i]]<-1
  ioM<-do.call(rbind,rep(list(io),nsamp))
  ind.ioM <- which(t(ioM)!=0)
  DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])
  
  ind.eta<-saemixObject["model"]["indx.omega"]
  nb.etas<-length(ind.eta)
  omega.eta<-saemixObject["results"]["omega"][ind.eta,ind.eta]
  diag.omega<-mydiag(saemixObject["results"]["omega"])
  domega2<-do.call(cbind,rep(list((sqrt(mydiag(omega.eta)))* saemixObject["options"]$rw.ini), nb.etas))
  VK<-rep(c(1:nb.etas),2)
  phi<-array(data=0,dim=c(N,nb.parameters, nsamp))
  
  domega<-cutoff(mydiag(saemixObject["results"]["omega"][ind.eta,ind.eta]), .Machine$double.eps)
  omega.eta<-saemixObject["results"]["omega"][ind.eta,ind.eta]
  omega.eta<-omega.eta-mydiag(mydiag(saemixObject["results"]["omega"][ind.eta, ind.eta]))+mydiag(domega)
  chol.omega<-chol(omega.eta)
  pres<-saemixObject["results"]["respar"]
  
  # Preparing plots
  if(saemixObject["options"]$displayProgress) {
    plot.opt<-saemixObject["prefs"]
    plot.opt$xlab<-"Iteration"
    plot.opt<-replace.plot.options(plot.opt,...)
    change.ylab<-FALSE
    if(plot.opt$ylab!=saemixObject["prefs"]$ylab & length(plot.opt$which.par)==1) change.ylab<-TRUE
    change.main<-FALSE
    if(plot.opt$main!=saemixObject["prefs"]$main & length(plot.opt$which.par)==1) change.main<-TRUE
    np<-nb.etas
    if(length(plot.opt$mfrow)==0) {
      n1<-round(sqrt(np))
      n2<-ceiling(np/n1)
      if(n1>5 | n2>5) {
        n1<-3
        n2<-4
        #      cat("Changing the plot layout\n")
      }
      plot.opt$mfrow<-c(n1,n2)
    }
  }
  # Simulation MCMC
  # initialisation a phiM=estimation des parametres individuels
  mean.phiM<-do.call(rbind,rep(list(saemixObject["results"]["mean.phi"]),nsamp))
  phiM<-do.call(rbind,rep(list(saemixObject["results"]["cond.mean.phi"]),nsamp))
  etaM<-phiM[,ind.eta]-mean.phiM[,ind.eta]  
  psiM<-transphi(phiM,saemixObject["model"]["transform.par"])
  fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
  ind.exp<-which(saemix.model["error.model"]=="exponential")
  for(ityp in ind.exp) fpred[XM$ytype==ityp]<-log(cutoff(fpred[XM$ytype==ityp]))
  if (saemixObject["model"]["type"]=="structural"){  
    gpred<-error(fpred,pres,XM$ytype)
    DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)^2+log(gpred)
  } else {
    DYF[ind.ioM]<- -fpred
  }
  U.y<-colSums(DYF)
  phiMc<-phiM
  
  econd<-sdcond<-array(0,dim=c(nb.parameters,N*kmax,nsamp)) # parametres individuels
  ebar<-sdbar<-array(0,dim=c(nb.parameters,kmax, nsamp)) # moyenne des parametres individuels
  cat("Estimating the conditional mean and variance of the distribution of individual parameters\n")
  k<-1
  while(k<=kmax) { # Set a maximum nb of iterations
    
    for(u in 1:saemixObject["options"]$nbiter.mcmc[1]) { # 1er noyau
      etaMc<-matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
      phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
      psiMc<-transphi(phiMc,saemixObject["model"]["transform.par"])
      fpred<-saemixObject["model"]["model"](psiMc, IdM, XM)
      for(ityp in ind.exp) fpred[XM$ytype==ityp]<-log(cutoff(fpred[XM$ytype==ityp]))
      if (saemixObject["model"]["type"]=="structural"){  
        gpred<-error(fpred,pres,XM$ytype)
        DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)^2+log(gpred)
      } else {
        DYF[ind.ioM]<- -fpred
      }
      Uc.y<-colSums(DYF)
      deltau<-Uc.y-U.y
      ind<-which(deltau<(-1)*log(runif(NM)))
      etaM[ind,]<-etaMc[ind,]
      U.y[ind]<-Uc.y[ind]
    }
    
    U.eta<-0.5*rowSums(etaM*(etaM%*%solve(omega.eta)))
    
    # Second stage
    
    if(saemixObject["options"]$nbiter.mcmc[2]>0) {
      nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
      nrs2<-1
      for (u in 1:saemixObject["options"]$nbiter.mcmc[2]) {
        for(vk2 in 1:nb.etas) {
          etaMc<-etaM
          etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(NM*nrs2), ncol=nrs2)%*%mydiag(domega2[vk2,nrs2],nrow=1) # 2e noyau ? ou 1er noyau+permutation?
          phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
          psiMc<-transphi(phiMc,saemixObject["model"]["transform.par"])
          fpred<-saemixObject["model"]["model"](psiMc, IdM, XM)
          for(ityp in ind.exp) fpred[XM$ytype==ityp]<-log(cutoff(fpred[XM$ytype==ityp]))
          if (saemixObject["model"]["type"]=="structural"){  
            gpred<-error(fpred,pres,XM$ytype)
            DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)^2+log(gpred)
          } else {
            DYF[ind.ioM]<- -fpred
          }
          Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
          Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%solve(omega.eta)))
          deltu<-Uc.y-U.y+Uc.eta-U.eta
          ind<-which(deltu<(-1)*log(runif(NM)))
          etaM[ind,]<-etaMc[ind,]
          U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
          U.eta[ind]<-Uc.eta[ind]
          nbc2[vk2]<-nbc2[vk2]+length(ind)
          nt2[vk2]<-nt2[vk2]+NM
        }
      }
      domega2[,nrs2]<-domega2[,nrs2]*(1+saemixObject["options"]$stepsize.rw* (nbc2/nt2-saemixObject["options"]$proba.mcmc))
    }
    
    if(saemixObject["options"]$nbiter.mcmc[3]>0) {
      nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
      nrs2<-k%%(nb.etas-1)+2
      #    if(is.nan(nrs2)) nrs2<-1 # to deal with case nb.etas=1
      for (u in 1:saemixObject["options"]$nbiter.mcmc[3]) {
        if(nrs2<nb.etas) {
          vk<-c(0,sample(c(1:(nb.etas-1)),nrs2-1))
          nb.iter2<-nb.etas
        } else {
          vk<-0:(nb.etas-1)
          nb.iter2<-1
        }
        for(k2 in 1:nb.iter2) {
          vk2<-VK[k2+vk]
          etaMc<-etaM
          etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(NM*nrs2), ncol=nrs2)%*%mydiag(domega2[vk2,nrs2])
          phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
          psiMc<-transphi(phiMc,saemixObject["model"]["transform.par"])
          fpred<-saemixObject["model"]["model"](psiMc, IdM, XM)
          for(ityp in ind.exp) fpred[XM$ytype==ityp]<-log(cutoff(fpred[XM$ytype==ityp]))
          if (saemixObject["model"]["type"]=="structural"){  
            gpred<-error(fpred,pres,XM$ytype)
            DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)^2+log(gpred)
          } else {
            DYF[ind.ioM]<- -fpred
          }
          Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
          Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%solve(omega.eta)))
          deltu<-Uc.y-U.y+Uc.eta-U.eta
          ind<-which(deltu<(-log(runif(NM))))
          etaM[ind,]<-etaMc[ind,]
          U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
          U.eta[ind]<-Uc.eta[ind]
          nbc2[vk2]<-nbc2[vk2]+length(ind)
          nt2[vk2]<-nt2[vk2]+NM
        }
      }
      domega2[,nrs2]<-domega2[,nrs2]*(1+saemixObject["options"]$stepsize.rw* (nbc2/nt2-saemixObject["options"]$proba.mcmc))
    }
    
    phiM[,ind.eta]<-mean.phiM[,ind.eta]+etaM
    if(k==1) {
      eik<-array(t(phiM),dim=c(nb.parameters,N,nsamp))
      varik<-0*eik
      ebar[,k,]<-apply(eik,c(1,3),mean)
      ij<-(k-1)*N
      econd[,(ij+1):(ij+N),]<-eik
      sdcond[,(ij+1):(ij+N),]<-varik
    } else {
      eik1<-eik
      eik<-eik*(k-1)/k+array(t(phiM),dim=c(nb.parameters,N,nsamp))/k
      varik<-varik*(k-1)/k+(array(t(phiM),dim=c(nb.parameters,N,nsamp))**2)/k+ (eik1**2)*(k-1)/k-(eik**2)
      sdik<-sqrt(varik)
      ij<-(k-1)*N
      econd[,(ij+1):(ij+N),]<-eik
      sdcond[,(ij+1):(ij+N),]<-sdik
      ebar[,k,]<-apply(eik,c(1,3),mean)
      sdbar[,k,]<-apply(sdik,c(1,3),mean)
      if(k>=saemixObject["options"]$ipar.lmcmc) {
        ibeg<-max(1,k-saemixObject["options"]$ipar.lmcmc)
        ekmax<-(ebar[,k,]+abs(ebar[,k,])* saemixObject["options"]$ipar.rmcmc)
        ekmin<-(ebar[,k,]-abs(ebar[,k,])*saemixObject["options"]$ipar.rmcmc)
        if(nsamp==1) {
          vec<-apply(ebar[,ibeg:k,],1,max)
          vec2<-apply(ebar[,ibeg:k,],1,min)
          iek<-sum(vec>ekmax)+sum(vec2<ekmin)
          vec<-apply(sdbar[,ibeg:k,],1,max)
          vec2<-apply(sdbar[,ibeg:k,],1,min)
          sdek<-sum(vec>(sdbar[,k,]+abs(sdbar[,k,])* saemixObject["options"]$ipar.rmcmc))+ sum(vec2<(sdbar[,k,]-abs(sdbar[,k,])* saemixObject["options"]$ipar.rmcmc))
        } else {
          vec<-apply(ebar[,ibeg:k,],c(1,3),max)
          vec2<-apply(ebar[,ibeg:k,],c(1,3),min)
          iek<-sum(vec>ekmax)+sum(vec2<ekmin)
          if(ibeg==1) ibeg<-2
          vec<-apply(sdbar[,ibeg:k,],c(1,3),max)
          vec2<-apply(sdbar[,ibeg:k,],c(1,3),min)
          sdek<-sum(vec>(sdbar[,k,]+abs(sdbar[,k,])* saemixObject["options"]$ipar.rmcmc))+ sum(vec2<(sdbar[,k,]-abs(sdbar[,k,])* saemixObject["options"]$ipar.rmcmc))
        }
        if(is.na(sdek)) sdek<-0
        if(iek==0 & sdek==0) break
      }
    }
    if((k%%50)==50) {
      cat(".")
      if(saemixObject["options"]$displayProgress) {
        #        par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
        par(mfrow=plot.opt$mfrow)
        limy<-cbind(apply(cbind(apply(ebar[,1:k,],1,min),ekmin),1,min), apply(cbind(apply(ebar[,1:k,],1,max),ekmax),1,max))
        for(ipar in ind.eta) {
          laby<-saemixObject["model"]["name.modpar"][ipar]
          plot(1:k,ebar[ipar,1:k,1],type="n",xlab=plot.opt$xlab, ylab=laby,ylim=limy[ipar,],main=plot.opt$main)
          for(isamp in 1:nsamp) 
            lines(1:k,ebar[ipar,1:k,isamp],col=plot.opt$col,lty=plot.opt$lty, lwd=plot.opt$lwd)
          abline(h=apply(ebar[,k,],1,mean), col=plot.opt$ablinecol, lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd)
          abline(h=apply(ekmax,1,mean)[ipar], col=plot.opt$ablinecol, lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd)
          abline(h=apply(ekmin,1,mean)[ipar], col=plot.opt$ablinecol, lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd)
        }
      }
    }
    k<-k+1
  }
  cat("\n")
  if(k>=kmax)
    cat("Computing the empirical conditional mean and variance: maximum number of iterations reached without meeting convergence criterion (max.iter=",kmax,")\n") else cat("Convergence achieved in",k,"iterations\n")
  eta.cond<-matrix(0,nrow=dim(etaM)[1],ncol=nb.parameters)
  eta.cond[,ind.eta]<-etaM
  eta.cond<-array(t(eta.cond),dim=c(nb.parameters,N,nsamp))
  eta.cond<-t(apply(eta.cond,c(1,2),mean))
  cond.shrinkage<-100*(1-apply(eta.cond,2, var)/mydiag(saemixObject["results"]["omega"]))
  names(cond.shrinkage)<-paste("Sh.",names(cond.shrinkage),".%",sep="")
  resh.eik<-resh.varik<-array(0,dim=c(nrow=N,ncol=nb.parameters,nsamp))
  psi.samp<-resh.eik
  for(isamp in 1:nsamp) {
    resh.eik[,,isamp]<-t(eik[,,isamp])
    resh.varik[,,isamp]<-t(varik[,,isamp])
    psi.samp[,,isamp]<-transphi(resh.eik[,,isamp],saemixObject["model"]["transform.par"])
  }
  
  saemixObject["results"]["cond.mean.phi"]<-apply(resh.eik,c(1,2),mean)
  saemixObject["results"]["cond.mean.psi"]<-transphi(saemixObject["results"]["cond.mean.phi"],saemixObject["model"]["transform.par"])
  saemixObject["results"]["cond.var.phi"]<-apply(resh.varik,c(1,2),mean)
  saemixObject["results"]["cond.mean.eta"]<-eta.cond
  saemixObject["results"]["cond.shrinkage"]<-cond.shrinkage
  saemixObject["results"]["phi.samp"]<-resh.eik
  saemixObject["results"]["phi.samp.var"]<-resh.varik
  saemixObject["results"]["psi.samp"]<-psi.samp
  
  # ECO TODO: verifier ce qui se passe quand un parametre n'est pas estime => faut-il prendre la moyenne des valeurs ou grace a chol.omega qui vaut 0, ca marche ? Verifier d'ailleurs si on a chol.omega quand il manque un parametre (normalement oui)
  
  return(saemixObject)
}
