###########################  Fisher Information Matrix & LL by linearisation 	#############################

#' Computes the Fisher Information Matrix by linearisation
#' 
#' Estimate by linearisation the Fisher Information Matrix and the standard
#' error of the estimated parameters.
#' 
#' The inverse of the Fisher Information Matrix provides an estimate of the
#' variance of the estimated parameters theta. This matrix cannot be computed
#' in closed-form for nonlinear mixed-effect models; instead, an approximation
#' is obtained as the Fisher Information Matrix of the Gaussian model deduced
#' from the nonlinear mixed effects model after linearisation of the function f
#' around the conditional expectation of the individual Gaussian parameters.
#' This matrix is a block matrix (no correlations between the estimated fixed
#' effects and the estimated variances).
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' @return The function returns an updated version of the object saemix.fit in
#' which the following elements have been added: \describe{
#' \item{se.fixed:}{standard error of fixed effects, obtained as part of the
#' diagonal of the inverse of the Fisher Information Matrix (only when
#' fim.saemix has been run, or when saemix.options$algorithms\[2\] is 1)}
#' \item{se.omega:}{standard error of the variance of random effects, obtained
#' as part of the diagonal of the inverse of the Fisher Information Matrix
#' (only when fim.saemix has been run, or when the saemix.options$algorithms\[2\]
#' is 1)} 
#' \item{se.res:}{standard error of the parameters of the residual error
#' model, obtained as part of the diagonal of the inverse of the Fisher
#' Information Matrix (only when fim.saemix has been run, or when the
#' saemix.options$algorithms\[2\] is 1)} 
#' \item{fim:}{Fisher Information Matrix}
#' \item{ll.lin:}{ likelihood calculated by linearisation} }
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
#' # Running the main algorithm to estimate the population parameters
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
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="constant")
#' 
#' saemix.options<-list(algorithm=c(1,0,0),seed=632545,save=FALSE,save.graphs=FALSE)
#' 
#' # Not run (strict time constraints for CRAN)
#' # saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' 
#' # Estimating the Fisher Information Matrix using the result of saemix 
#' # & returning the result in the same object
#' # fim.saemix(saemix.fit)
#' 
#' 
#' @export fim.saemix
fim.saemix<-function(saemixObject) {
  
  ind_par_struct = c(1,2,3,4,5,6)
  ntot_obs_longi = length(which(saemix.data@data$ytype==1 |saemix.data@data$ytype==2 |saemix.data@data$ytype==3))
  ind_obs_longi = which(saemix.data@data$ytype==1 |saemix.data@data$ytype==2 |saemix.data@data$ytype==3)
  
  # Estimate the Fisher Information Matrix and the s.e. of the estimated parameters  
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  xind<-saemix.data["data"][,saemix.data["name.predictors"],drop=FALSE]
  yobs<-saemix.data["data"][,saemix.data["name.response"]]
  
  #  covariance.model<-0*saemix.model["covariance.model"]
  covariance.model<-saemix.model["covariance.model"]
  omega<-saemix.res["omega"]
  omega.null<-0*omega
  #  diag(covariance.model)<-mydiag(saemix.model["covariance.model"])
  #  omega<-0*saemix.res["omega"] # Why use only diag(omega) ???
  #  diag(omega)<-mydiag(saemix.res["omega"])
  hat.phi<-saemix.res["cond.mean.phi"]
  nphi<-dim(hat.phi)[2]
  nomega<-sum(covariance.model[lower.tri(covariance.model,diag=TRUE)])
  if(length(grep("structural",saemix.model["modeltype"]))>0) {
    nres<-length(saemix.res["indx.res"])
  } else{
    nres <- 0
  }
  #nytype<-length(unique(saemix.data["data"]["ytype"])) pas bon comme Ã§a....
  nytype<-nrow(unique(saemix.data["data"]["ytype"]))
  dphi<-cutoff(abs(colMeans(hat.phi))*1e-4,1e-10)
  coefphi<-c(0,-1,1)
  
  F<-array(data=0,dim=c(saemix.data["ntot.obs"],nphi,length(coefphi)))
  #F<-array(data=0,dim=c(length(saemix.data@data[["ytype"]][saemix.data@data[["ytype"]]==1]),nphi,length(coefphi)))
  F2<-array(data=0,dim=c(length(saemix.data@data[["ytype"]][saemix.data@data[["ytype"]]%in%c(1,2,3)]),6,length(coefphi))) # 6 par longi !! a adapter !!!! 
  gs<-matrix(0,saemix.data["ntot.obs"],4)
  #gs<-matrix(0,length(saemix.data@data[["ytype"]][saemix.data@data[["ytype"]]==1]),4)
  etype.exp<-which(saemix.model["error.model"]=='exponential')
  
  ytype = saemix.data@data[["ytype"]]
  xind$ytype = ytype
  
  for (l in 1:length(coefphi)) {
    for (j in 1:nphi) {
      phi<-hat.phi
      phi[,j]<-phi[,j]+coefphi[l]*dphi[j]
      psi<-transphi(phi,saemix.model["transform.par"])
      f <- saemix.model["model"](psi, saemix.data["data"][,"index"],xind)
      for(ityp in etype.exp) f[saemix.data["data"][,saemix.data["name.ytype"]]==ityp]<-log(cutoff(f[saemix.data["data"][,saemix.data["name.ytype"]]==ityp]))    
      F[,j,l]<-f
    }
    F2[,,l] = F[which(ytype %in% c(1,2,3)),1:6,l]  ## a adapter aussi
  }
  
  
  ind.covariates<-which(saemix.model["betaest.model"]>0)
  f0<-F[,1,1]
  # g0<-cutoff(saemix.res["respar"][1]+saemix.res["respar"][2]*abs(f0))
  if(length(grep("structural",saemix.model["modeltype"]))>0) {
    g0<-error(f0,saemix.res@respar,saemix.data["data"][["ytype"]]) 
  }
  g0 = g0[ind_obs_longi]
  f0 = f0[ind_obs_longi]
  #  DF<-(F[,,3]-F[,,2])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi))/2 
  #DF<-(F[,,3]-F[,,1])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi)) #gradient of f (changed from F[,,2] to F[,,1])
  DF<-(F2[,,3]-F2[,,1])/matrix(rep(dphi[ind_par_struct],each=ntot_obs_longi), ncol=length(dphi[ind_par_struct]))
  #z<-matrix(0,saemix.data["ntot.obs"],1)
  z<-matrix(0,ntot_obs_longi,1)
  
  invVi<-Gi<-list() # Individual variance matrices
  j2<-0
  for (i in 1:saemix.data["N"]) {
    ni<-saemix.data["nind.obs"][i]-1 # -1 pour le temps de survie 
    j1<-j2+1
    j2<-j2+ni
    z[j1:j2]<-yobs[ind_obs_longi][j1:j2] - f0[j1:j2] + DF[j1:j2,,drop=FALSE]%*%hat.phi[i,ind_par_struct]
    if(length(grep("structural",saemix.model["modeltype"]))>0) {
      Vi<- DF[j1:j2,,drop=FALSE] %*% omega[ind_par_struct,ind_par_struct] %*% t(DF[j1:j2,,drop=FALSE]) + mydiag((g0[j1:j2])^2, nrow=ni)
    } else{
      Vi<- DF[j1:j2,,drop=FALSE] %*% t(DF[j1:j2,,drop=FALSE])+ mydiag(1, nrow=ni)
    }
    #    invVi[[i]]<-solve(Vi[[i]])
    # Invert avoiding numerical problems
    # Invert avoiding numerical problems
    Gi[[i]]<-round(Vi*1e12)/1e12
    VD<-try(eigen(Gi[[i]]))
    if(inherits(VD,"try-error") || det(Gi[[i]])==0) {
      cat("Unable to compute the FIM by linearisation.\n") # si matrice de variance non inversible
      stop()
    }
    D<-Re(VD$values)
    V<-Re(VD$vectors)
    invVi[[i]] <- V%*%mydiag(1/D,nrow=length(D))%*%t(V)
  }
  
  # ECO ici modifie car role de covariate.estim pas clair
  # covariate.estim=si un parametre (et ses covariables associees) sont estimees ou non
  #  covariate.estim<-matrix(rep(saemix.model["fixed.estim"], dim(saemix.model["betaest.model"])[1]),byrow=TRUE, ncol=length(saemix.model["fixed.estim"]))*saemix.model["betaest.model"] # 29/05/20
  
  covariate.estim<-matrix(saemix.model["betaest.model"][ind_par_struct],nrow=1)
  covariate.estim[1,]<-saemix.model["fixed.estim"][ind_par_struct]
  
  j<-which(saemix.model["betaest.model"][ind_par_struct]>0)
  ind.fixed.est<-(covariate.estim[j]>0)
  #npar<-sum(ind.fixed.est)
  npar = sum(ind.fixed.est)
  # Tracking indices for covariances
  myidx.omega<-c()
  myidx.cor<-c()
  name.rand1<-name.rand2<-c()
  myidx.track<-NULL
  ipar<-npar
  for(iom in 1:dim(covariance.model[ind_par_struct,ind_par_struct])[1]) {
    for(jom in iom:dim(covariance.model[ind_par_struct,ind_par_struct])[1]) {
      if(covariance.model[iom,jom]==1) {
        ipar<-ipar+1
        myidx.track<-rbind(myidx.track,c(ipar,iom,jom))
        if(iom==jom) {
          myidx.omega<-c(myidx.omega,ipar)
          name.rand1<-c(name.rand1,paste("Var",saemixObject@model@name.modpar[iom],sep="."))
          name.rand2<-c(name.rand2,paste("SD",saemixObject@model@name.modpar[iom],sep="."))
        }
        else {
          myidx.cor<-c(myidx.cor,ipar)
          name.rand1<-c(name.rand1,paste("Cov",saemixObject@model@name.modpar[iom],saemixObject@model@name.modpar[jom],sep="."))
          name.rand2<-c(name.rand2,paste("Corr",saemixObject@model@name.modpar[iom],saemixObject@model@name.modpar[jom],sep="."))
        }
      }
    }
  }
  if(length(myidx.cor)>0) {
    track.var<-myidx.track[myidx.track[,1] %in% myidx.omega,]
    for(i in myidx.cor) {
      ij<-which(myidx.track[,1]==i)
      myidx.track[ij,2]<-track.var[track.var[,2]==myidx.track[ij,2],1]
      myidx.track[ij,3]<-track.var[track.var[,2]==myidx.track[ij,3],1]
    }
  }
  if(length(grep("structural",saemix.model["modeltype"]))>0) 
    namallpar<-c(saemixObject@results@name.fixed[ind_par_struct],name.rand1, saemixObject@results@name.sigma[saemixObject@results@indx.res], name.rand2) else
      namallpar<-c(saemixObject@results@name.fixed,name.rand1, name.rand2)
  
  # hw=waitbar(1,'Estimating the population parameters (SAEM). Wait...')
  
  #ll.lin<- -0.5*saemix.data["ntot.obs"]*log(2*pi)
  ll.lin<- -0.5*ntot_obs_longi*log(2*pi)
  j2<-0
  #  indMF<-list() # Individual FIM
  MF<-matrix(0,nrow=(npar+nomega+nres),ncol=(npar+nomega+nres))
  for (i in 1:saemix.data["N"]) {
    #waitbar(i/N,hw)
    ni<-saemix.data["nind.obs"][i]-1 # -1 obs de survie 
    j1<-j2+1
    j2<-j2+ni
    yi<-yobs[ind_obs_longi][j1:j2]
    DFi<-DF[j1:j2,,drop=FALSE]
    f0i<-f0[j1:j2]
    if(length(grep("structural",saemix.model["modeltype"]))>0) {
      g0i<-g0[j1:j2]
    }
    zi<-z[j1:j2]
    Ai<-kronecker(diag(6),as.matrix(saemix.model["Mcovariates"][i,]))     #### ici nombre de param longi !!! à adapter 
    Ai<-Ai[,ind.covariates[ind_par_struct],drop=FALSE]
    DFAi<-DFi%*%Ai
    Dzi<-zi-DFAi%*%saemix.res["betas"][ind_par_struct]
    
    # Derivatives of Vi=var(yi) for subject i, w/r to lambda (FO approximation, neglecting dVi/dmu)
    DV<-list()
    for(ipar in 1:npar) {
      DV[[ipar]]<-matrix(0,ncol=ni,nrow=ni)
    }
    for(iom in 1:dim(covariance.model[ind_par_struct,ind_par_struct])[1]) {
      for(jom in iom:dim(covariance.model[ind_par_struct,ind_par_struct])[1]) {
        if(covariance.model[iom,jom]==1) {
          ipar<-ipar+1
          domega<-omega.null[ind_par_struct,ind_par_struct]
          domega[iom,jom]<-domega[jom,iom]<-1 
          #          if(iom==jom) domega[iom,jom]<-1*sqrt(omega[iom,jom]) else domega[iom,jom]<-1 # if parameterised in omega and not omega2,
          if(length(grep("structural",saemix.model["modeltype"]))>0) {
            DV[[ipar]]<-DFi %*% domega %*% t(DFi)
          } else {
            DV[[ipar]]<-DFi %*% t(DFi)
          }
        }
      }
    }
    # for(ipar in 1:nomega) {
    #   domega<-omega.null
    #   domega[ipar,ipar]<-sqrt(omega[ipar,ipar])*2
    #   DV[[ipar+npar]] <- DFi %*% t(DFi)
    # }
    
    if(length(grep("structural",saemix.model["modeltype"]))>0) {
      for(ipar.res in 1:(2*3)) {                                        #### 3 mod longi !!! à adapter 
        if(!is.na(match(ipar.res,saemix.res@indx.res))) {
          ipar<-ipar+1
          if(ipar.res%%2 == 1) DV[[ipar]]<-mydiag(2*g0i, nrow=ni) else DV[[ipar]]<-mydiag(2*g0i*f0i, nrow=ni)
        }
      }
    }
    # for(ipar.res in 1:(2*nytype)) {
    #   if(!is.na(match(ipar.res,saemix.res@indx.res))) {
    #     ipar<-ipar+1
    #     if (saemixObject["model"]["modeltype"]=="structural"){
    #       if(ipar.res%%2 == 1) DV[[ipar]]<-mydiag(2*g0i, nrow=ni) else DV[[ipar]]<-mydiag(2*g0i*f0i, nrow=ni)
    #     } else{
    #       DV[[ipar]]<-mydiag(0, nrow=ni)
    #     }
    #   }
    # }
    #    blocA <- t(DFAi) %*% invVi[[i]] %*% DFAi
    if (sum(ind.fixed.est)>0) {
      DFAiest<-DFAi[,ind.fixed.est,drop=FALSE]
      blocA<-t(DFAiest)%*% invVi[[i]] %*%DFAiest
    } else blocA<-NULL
    
    # blocAbis<-matrix(0,ncol=(npar),nrow=(npar))
    # for(ii in 1:npar) {
    #   for(ij in 1:npar) {
    #     blocAbis[ii,ij]<-DFi[,ii] %*% invVi[[i]] %*% DFi[,ij]
    #   }
    # }
      blocB<-matrix(0,ncol=(nomega+nres),nrow=(nomega+nres))
    for(ij in 1:(nomega+nres)) { # columns
      for(ii in 1:(nomega+nres)) { # lines, so that blocB is ordered according to c(covariance.model)
        blocB[ii,ij]<-sum(diag(DV[[ii+npar]] %*% invVi[[i]] %*% DV[[ij+npar]] %*% invVi[[i]] ))/2
      }
    }
    blocC<-matrix(0,ncol=(npar),nrow=(nomega+nres))
    MFi <-rbind( cbind(blocA,t(blocC)), cbind(blocC, blocB))
    MF <- MF+MFi
    #    FIMi[[i]]<-MFi
    ll.lin <- ll.lin - 0.5*log(det(Gi[[i]])) - 0.5*t(Dzi)%*% invVi[[i]] %*%Dzi 
  }
  
  for(ityp in etype.exp) ll.lin<-ll.lin-sum(yobs[saemix.data["data"][,saemix.data["name.ytype"]]==ityp])
  
  
  if (sum(ind.fixed.est)>0) {
    #Mparam<-matrix(0,dim(saemix.model["betaest.model"])[1], dim(saemix.model["betaest.model"])[2])
    Mparam<-matrix(0,dim(saemix.model["betaest.model"])[1], length(ind_par_struct))
    Mparam[1,]<-saemix.model["transform.par"][ind_par_struct]
    Mtp<-Mparam[saemix.model["betaest.model"][ind_par_struct]>0]    
    Mtp<-Mtp[ind.fixed.est][ind_par_struct]
    dbetas <- dtransphi(saemix.res["betas"][ind.fixed.est][ind_par_struct],Mtp)
    Mupth<-mydiag(1/dbetas,nrow=length(dbetas))
    Fmu<-MF[1:npar,1:npar]
    Fth<-t(Mupth)%*%Fmu%*%Mupth
    MF[1:npar,1:npar]<-Fth
    # Individual FIM
    # for(i in 1:saemix.data["N"]) {
    #   Fmui<-t(Mupth) %*% indMF[[i]][1:npar,1:npar] %*% Mupth
    #   indMF[[i]][1:npar,1:npar]<-Fmui
    # }
    Cth<-try(solve(Fth))
    if(inherits(Cth,"try-error")) {
      if(saemixObject@options$warnings) cat("Error computing the Fisher Information Matrix: singular system.\n")
      Cth<-NA*Fth
    }
  } else {
    Cth<-NULL
  }
  fim<-MF
  
  sTHest<-sqrt(mydiag(Cth))
  #sTH<-matrix(0,1,length(saemix.res["betas"]))
  sTH<-rep(0,length(saemix.res["betas"][ind_par_struct]))
  sTH[ind.fixed.est[ind_par_struct]]<-sTHest
  se.fixed<-sTH
  
  FO<-MF[-c(1:npar),-c(1:npar)]
  CO<-try(solve(FO))
  if(inherits(CO,"try-error")) {
    CO<-NA*FO
    if(saemixObject@options$warnings) cat("Error computing the Fisher Information Matrix: singular system.\n")
  }
  sO<-sqrt(mydiag(CO))
  se.omega<-matrix(0,6,1)                                                          #### 6 par !!! à adapter !!!! 
  se.sdcor<-se.cov<-matrix(0,6,6)                                        #### idem 
  se.omega[saemix.model["indx.omega"]]<-sO[myidx.omega-npar]
  se.res<-matrix(0,2*3,1)                                                       #### 3 mod longi, à adapter !!!
  if(length(grep("structural",saemix.model["modeltype"]))>0) se.res[saemix.res["indx.res"]]<-sO[(nomega+1):length(sO)]    
  # Table with SE, CV and confidence intervals
  estpar<-c(saemixObject@results@fixed.effects[ind_par_struct])
  estSE<-c(se.fixed)
  est1<-est2<-se1<-se2<-c()
  if(length(myidx.cor)>0) {   # Ã  adapter car pas fait ici...
    ipar<-npar
    for(iom in 1:6){                                                    ### 6 marq, à adapter !! 
      for(jom in iom:6) {                                               ### idem 
        if(covariance.model[iom,jom]==1) {
          ipar<-ipar+1
          se.cov[iom,jom]<-se.cov[jom,iom]<-sO[(ipar-npar)]
          est1<-c(est1,omega[iom,jom])
          se1<-c(se1,sO[ipar-npar])
          if(iom==jom) {
            se.sdcor[iom,iom]<-sqrt(CO[iom,iom])/2/sqrt(omega[iom,iom])
            est2<-c(est2,sqrt(omega[iom,iom]))
            se2<-c(se2,se.sdcor[iom,iom])
          } else { # compute correlation and SE on correlation using the delta-method
            ebet<-c(omega[iom,jom],omega[iom,iom],omega[jom,jom])
            varbet<-CO[(myidx.track[myidx.track[,1]==ipar,]-npar),(myidx.track[myidx.track[,1]==ipar,]-npar)]
            rho<-ebet[1]/sqrt(ebet[2]*ebet[3])
            debet<-c(1/sqrt(ebet[2]*ebet[3]), -ebet[1]/(ebet[2]**(3/2))/sqrt(ebet[3])/2, -ebet[1]/(ebet[3]**(3/2))/sqrt(ebet[2])/2)
            se.sdcor[iom,jom]<-se.sdcor[jom,iom]<-t(debet) %*% varbet %*% debet
            est2<-c(est2,rho)
            se2<-c(se2,se.sdcor[iom,jom])
          }
        }
      }
    }
    if(length(grep("structural",saemix.model["modeltype"]))>0) estpar<-c(estpar,est1,saemixObject@results@respar[saemixObject@results@indx.res],est2) else
      estpar<-c(estpar,est1,est2)
    if(length(grep("structural",saemix.model["modeltype"]))>0) estSE<-c(estSE,se1,se.res[saemixObject@results@indx.res],se2) else estSE<-c(estSE,se1,se2)
  } else {
    diag(se.cov)<-se.omega
    if(length(grep("structural",saemix.model["modeltype"]))>0)
      estpar<-c(estpar,diag(omega[ind_par_struct,ind_par_struct])[saemixObject@results@indx.omega],saemixObject@results@respar[saemixObject@results@indx.res], sqrt(diag(omega)[saemixObject@results@indx.omega])) else estpar<-c(estpar,diag(omega)[saemixObject@results@indx.omega], sqrt(diag(omega)[saemixObject@results@indx.omega])) 
      if(length(grep("structural",saemix.model["modeltype"]))>0)
        estSE<-c(estSE,se.omega[saemixObject@results@indx.omega],se.res[saemixObject@results@indx.res],se.omega[saemixObject@results@indx.omega]/2/sqrt(diag(omega[ind_par_struct,ind_par_struct])[saemixObject@results@indx.omega])) else estSE<-c(estSE,se.omega[saemixObject@results@indx.omega],se.omega[saemixObject@results@indx.omega]/2/sqrt(diag(omega)[saemixObject@results@indx.omega]))
  }
  conf.int<-data.frame(name=namallpar, estimate=estpar, se=estSE)
  conf.int$cv<-100*conf.int$se/conf.int$estimate
  conf.int$lower<-conf.int$estimate - 1.96*conf.int$se
  conf.int$upper<-conf.int$estimate + 1.96*conf.int$se
  
  ev1 = saemix.data@data$id[saemix.data@data$ytype==4 & saemix.data@data$obs==1]
  ev2 = saemix.data@data$id[saemix.data@data$ytype==4 & saemix.data@data$obs==2]
  ev0 = saemix.data@data$id[saemix.data@data$ytype==4 & saemix.data@data$obs==0]
  Ti = saemix.data@data$time[saemix.data@data$ytype==4]
  
  p1 = saemix.res@fixed.effects[7]
  g1 = saemix.res@fixed.effects[8]
  alpha1 = saemix.res@fixed.effects[9]
  alpha2 = saemix.res@fixed.effects[10]
  alpha3 = saemix.res@fixed.effects[11]
  
  ##################### Pour p1 
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = (-1/p1**2)+((1-exp(-g1*Ti1))/(1-p1*(1-exp(-g1*Ti1))))**2
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) g1*exp(-g1*t)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))*2*(1-exp(-g1*t))/(1-p1*(1-exp(-g1*t)))**3
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1 = function(t) g1*exp(-g1*t)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))/((1-p1*(1-exp(-g1*t)))**2)
    f_2 = function(t) g1*exp(-g1*t)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))*2*(1-exp(-g1*t))/((1-p1*(1-exp(-g1*t)))**3)
    
    fcur = f(ti)
    fprime = integrate(f_1,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprime^2*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-(fprime*exp(-fcur))**2)/((exp(-fcur)+exp(-H2)-1)**2))
  }
  
  var_p1 = -(contr1+sum(contr2))
  
  
  ##################### Pour g1 
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = (-1/g1**2)+(p1*Ti1*Ti1*exp(-g1*Ti1)*(p1-1))/(1-p1*(1-exp(-g1*Ti1)))**2
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) ((p1-1)*t*exp(g1*t)*((p1-1)*exp(g1*t)*(g1*t-2)+p1*(g1*t+2))*p1*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/((p1-1)*exp(g1*t)-p1)**3
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1 = function(t) p1*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))*((p1-1)*exp(g1*t)*(g1*t-1)+p1)/((p1-(p1-1)*exp(g1*t))**2)
    f_2 = function(t) ((p1-1)*t*exp(g1*t)*((p1-1)*exp(g1*t)*(g1*t-2)+p1*(g1*t+2))*p1*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/((p1-1)*exp(g1*t)-p1)**3
    
    fcur = f(ti)
    fprime = integrate(f_1,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprime^2*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-(fprime*exp(-fcur))**2)/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  var_g1 = -(contr1+sum(contr2))
  
  
  ##################### Pour alpha1 
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = 0
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) (g1*p1*exp(-g1*t)*(b01+b11*t-5.78)**2 *exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1 = function(t) (g1*p1*exp(-g1*t)*(b01+b11*t-5.78)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_2 = function(t) (g1*p1*exp(-g1*t)*(b01+b11*t-5.78)**2 *exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    
    fcur = f(ti)
    fprime = integrate(f_1,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprime^2*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-(fprime*exp(-fcur))**2)/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  var_alpha1 = -(contr1+sum(contr2))
  
  
  
  ##################### Pour alpha2
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = 0
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) (g1*p1*exp(-g1*t)*(b02+b12*t-7.42)**2 *exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1 = function(t) (g1*p1*exp(-g1*t)*(b02+b12*t-7.42)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_2 = function(t) (g1*p1*exp(-g1*t)*(b02+b12*t-7.42)**2 *exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    
    fcur = f(ti)
    fprime = integrate(f_1,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprime^2*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-(fprime*exp(-fcur))**2)/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  var_alpha2 = -(contr1+sum(contr2))
  
  
  ##################### Pour alpha3
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = 0
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) (g1*p1*exp(-g1*t)*(b03+b13*t-2.89)**2 *exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1 = function(t) (g1*p1*exp(-g1*t)*(b03+b13*t-2.89)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_2 = function(t) (g1*p1*exp(-g1*t)*(b03+b13*t-2.89)**2 *exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    
    fcur = f(ti)
    fprime = integrate(f_1,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprime^2*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-(fprime*exp(-fcur))**2)/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  var_alpha3 = -(contr1+sum(contr2))
  
  ##################################### dérivées croisées
  ##################### p1, g1 
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = Ti1*exp(-g1*Ti1)/((1-p1*(1-exp(-g1*Ti1)))**2)
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) ((g1*p1*t+(p1-1)*exp(g1*t)*(g1*t-1)+p1)*exp(g1*t+alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/((p1-1)*exp(g1*t)-p1)**3
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1p = function(t) g1*exp(-g1*t)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))/((1-p1*(1-exp(-g1*t)))**2)
    f_1g = function(t) p1*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))*((p1-1)*exp(g1*t)*(g1*t-1)+p1)/((p1-(p1-1)*exp(g1*t))**2)
    f_2 = function(t) ((g1*p1*t+(p1-1)*exp(g1*t)*(g1*t-1)+p1)*exp(g1*t+alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/((p1-1)*exp(g1*t)-p1)**3
    
    fcur = f(ti)
    fprimep = integrate(f_1p,0,ti)$value
    fprimeg = integrate(f_1g,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprimep*fprimeg*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-fprimep*fprimeg*(exp(-fcur)**2))/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  covar_p1g1 = -(contr1+sum(contr2))
  
  
  ##################################### dérivées croisées 
  ##################### p1, alpha1
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = 0
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) ((g1*(b01+b11*t-5.78)*exp(g1*t+alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))))/(p1-(p1-1)*exp(g1*t))**2
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1p = function(t) g1*exp(-g1*t)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))/((1-p1*(1-exp(-g1*t)))**2)
    f_1a = function(t) (g1*p1*exp(-g1*t)*(b01+b11*t-5.78)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_2 = function(t) ((g1*(b01+b11*t-5.78)*exp(g1*t+alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))))/(p1-(p1-1)*exp(g1*t))**2
    
    fcur = f(ti)
    fprimep = integrate(f_1p,0,ti)$value
    fprimea = integrate(f_1a,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprimep*fprimea*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-fprimep*fprimea*(exp(-fcur)**2))/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  covar_p1alpha1 = -(contr1+sum(contr2))
  
  ##################### p1, alpha2
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = 0
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) ((g1*(b02+b12*t-7.42)*exp(g1*t+alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))))/(p1-(p1-1)*exp(g1*t))**2
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1p = function(t) g1*exp(-g1*t)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))/((1-p1*(1-exp(-g1*t)))**2)
    f_1a = function(t) (g1*p1*exp(-g1*t)*(b02+b12*t-7.42)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_2 = function(t) ((g1*(b02+b12*t-7.42)*exp(g1*t+alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))))/(p1-(p1-1)*exp(g1*t))**2
    
    fcur = f(ti)
    fprimep = integrate(f_1p,0,ti)$value
    fprimea = integrate(f_1a,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprimep*fprimea*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-fprimep*fprimea*(exp(-fcur)**2))/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  covar_p1alpha2 = -(contr1+sum(contr2))
  
  
  ##################### p1, alpha3
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = 0
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) ((g1*(b01+b11*t-5.78)*exp(g1*t+alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))))/(p1-(p1-1)*exp(g1*t))**2
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1p = function(t) g1*exp(-g1*t)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))/((1-p1*(1-exp(-g1*t)))**2)
    f_1a = function(t) (g1*p1*exp(-g1*t)*(b03+b13*t-2.89)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_2 = function(t) ((g1*(b03+b13*t-2.89)*exp(g1*t+alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))))/(p1-(p1-1)*exp(g1*t))**2
    
    fcur = f(ti)
    fprimep = integrate(f_1p,0,ti)$value
    fprimea = integrate(f_1a,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprimep*fprimea*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-fprimep*fprimea*(exp(-fcur)**2))/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  covar_p1alpha3 = -(contr1+sum(contr2))
  
  
  ##################################### dérivées croisées 
  ##################### g1, alpha1
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = 0
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) (p1*((p1-1)*exp(g1*t)*(g1*t-1)+p1)*(b01+b11*t-5.78)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(p1-(p1-1)*exp(g1*t))**2
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1g = function(t) p1*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))*((p1-1)*exp(g1*t)*(g1*t-1)+p1)/((p1-(p1-1)*exp(g1*t))**2)
    f_1a = function(t) (g1*p1*exp(-g1*t)*(b01+b11*t-5.78)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_2 = function(t) (p1*((p1-1)*exp(g1*t)*(g1*t-1)+p1)*(b01+b11*t-5.78)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(p1-(p1-1)*exp(g1*t))**2
    
    fcur = f(ti)
    fprimeg = integrate(f_1g,0,ti)$value
    fprimea = integrate(f_1a,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprimeg*fprimea*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-fprimeg*fprimea*(exp(-fcur)**2))/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  covar_g1alpha1 = -(contr1+sum(contr2))
  
  
  ##################### g1, alpha2
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = 0
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) (p1*((p1-1)*exp(g1*t)*(g1*t-1)+p1)*(b02+b12*t-7.42)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(p1-(p1-1)*exp(g1*t))**2
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1g = function(t) p1*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))*((p1-1)*exp(g1*t)*(g1*t-1)+p1)/((p1-(p1-1)*exp(g1*t))**2)
    f_1a = function(t) (g1*p1*exp(-g1*t)*(b02+b12*t-7.42)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_2 = function(t) (p1*((p1-1)*exp(g1*t)*(g1*t-1)+p1)*(b02+b12*t-7.42)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(p1-(p1-1)*exp(g1*t))**2
    
    fcur = f(ti)
    fprimeg = integrate(f_1g,0,ti)$value
    fprimea = integrate(f_1a,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprimeg*fprimea*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-fprimeg*fprimea*(exp(-fcur)**2))/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  covar_g1alpha2 = -(contr1+sum(contr2))
  
  
  ##################### g1, alpha3
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = 0
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) (p1*((p1-1)*exp(g1*t)*(g1*t-1)+p1)*(b03+b13*t-2.89)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(p1-(p1-1)*exp(g1*t))**2
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1g = function(t) p1*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89))*((p1-1)*exp(g1*t)*(g1*t-1)+p1)/((p1-(p1-1)*exp(g1*t))**2)
    f_1a = function(t) (g1*p1*exp(-g1*t)*(b03+b13*t-2.89)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_2 = function(t) (p1*((p1-1)*exp(g1*t)*(g1*t-1)+p1)*(b03+b13*t-2.89)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(p1-(p1-1)*exp(g1*t))**2
    
    fcur = f(ti)
    fprimeg = integrate(f_1g,0,ti)$value
    fprimea = integrate(f_1a,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprimeg*fprimea*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-fprimeg*fprimea*(exp(-fcur)**2))/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  covar_g1alpha3 = -(contr1+sum(contr2))
  
  
  ##################### alpha1, alpha2
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = 0
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) (p1*g1*exp(-g1*t)*(b01+b11*t-5.78)*(b02+b12*t-7.42)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1a1 = function(t)  (g1*p1*exp(-g1*t)*(b01+b11*t-5.78)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_1a2 = function(t) (g1*p1*exp(-g1*t)*(b02+b12*t-7.42)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_2 = function(t) (p1*g1*exp(-g1*t)*(b01+b11*t-5.78)*(b02+b12*t-7.42)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    
    fcur = f(ti)
    fprimeg = integrate(f_1g,0,ti)$value
    fprimea = integrate(f_1a,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprimeg*fprimea*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-fprimeg*fprimea*(exp(-fcur)**2))/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  covar_alpha1alpha2 = -(contr1+sum(contr2))
  
  
  ##################### alpha1, alpha3
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = 0
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) (p1*g1*exp(-g1*t)*(b01+b11*t-5.78)*(b03+b13*t-2.89)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1a1 = function(t)  (g1*p1*exp(-g1*t)*(b01+b11*t-5.78)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_1a3 = function(t) (g1*p1*exp(-g1*t)*(b03+b13*t-2.89)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_2 = function(t) (p1*g1*exp(-g1*t)*(b01+b11*t-5.78)*(b03+b13*t-2.89)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    
    fcur = f(ti)
    fprimeg = integrate(f_1g,0,ti)$value
    fprimea = integrate(f_1a,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprimeg*fprimea*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-fprimeg*fprimea*(exp(-fcur)**2))/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  covar_alpha1alpha3 = -(contr1+sum(contr2))
  
  
  ##################### alpha2, alpha3
  # si ev = 1
  
  b01i = saemix.res@map.psi[,1][ev1]
  b11i = saemix.res@map.psi[,2][ev1]
  b02i = saemix.res@map.psi[,3][ev1]
  b12i = saemix.res@map.psi[,4][ev1]
  b03i = saemix.res@map.psi[,5][ev1]
  b13i = saemix.res@map.psi[,6][ev1]
  Ti1 = Ti[ev1]
  
  part1 = 0
  
  part2 = c()
  for (i in 1:length(ev1)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti1[i]
    f_2 = function(t) (p1*g1*exp(-g1*t)*(b02+b12*t-7.42)*(b03+b13*t-2.89)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    part2 = c(part2,integrate(f_2,0,ti)$value)
  }
  
  contr1 = sum(part1+part2)
  
  # si ev = 0
  
  b01i = saemix.res@map.psi[,1][ev0]
  b11i = saemix.res@map.psi[,2][ev0]
  b02i = saemix.res@map.psi[,3][ev0]
  b12i = saemix.res@map.psi[,4][ev0]
  b03i = saemix.res@map.psi[,5][ev0]
  b13i = saemix.res@map.psi[,6][ev0]
  Ti0 = Ti[ev0]
  
  contr2 = c()
  for (i in 1:length(ev0)){
    b01 = b01i[i]
    b11 = b11i[i]
    b02 = b02i[i]
    b12 = b12i[i]
    b03 = b03i[i]
    b13 = b13i[i]
    ti = Ti0[i]
    f = function(t) p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1)))
    f_1a2 = function(t)  (g1*p1*exp(-g1*t)*(b02+b12*t-7.42)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_1a3 = function(t) (g1*p1*exp(-g1*t)*(b03+b13*t-2.89)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    f_2 = function(t) (p1*g1*exp(-g1*t)*(b02+b12*t-7.42)*(b03+b13*t-2.89)*exp(alpha1*(b01+b11*t-5.78)+alpha2*(b02+b12*t-7.42)+alpha3*(b03+b13*t-2.89)))/(1-p1*(1-exp(-g1*t)))
    
    fcur = f(ti)
    fprimeg = integrate(f_1g,0,ti)$value
    fprimea = integrate(f_1a,0,ti)$value
    fsec = integrate(f_2,0,ti)$value
    
    F1 =  1-exp(-p1*g1*exp(alpha1*(b01-5.78)+alpha2*(b02-7.42)+alpha3*(b03-2.89))/(1-p1)*(exp((alpha1*b11+alpha2*b12+alpha3*b13-g1)*T)/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1*exp(-g1*T)/(1-p1))-1/(alpha1*b11+alpha2*b12+alpha3*b13-g1)*hyperg_2F1(1,1-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,2-(alpha1*b11+alpha2*b12+alpha3*b13)/g1,-p1/(1-p1))))
    H2 = -log(1-(1-F1)*(1-exp(-ti/10)))
    contr2 = c(contr2,((-fsec*exp(-fcur)+fprimeg*fprimea*exp(-fcur))*(exp(-fcur)+exp(-H2)-1)-fprimeg*fprimea*(exp(-fcur)**2))/(exp(-fcur)+exp(-H2)-1)**2)
  }
  
  covar_alpha2alpha3 = -(contr1+sum(contr2))
  
  
  
  mat = matrix(c(var_p1,covar_p1g1,covar_p1alpha1,covar_p1alpha2,covar_p1alpha3,
                 covar_p1g1,var_g1,covar_g1alpha1,covar_g1alpha2,covar_g1alpha3,
                 covar_p1alpha1,covar_g1alpha1,var_alpha1,covar_alpha1alpha2,covar_alpha1alpha3,
                 covar_p1alpha2,covar_g1alpha2,covar_alpha1alpha2,var_alpha2,covar_alpha2alpha3,
                 covar_p1alpha3,covar_g1alpha3,covar_alpha1alpha3,covar_alpha2alpha3,var_alpha3),byrow = T,nrow = 5)
  matinv = solve(mat)
  
  se_p1 = sqrt(matinv[1,1])
  se_g1 = sqrt(matinv[2,2])
  se_alpha1 = sqrt(matinv[3,3])
  se_alpha2 = sqrt(matinv[4,4])
  se_alpha3 = sqrt(matinv[5,5])
  
  
  saemix.res["se.fixed"]<-c(se.fixed,se_p1,se_g1,se_alpha1)
  saemix.res["se.omega"]<-c(se.omega)
  saemix.res["se.cov"]<-se.cov
  if(length(grep("structural",saemix.model["modeltype"]))>0) saemix.res["se.respar"]<-c(se.res)
  saemix.res["conf.int"]<-conf.int
  saemix.res["ll.lin"]<-c(ll.lin )
  saemix.res["fim"]<-fim
  saemix.res["aic.lin"]<-(-2)*saemix.res["ll.lin"]+ 2*saemix.res["npar.est"]
  saemix.res["bic.lin"]<-(-2)*saemix.res["ll.lin"]+ log(saemix.data["N"])*saemix.res["npar.est"]
  saemix.res["bic.covariate.lin"]<-(-2)*saemix.res["ll.lin"]+ log(saemix.data["N"])*saemix.res["nbeta.random"]+log(sum(saemix.data["nind.obs"]))*saemix.res["nbeta.fixed"]
  
  ##################################
  #delete(hw)
  saemixObject["results"]<-saemix.res
  return(saemixObject)
  #  return(list(ll.lin,fim,DFi, Dzi, invVi))
}




### paramÃ¨tres de survie :
#b0 = saemix.res@fixed.effects[1]
#b1 = saemix.res@fixed.effects[2]
#h0 = saemix.res@fixed.effects[3]
#alpha = saemix.res@fixed.effects[4]
#Ti = saemix.data@data$time[saemix.data@data$ytype==2]
#Ti = Ti[which(Ti!=0)]
#deltai = saemix.data@data$obs[saemix.data@data$ytype==2]

#var_h0 = 1/(sum(deltai/h0**2))
#var_alpha = 1/sum((2*h0)/(b1*alpha**3)*(exp(alpha*(b0+b1*Ti))+exp(alpha*b0))-2*(h0/(b1*alpha**2))*((b0+b1*Ti)*exp(alpha*(b0+b1*Ti))+b0*exp(alpha*b0))+
#                    (h0/(b1*alpha))*(((b0+b1*Ti)**2)*exp(alpha*(b0+b1*Ti))+b0**2*exp(alpha*b0)))

#se_alpha = sqrt(var_alpha)
#se_h0 = sqrt(var_h0)

#b0i = saemix.res@map.psi[,1]
#b1i = saemix.res@map.psi[,2]

#var_alpha = 1/sum((2*h0)/(b1i*alpha**3)*(exp(alpha*(b0i+b1i*Ti))+exp(alpha*b0i))-2*(h0/(b1i*alpha**2))*((b0i+b1i*Ti)*exp(alpha*(b0i+b1i*Ti))+b0i*exp(alpha*b0i))+
#                    (h0/(b1i*alpha))*(((b0i+b1i*Ti)**2)*exp(alpha*(b0i+b1i*Ti))+b0i**2*exp(alpha*b0i)))

#se_alpha = sqrt(var_alpha)
