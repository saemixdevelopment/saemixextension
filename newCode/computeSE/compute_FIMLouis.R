# WIP... 
# Compute LLi (complete likelihood)

########## Move to aux
# Function to convert a symmetric matrix to a vector
convert.cov2vec <- function(xmat, removeZero = FALSE) {
  # Convert a symmetric non-negative matrix into a vector (no check)
  if(removeZero) { # keep only non-zero lines/columns
    x <- which(abs(mydiag(xmat))>.Machine$double.xmin)
    xmat <- xmat[x, x, drop=FALSE]
  }
  # return first the diagonal then the off-diagonal elements
  xoff<-xmat[lower.tri(xmat)]
  if(removeZero) c(mydiag(xmat), xoff[abs(xoff)>.Machine$double.xmin]) else c(mydiag(xmat), xoff)
}

# Function to create a symmetric matrix with specified variances and covariance terms
convert.vec2cov <- function(size, variances, indx.var, covariances=NULL, indx.covar=NULL) {
  # input: 
  ## size: matrix size
  ## variances: vector of elements on the diagonal
  ## indx.var: indexes of the elements
  ## covariances: vector of off-diagonal elements
  ## indx.covar: indexes of the off-diagonal elements
  xmat<-matrix(data=0, nrow=size, ncol= size)
  diag(xmat)[indx.var] <- variances
  if(!is.null(covariances)) {
    xmat[lower.tri(xmat)][indx.covar]<-covariances
    xmat[upper.tri(xmat)][indx.covar]<-covariances
  }
  return(xmat)
} 

########## 

compute.completeLL <- function(theta, phi, args, Dargs, DYF) {
  # input
  ## theta: vector of parameters (fixed+random+residual)
  ## phi: individual parameters (Nxnpar)
  ## args: a list containing
  ### ind.ioM: which elements of DYF are to be filled
  ### MCOV: parameter matrix
  ### NCOV: design matrix
  ### indices used to convert omega to a vector and back
  ### indices 
  ### structural model, model type, parameter transformation, type of error model (ytype, etype.exp)
  ### data (XM and IdM)
  ## args: 
  ## DYF: a matrix to hold intermediate computations
    # compute.LLtheta
  if(length(args$idx.covar)>0) 
    omega <- convert.vec2cov(size=args$npar, variances=theta[args$idx.omega], indx.var=args$indx.omega, covariances = theta[args$idx.covar ], indx.covar = args$indx.ltri.covar) else
      omega <- convert.vec2cov(size=args$indx.omega, variances=theta[args$idx.omega], indx.omega=args$indx.omega)
    parfix <- theta[1:args$nfix]
    parfix[args$idx.mu]<-transpsi(t(parfix[args$idx.mu]), Dargs$transform.par)
    args$MCOV[args$idx.mcov]<-parfix
    mean.phi <- args$NCOV %*% args$MCOV
    lltheta <- compute.LLtheta(mu=mean.phi[,args$indx.omega,drop=FALSE], 
                               omega=omega[args$indx.omega, args$indx.omega, drop=FALSE], 
                               phi=phi[,args$indx.omega,drop=FALSE]) 
    # compute.LLy
    pres<-c(0,0)
    pres[args$indx.res]<-theta[args$idx.res]
    lly <- compute.LLy(phi,args,Dargs,DYF,pres)
    return((lly+lltheta))
}

# ToDo: check dimensions when phi is given for all subjects
compute.LLtheta <- function(mu, omega, phi) {
  # compute LLtheta over the parameters with variability
  # mu (= mu+betaCOV) with corresponding IIV matrix omega
  # input
  ## mu: population parameters (matrix of size N*saemixObject@results@indx.omega)
  ## omega: variance-covariance matrix (square matrix of size saemixObject@results@indx.omega)
  ## phi: individual parameters (same dimension as mu)
  invOmega <- solve(omega)
  diff.phi <- (phi-mu)
  lldiff <- c()
  for(i in 1:dim(diff.phi)[1]) {
    lldiff<-c(lldiff,
              diff.phi[i,,drop=FALSE] %*% invOmega %*% t(diff.phi[i,,drop=FALSE]))
    
  }
  return((lldiff-log(det(omega)))/2)
#  -log(det(omega))/2-(phi-mu) %*% solve(omega) %*% t(phi-mu) /2 # ici le mu représente mu+betaCOV
}

# when omega diagonal
# -log(sqrt(omega))-(phiM2-(mu))**2 /(2*omega) # ici le mu représente mu+betaCOV



# LLy: compute.LLy
compute.LLy<-function(phiM,args,Dargs,DYF,pres) {
  # input
  ## phiM: individual parameters
  ## Dargs: 
  ### structural model, model type, parameter transformation, type of error model (ytype, etype.exp)
  ### data (XM and IdM)
  ## args: 
  ### ind.ioM: which elements of DYF are to be filled
  ## pres: residual error model (when y Gaussian)
  ## DYF: matrix storing the likelihood for the individual observations
  # returns
  ## U: vector with the contribution of each individual to the likelihood
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


# juste pour le calcul des SE, meme fonction mais moyennee sur les phiM (cas où il y a plusieurs chaînes)
compute.LLy.multi2<-function(phiM,args,Dargs,DYF,pres) {
  DYF = DYF[,1:Dargs$N]
  args$ind.ioM =  args$ind.ioM[1:length(Dargs$yobs)]
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM[1:length(Dargs$yobs)],Dargs$XM[1:length(Dargs$yobs),])
  ytype = Dargs[["XM"]][["ytype"]][1:length(Dargs$yobs)]
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  gpred<-error(fpred,pres,Dargs$XM$ytype[1:length(Dargs$yobs)])
  for(itype in 1:length(Dargs$modeltype)) {
    if (Dargs$modeltype[itype]=="structural"){
      DYF[args$ind.ioM][ytype==itype]<-0.5*((Dargs$yM[1:length(Dargs$yobs)][ytype==itype]-fpred[ytype==itype])/gpred[ytype==itype])**2+log(gpred[ytype==itype])
    } else {
      DYF[args$ind.ioM][ytype==itype]<- -fpred[ytype==itype]
    }
  }
  U<-colSums(DYF)
  return(U)
}


# LLphi: portion of the complete log-likelihood for the parameter
compute.LLphi<-function(phiM, theta, omega) {
  # input
  ## phiM: individual parameters (a matrix with nsuj row and npar parameters)
  ## theta: fixed effects (mu and beta)
  ## omega: variance-covariance matrix
  # returns
  ## llphi: vector with the contribution of each subject's individual parameters to LLi
  
  
  
  return(llphi)
}

gradient.LLphi<-function(phiM,args,Dargs) {
  # returns a vector with the contribution of each subject's individual parameter to d_theta(LLi)
  
  return(dLLphi)
}
hessian.LLphi<-function(phiM,args,Dargs) {
  # returns a vector with the contribution of each subject's individual parameter to d2_theta(LLi)
  
  return(d2LLphi)
}

initialiseComputationLL <- function(saemixObject) {
  # extract elements 
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
  if (saemixObject["model"]["modeltype"]=="structural"){
    nres<-length(saemix.res["indx.res"])
  } else{
    nres <- 0
  }
  nytype<-length(unique(saemix.data["data"]["ytype"]))
  dphi<-cutoff(abs(colMeans(hat.phi))*1e-4,1e-10)
  coefphi<-c(0,-1,1)
  
  F<-array(data=0,dim=c(saemix.data["ntot.obs"],nphi,length(coefphi)))
  gs<-matrix(0,saemix.data["ntot.obs"],4)
  etype.exp<-which(saemix.model["error.model"]=='exponential')
  
  for (l in 1:length(coefphi)) {
    for (j in 1:nphi) {
      phi<-hat.phi
      phi[,j]<-phi[,j]+coefphi[l]*dphi[j]
      psi<-transphi(phi,saemix.model["transform.par"])
      f <- saemix.model["model"](psi, saemix.data["data"][,"index"],xind)
      for(ityp in etype.exp) f[saemix.data["data"][,saemix.data["name.ytype"]]==ityp]<-log(cutoff(f[saemix.data["data"][,saemix.data["name.ytype"]]==ityp]))    
      F[,j,l]<-f
    }
  }
  #  DF<-(F[,,3]-F[,,2])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi))/2 
  DF<-(F[,,3]-F[,,1])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi)) #gradient of f (changed from F[,,2] to F[,,1])
  z<-matrix(0,saemix.data["ntot.obs"],1)
  
  ind.covariates<-which(saemix.model["betaest.model"]>0)
  f0<-F[,1,1]
  
  # Covariance model
  omega <- saemix.res["omega"]
  idx.omega<-which(omega[lower.tri(omega, diag=T)]>0)
  
  # Design matrices
  Mcov <- vector(mode='list', length=saemix.data["N"])
  j2<-0
  for (i in 1:saemix.data["N"]) {
    ni<-saemix.data["nind.obs"][i]
    j1<-j2+1
    j2<-j2+ni
    z[j1:j2]<-yobs[j1:j2] - f0[j1:j2] + DF[j1:j2,,drop=FALSE]%*%hat.phi[i,]
    Ai<-kronecker(diag(nphi),as.matrix(saemix.model["Mcovariates"][i,]))
    Mcov[[i]]<-Ai[,ind.covariates,drop=FALSE]
  }
  ndat.exp<-0
  for(ityp in etype.exp) ndat.exp<-ndat.exp+sum(yobs[saemix.data["data"][,saemix.data["name.ytype"]]==ityp])
  
  theta<-c(saemix.res@fixed.effects, omega[lower.tri(omega,diag=T)][idx.omega])
  if(saemix.model@modeltype=="structural") theta<-c(theta, saemix.res["respar"][saemix.res["indx.res"]])
  
  return(list(theta=theta, modeltype=saemix.model@modeltype, etype=saemix.data["data"][["ytype"]], f0=f0, DF0=DF, z=z, Mcov=Mcov, nfix=length(saemix.res@betas), nomega=dim(omega)[1], nind.obs=saemix.data["nind.obs"], ndat.exp=ndat.exp, idx.fix=saemix.res["indx.fix"], idx.omega=idx.omega, idx.res=saemix.res["indx.res"], tr.fix=saemix.model["transform.par"]))
}
