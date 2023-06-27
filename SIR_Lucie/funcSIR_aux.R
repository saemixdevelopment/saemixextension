#########################################################################
####################### FUNCTIONS used in SaemixSIR #####################
#########################################################################


###### Index of covariance parameters in omega ######
indx.covomega <- function(SaemixObject){
  covariance.model <- SaemixObject['model']['covariance.model']
  omega <- SaemixObject['results']['omega']
  vec <- omega[lower.tri(covariance.model, diag=T)]
  omega[upper.tri(covariance.model)] <- 0
  indx.covomega <- c()
  for (i in vec){
    if (i!=0){
      ind <- which(omega==i, arr.ind=T)
      if (ind[1,1]!=ind[1,2]){
        indx.covomega <- rbind(indx.covomega, ind)
      }
    }
  }
  colnames(indx.covomega) <- c()
  row.names(indx.covomega) <- c()
  return(indx.covomega) # returns a matrix, each row contains the 2 index of covariances (in the omega matrix)
}

###### Names of covariance parameters ######
name.covparam <- function(SaemixObject){
  indx.covomega<- indx.covomega(SaemixObject)
  if (length(indx.covomega)==0){
    return()
  }
  name.fixed <- SaemixObject['model']['name.fixed']
  name <- c()
  for (i in 1:nrow(indx.covomega)){
    cov <- indx.covomega[i,]
    newname <- paste('cov(',name.fixed[cov[2]],',', name.fixed[cov[1]],')', sep='')
    name <- c(name, newname)
  }
  return(name)
}

###### Names of parameters in omega (random and covariance) ######

name.covomega <- function(SaemixObject){
  name.covparam <- name.covparam(SaemixObject)
  name.omega <- SaemixObject['model']['name.random']
  covariance.model <- SaemixObject['model']['covariance.model']
  omega <- SaemixObject['results']['omega']
  vec <- omega[lower.tri(covariance.model, diag=T)]
  omega[upper.tri(covariance.model)] <- 0
  name <- c()
  o <- 1
  c <- 1
  for (i in vec){
    if (i!=0){
      ind <- which(omega==i, arr.ind=T)
      if (ind[1,1]!=ind[1,2]){
        name <- c(name, name.covparam[c])
        c <- c+1
      }
      else {
        name <- c(name, name.omega[o])
        o <- o+1
      }
    }
  }
  return(name)
}




###### Estimated Parameters of SaemixObject into vector ######
estpar.vector <- function(SaemixObject){
  indx.fix <- SaemixObject['results']['indx.fix']
  indx.cov <- SaemixObject['results']['indx.cov']
  indx.omega <- SaemixObject['results']['indx.omega']
  indx.res <- SaemixObject['results']['indx.res']

  indx.fixed <- sort(c(indx.fix, indx.cov))
  thetamu<- SaemixObject["results"]["fixed.effects"][indx.fixed]
  
  omega <- SaemixObject['results']['omega']
  ltriomega <- omega[lower.tri(omega, diag=T)]
  ind <- which(ltriomega!=0)
  omegamu <- ltriomega[ind] 
  resparmu <- SaemixObject['results']['respar'][indx.res]
  mu  <- c(thetamu,omegamu,resparmu)
  mu <- mu[!is.na(mu)]
  return(mu)
  
}

###### Estimated Parameters of SaemixObject from vector to different matrix ######
estpar.fromvect.tomatrix <- function(newparam, SaemixObject){
  npar.est <- length(estpar.vector(SaemixObject))
  indx.fix <- SaemixObject['results']['indx.fix']
  indx.cov <- SaemixObject['results']['indx.cov']
  indx.omega <- SaemixObject['results']['indx.omega']
  indx.res <- SaemixObject['results']['indx.res']
  param <- newparam
  indx.fixed <- sort(c(indx.fix, indx.cov))
  
  newfixed.effects <- newparam[indx.fixed]
  param <- param[-c(indx.fixed)]
  
  covariance.model <- SaemixObject['model']['covariance.model']
  vec <- covariance.model[lower.tri(covariance.model, diag=T)]
  for (i in 1:length(vec)){
    if (vec[i]!=0){
      vec[i]<- param[1]
      param <- param[-1]
    }
  }
  p <- ncol(covariance.model)
  newomega <- diag(p) 
  newomega[lower.tri(newomega, diag=TRUE)] <- vec
  newomega <- newomega + t(newomega) - diag(diag(newomega))
  
  newres <- c(0,0)
  if(length(param)!=0) newres[indx.res] <- param
  
  return(list(newfixed.effects=newfixed.effects, newomega=newomega, newres=newres))
}



################ ll.lin ##################

lllin.saemix<-function(saemixObject) {
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  xind<-saemix.data["data"][,saemix.data["name.predictors"],drop=FALSE]
  yobs<-saemix.data["data"][,saemix.data["name.response"]]
  
  #  covariance.model<-0*saemix.model["covariance.model"]
  covariance.model<-saemix.model["covariance.model"]
  omega<-saemix.res["omega"]

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
  
  
  ind.covariates<-which(saemix.model["betaest.model"]>0)
  f0<-F[,1,1]
  # g0<-cutoff(saemix.res["respar"][1]+saemix.res["respar"][2]*abs(f0))
  if (saemixObject["model"]["modeltype"]=="structural"){
    g0<-error(f0,saemix.res@respar,saemix.data["data"]["ytype"]) 
  }
  
  #  DF<-(F[,,3]-F[,,2])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi))/2 
  DF<-(F[,,3]-F[,,1])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi)) #gradient of f (changed from F[,,2] to F[,,1])
  z<-matrix(0,saemix.data["ntot.obs"],1)
  
  invVi<-Gi<-list() # Individual variance matrices
  j2<-0
  for (i in 1:saemix.data["N"]) {
    ni<-saemix.data["nind.obs"][i]
    j1<-j2+1
    j2<-j2+ni
    z[j1:j2]<-yobs[j1:j2] - f0[j1:j2] + DF[j1:j2,,drop=FALSE]%*%hat.phi[i,]
    if (saemixObject["model"]["modeltype"]=="structural"){
      Vi<- DF[j1:j2,,drop=FALSE] %*% omega %*% t(DF[j1:j2,,drop=FALSE]) + mydiag((g0[j1:j2])^2, nrow=ni)
    } else{
      Vi<- DF[j1:j2,,drop=FALSE] %*% t(DF[j1:j2,,drop=FALSE])+ mydiag(1, nrow=ni)
    }
    #    invVi[[i]]<-solve(Vi[[i]])
    # Invert avoiding numerical problems
    Gi[[i]]<-round(Vi*1e10)/1e10
    VD<-try(eigen(Gi[[i]]))
    if(class(VD)=="try-error") {
      cat("Unable to compute the FIM by linearisation.\n")
      stop()
    
    }
    D<-Re(VD$values)
    V<-Re(VD$vectors)
    invVi[[i]] <- V%*%mydiag(1/D,nrow=length(D))%*%t(V)
  }
  
  
  
  ll.lin<- -0.5*saemix.data["ntot.obs"]*log(2*pi)
  j2<-0
  
  for (i in 1:saemix.data["N"]) {
    #waitbar(i/N,hw)
    ni<-saemix.data["nind.obs"][i]
    j1<-j2+1
    j2<-j2+ni
    yi<-yobs[j1:j2]
    DFi<-DF[j1:j2,,drop=FALSE]
    f0i<-f0[j1:j2]
    if (saemixObject["model"]["modeltype"]=="structural"){
      g0i<-g0[j1:j2]
    }
    zi<-z[j1:j2]
    Ai<-kronecker(diag(nphi),as.matrix(saemix.model["Mcovariates"][i,]))
    Ai<-Ai[,ind.covariates,drop=FALSE]
    DFAi<-DFi%*%Ai
    Dzi<-zi-DFAi%*%saemix.res["betas"]
    
    
    ll.lin <- ll.lin - 0.5*log(det(Gi[[i]])) - 0.5*t(Dzi)%*% invVi[[i]] %*%Dzi 
  }
  
  for(ityp in etype.exp) ll.lin<-ll.lin-sum(yobs[saemix.data["data"][,saemix.data["name.ytype"]]==ityp])
  
  
  saemix.res["ll.lin"]<-c(ll.lin )
  
  return(ll.lin)
  
}



