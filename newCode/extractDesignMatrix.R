# Combine saemix data and model to get the model structure and design matrix

extractModelStructure <- function(saemix.data, saemix.model, level=1, verbose=FALSE) {
  if(is.numeric(level)) {
    if(level<=length(names(saemix.model@var.model))) level<-names(saemix.model@var.model)[level] else level<-names(saemix.model@var.model)[1]
  }
  var.model<-saemix.model@var.model[[level]]
  covariate.model<-saemix.model@covariate.model
  if(length(covariate.model)==0) dim(covariate.model)[2]<-saemix.model@npar
  covariate.model.fix<-saemix.model@covariate.model.fix
  if(dim(covariate.model.fix)[2]==0) dim(covariate.model.fix)[2]<-saemix.model@npar
  betaest.model<-rbind(rep(1,saemix.model@npar), covariate.model)
  name.variable<-var.model@variable
  data <- saemix.data@data

  #  if(is.character(level)) level<-match(level, names(object@var.model))
  index.omega.var <- var.model@index.omega.var
  index.omega.novar <- var.model@index.omega.novar
  idxmat.omega <- var.model@idxmat.omega
  # Matching covariates in the model and the data
  ivar<-data[, name.variable]
  if(length(covariate.model)==0) {
    idx<-which(rownames(covariate.model) %in% saemix.data@name.covariates)
    covariate.model<-covariate.model[idx,]
    covariate.model.fix<-covariate.model.fix[idx,]
  }
  N<-length(unique(ivar))
  # Covariate model & design matrix
  ## Mcovariate: nrow=N, ncol=1 (id) + nb of covariates in the model
  ### for each subject, 1 in the ID column and the values of the covariates for that subject 
  ## TODO: here or before ? apply transformations before creating Mcovariates [see SaemixcovariateTransform.R]
  ### continuous covariates => replace with value
  ### binary covariates => transform to 0/1
  ### categorical covariates => create (ncat-1) columns
  ### also transform covariate.model (may add rows eg for 2 dummy variables) and adjust everything
  if(length(covariate.model)==0) {
    tab<-data.frame(ivar=ivar) 
    Mcovariates<-data.frame(ivar=rep(1,N))
  } else {
    tab<-data.frame(ivar=ivar,data[, rownames(covariate.model),drop=FALSE])
    temp2<-unique(tab)
    temp<-tab[!duplicated(ivar),,drop=FALSE]
    if(dim(temp)[1]!=dim(temp2)[1]) {
      if(verbose) messsage("Some covariates have time-varying values; only the first is taken into account in the current version of the algorithm.\n")
    }
    Mcovariates<-data.frame(ivar=rep(1,N),temp[,2:dim(temp)[2]])
  }
  colnames(Mcovariates)[1]<-name.variable
  
  #  mu.start<-saemix.model@mu.start # initial fixed effects (original parametrization)
  muphi.start<-transpsi(matrix(saemix.model@mu.start,nrow=1),saemix.model["transform.par"]) # starting value for mu on the scale of phi (Gaussian parametrization)
  if(length(covariate.model)>0) {
    beta.start<-covariate.model
    beta.start[covariate.model==1]<- saemix.model@beta.start
    matrixfixedpar.start<-rbind(muphi.start, beta.start)
    fixedpar.start<-matrixfixedpar.start[which(betaest.model>0)]
    fixedpar.start<-matrix(fixedpar.start,ncol=1)
  } else {
    matrixfixedpar.start<-matrix(muphi.start, nrow=1)
    fixedpar.start<-t(matrixfixedpar.start)
  }
  # covariate.estim <-1-covariate.model.fix
  # betas.ini => fixedpar.start
  # ind.covariates => idxmat.fixedpar
  idxmat.fixedpar<-which(saemix.model["betaest.model"]==1)
  nb.fixedpar<-length(idxmat.fixedpar) # total nb of fixed parameters (mu and beta), was nb.betas, should be length(fixedpar.start)
  
  # the covariates in matrix form to compute the vector of parameters as phi^p = mu^p + sum_q beta_q^p cov_q
  ### order: mu_1 for par1, beta_11, beta_12,... for par1,  mu_2 for par2, beta_21, beta_22,... for par2, ...
  ## LCOV, MCOV: nrow=nb of mu+nb of betas; ncol=nb of parameters in the model (=modpar, eg ka, cl)
  ### LCOV: design matrix to pass from the vector of fixed parameters (mu+beta) to the model parameters 
  ### MCOV: values of fixed parameters in matrix form (same structure)
  ## COV: nrow=N, ncol=nb of mu+beta
  ### when a column corresponds to a muphi, it is filled with 1
  ### when a column corresponds to a beta, the values are those of the corresponding covariate for each subject
  ## COV2: nrow=ncol=nb of mu+beta obtained as t(COV).COV (initialised, not used here)
  ## COV1: nrow=N, ncol=nb of columns of COV associated to fixed parameters in the model (???)
  ## dstatCOV: nrow=N, ncol=nb of parameters in the model
  ### formed as the product of COV.MCOV, removing the columns corresponding to fixed parameters
  ### for each parameter, 
  nb.parameters<-saemix.model@npar
  LCOV<-matrix(data=0,nrow=nb.fixedpar,ncol=nb.parameters)
  colnames(LCOV)<-saemix.model@name.modpar
  j1<-1
  COV<-matrix(nrow=dim(Mcovariates)[1],ncol=0)
  name.fixpar <-c()
  #  pfix<-matrix(data=0,nrow=1,ncol=nb.parameters)
  # mean.phi now initialised later
#  mean.phi<-matrix(data=0,nrow=N,ncol=nb.parameters) # population parameters for subject i = mu + sum_q beta_q cov_i,q
  # indx.betaI and indx.betaC (index of mu and of beta) renamed as index.mu and index.beta
  index.mu<-index.beta<-c()
  jpar<-1
  for(j in 1:nb.parameters) {
    name.fixpar<-c(name.fixpar, paste0("muphi.",saemix.model@name.modpar[j]))
    index.mu <- c(index.mu, jpar)
    jpar<-jpar+1
    jcov<-which(betaest.model[,j]==1)
    if(length(jcov)>1) {
      index.beta <- c(index.beta, jpar:(jpar+length(jcov)-2))
      jpar<-jpar+length(jcov)-1
      name.fixpar<-c(name.fixpar, paste0("beta.",saemix.model@name.modpar[j],".",rownames(betaest.model)[jcov[jcov>1]]))
    }
    lambdaj<-matrixfixedpar.start[jcov,j]
    aj<-as.matrix(Mcovariates[,jcov])
    COV<-cbind(COV,aj)
    nlj<-length(lambdaj)
    j2<-j1+nlj-1
    LCOV[j1:j2,j]<-matrix(data=1,nrow=nlj,ncol=1)
    j1<-j2+1
#    if(length(jcov)<=1) mean.phi[,j]<-aj*lambdaj else mean.phi[,j]<-aj%*%lambdaj
    #    pfix[j]<-length(lambdaj)
  }
  #   pfix<-colSums(betaest.model) # ? needed ?
  # indx.betaI and indx.betaC (index of mu and of beta) renamed as index.mu and index.beta
  # idxmat.mu and idxmat.beta
  # index.fixedpar.fix which of these parameters are fixed, see if still needed [was ind.fix1]
  if(length(covariate.model)>0) {
    idxmat.mu<-(c(1:nb.parameters)-1)*dim(betaest.model)[1]+1
    idxmat.beta<-which(betaest.model==1)
    idxmat.beta<-setdiff(idxmat.beta, idxmat.mu)
  } else {
    idxmat.mu<-c(1:nb.parameters)
    idxmat.beta<-c()
  }
  rownames(LCOV)<-name.fixpar
  rownames(COV)<-1:dim(COV)[1]
  MCOV<-LCOV
  
  COV2<-t(COV)%*%COV
  # j.covariate => idxmat.mcov.par: index of elements in LCOV present in the model
  idxmat.mcov.par<-which(LCOV==1)
  MCOV[idxmat.mcov.par]<-fixedpar.start
  betas<-fixedpar.start
  
  # New in 4.0: covariate parameters may be fixed, so we now need to define estimated and non-estimated parameters
  # TODO: check if impact on COV1...
  betaest.model.fix <- rbind(saemix.model@mu.fix, covariate.model.fix)
  index.fixedpar.fix<-which(betaest.model.fix[idxmat.fixedpar]==1)  # fixed fixed effects (elements of mu+beta fixed to their starting value) [ind.fix0]
  index.fixedpar.estim<-which(betaest.model.fix[idxmat.fixedpar]==0) # estimated fixed effects (estimated elements of mu+beta) [ind.fix1]
  COV1 <- COV[,index.fixedpar.estim,drop=FALSE]
  dstatCOV<-COV[,index.fixedpar.fix,drop=FALSE]%*%MCOV[index.fixedpar.fix,,drop=FALSE]
  
  # Index of fixed parameters (mu+beta) on parameters with IIV that are estimated
  # ind.fix11 => index.fixedpariiv.estim (unused for the moment)
#  covariate.estim<-betaest.model*(1-betaest.model.fix)
  fixedpar.model <- betaest.model*(1-betaest.model.fix)
  betaest.model.estimnoIIV <- betaest.model.estimIIV <- fixedpar.model
  if(length(index.omega.novar)>0) betaest.model.estimIIV[,index.omega.novar]<-0
  index.fixedpariiv.estim<-which(betaest.model.estimIIV[idxmat.fixedpar]==1)
#  index.fixedpariiv.estim<-which(betaest.model.estimIIV[idxmat.fixedpar]==0)
  
  # Index of fixed parameters (mu+beta) on parameters without IIV that are estimated
  # ind.fix10 => index.fixedparnoiiv.estim (index in 1:nb(mu+beta))
  # j0.covariate => idxmat.mcov.fixedpar.optim (index in the matrix MCOV (or LCOV))
  if(length(index.omega.var)>0) betaest.model.estimnoIIV[,index.omega.var]<-0
  index.fixedparnoiiv.estim<-which(betaest.model.estimnoIIV[idxmat.fixedpar]==1)
#  index.fixedparnoiiv.estim<-which(betaest.model.estimnoIIV[idxmat.fixedpar]==0)
  
  MCOV0<-MCOV[index.fixedparnoiiv.estim,index.omega.novar,drop=FALSE]
  COV0<-COV[,index.fixedparnoiiv.estim, drop=FALSE]
  idxmat.mcov.fixedpar.optim <- which(LCOV[index.fixedparnoiiv.estim, index.omega.novar]==1)
  # flag.fmin is set to 1 if at least one mu in the model has no IIV and is estimated
  flag.fmin <- as.integer(sum(betaest.model.estimnoIIV[1,])>0)
  
  # Which parameters have variability
  if(length(index.omega.novar)>0) {
    xmat<-fixedpar.model[,index.omega.novar,drop=FALSE]
#    if(is.null(dim(xmat))) xmat<-matrix(xmat,ncol=length(index.omega.novar))
    i0.temp<-which(colSums(xmat)==0)
    ind0.eta<-index.omega.novar[i0.temp] # ind0.eta: index of parameters without IIV or covariates
  } else ind0.eta<-c()
  if(length(ind0.eta)>0) { # ind.eta: index of parameters with IIV
    idx<-1:nb.parameters
    ind.eta<-idx[-c(ind0.eta)]
  } else ind.eta<-1:nb.parameters
  nb.etas<-length(ind.eta)

  # WIP - 21/06/22 to remove eventually (here to check for compatibility with saemix 3.0)
  # Uargs<-list(nchains=saemix.options$nb.chains,nb.parameters=nb.parameters, nb.betas=nb.betas, nb.etas=nb.etas, 
  #             nb.parest=nb.parest,indx.betaC=indx.betaC, indx.betaI=indx.betaI, ind.res=ind.res,indest.omega=indest.omega,
  #             i0.omega2=i0.omega2, i1.omega2=i1.omega2,j.covariate=j.covariate, j0.covariate=j0.covariate,
  #             ind.fix10=ind.fix10, ind.fix11=ind.fix11, ind.fix1=ind.fix1, ind.fix0=ind.fix0,
  #             MCOV0=MCOV0, COV=COV, COV0=COV0, COV1=COV1, LCOV=LCOV, COV2=COV2, dstatCOV=dstatCOV, 
  #             Mcovariates=Mcovariates, ind.ioM=ind.ioM)
#  nb.parest<-sum(covariate.estim)+ sum( saemix.model@var.model[[1]]@omega.model[upper.tri( saemix.model@var.model[[1]]@omega.model, diag=TRUE)])+length(pres) # never used
  
  Uargs<-list(nb.parameters=nb.parameters, nb.betas=nb.fixedpar, nb.etas=nb.etas, 
              indx.betaC=as.integer(index.beta), indx.betaI=as.integer(index.mu), 
              indest.omega=idxmat.omega, i0.omega2=index.omega.novar, i1.omega2=index.omega.var, 
              j.covariate=idxmat.mcov.par, j0.covariate=idxmat.mcov.fixedpar.optim,
              ind.fix10=index.fixedparnoiiv.estim, ind.fix11=index.fixedpariiv.estim, 
              ind.fix1=index.fixedpar.estim, ind.fix0=index.fixedpar.fix,
              MCOV0=MCOV0, COV=COV, COV0=COV0, COV1=COV1, LCOV=LCOV, COV2=COV2, dstatCOV=dstatCOV, 
              Mcovariates=Mcovariates) #, nb.parest=nb.parest, ind.res=ind.res, nchains=saemix.options$nb.chains, ind.ioM=ind.ioM)
  
  # WIP - 21/06/2022
  xvar <- new(Class="SaemixIndivModel", var.model=var.model, covariate.model=covariate.model, covariate.model.fix=covariate.model.fix)
  xvar@Mcovariates <- as.matrix(Mcovariates)
  xvar@N <- nrow(Mcovariates)
  xvar@COV <- COV
  xvar@LCOV <- LCOV
  xvar@MCOV <- MCOV
  xvar@COV0 <- COV0
  xvar@MCOV0 <- MCOV0
  xvar@COV2 <- COV2
  xvar@dstatCOV <- dstatCOV
  xvar@fixedpar <- c(fixedpar.start) # (mu, beta)
  xvar@name.fixedpar <- name.fixpar # names of (mu, beta)
  xvar@nb.fixedpar <- nb.fixedpar # nb(mu, beta)
  if(FALSE) { # To remove, already filled in when calling the class constructor
    xvar@nb.modpar <- nb.parameters # nb of model parameters
    xvar@index.omega.var <- index.omega.var ## i1.omega2
    xvar@index.omega.novar <- index.omega.novar ## i0.omega2
  }
  xvar@index.eta <- ind.eta
  xvar@nb.etas <- nb.etas
  xvar@index.fixedpar.fix <- index.fixedpar.fix ## ind.fix1
  xvar@index.fixedpar.estim <- index.fixedpar.estim  ## ind.fix0
  xvar@index.fixedpariiv.estim <- index.fixedpariiv.estim ## ind.fix11
  xvar@index.fixedparnoiiv.estim <- index.fixedparnoiiv.estim ## ind.fix10
  xvar@idxmat.mcov.fixedpar.optim <- idxmat.mcov.fixedpar.optim ## j0.covariate
  xvar@idxmat.mu <- as.integer(idxmat.mu)
  xvar@idxmat.beta <- as.integer(idxmat.beta)
  xvar@index.mu <- as.integer(index.mu) # indx.betaI
  xvar@index.beta <- as.integer(index.beta) # indx.betaC
  
# needed ?
  # xvar@betaest.model <- betaest.model
  # xvar@betaest.model.fix <- betaest.model.fix
  
  return(list(indiv.model = xvar, flag.fmin=flag.fmin, Uargs=Uargs))
#  return(list(indiv.model = xvar, flag.fmin=flag.fmin, Uargs=Uargs, mean.phi=mean.phi))
}

# call:
# initialisePhi(saemix.data, saemix.model, mean.phi=mean.phi, nb.chains=saemix.options$nb.chains, verbose=saemix.options$warnings)

simulatePhiIOV <- function(saemix.data, saemix.model, nb.chains=1, verbose=FALSE) {
  eta<-NULL
  for(ivarlev in 1:saemix.model@nvarlevel) {
    etaMc<-matrix(data=0, nrow=saemix.model@ind.model[[ivarlev]]@N * nb.chains, ncol=saemix.model@npar)
    etaMc[,saemix.model@ind.model[[ivarlev]]@index.omega.var]<- -0.5*matrix(rnorm(saemix.model@ind.model[[ivarlev]]@N*saemix.model@ind.model[[ivarlev]]@nb.etas), ncol=saemix.model@ind.model[[ivarlev]]@nb.etas)%*% saemix.model@ind.model[[ivarlev]]@chol.omega
    eta<-cbind(eta,
               repmatrix(etaMc, saemix.model@nphirep[[ivarlev]]))
  }
}

initialisePhi <- function(saemix.data, saemix.model, nb.chains=1, verbose=FALSE) {
  # using several Markov chains
  chdat<-new(Class="SaemixRepData",data=saemix.data, nb.chains=nb.chains)
  NM<-chdat["NM"]
  IdM<-chdat["dataM"]$IdM
  yM<-chdat["dataM"]$yM
  XM<-chdat["dataM"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
  # structure of data: DYF, ind.ioM, is.lpdfM
  io<-matrix(data=0,nrow=saemix.data["N"],ncol=max(saemix.data["nind.obs"]))
  for(i in 1:saemix.data["N"])
    io[i,1:saemix.data["nind.obs"][i]]<-1
  ioM<-do.call(rbind,rep(list(io),nb.chains))
  ind.ioM <- which(t(ioM)!=0)
  DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])
  is.lpdf<-rep(0,saemix.data["ntot.obs"])
  which.lpdf<-which(saemix.data@outcome!="continuous")
  if(length(which.lpdf)>0) 
    is.lpdf[saemix.data@data$ytype %in% which.lpdf]<-1
  is.lpdfM<-rep(is.lpdf,nb.chains)
  
  itest.phi<-1:NM
  ltest.phi<-length(itest.phi)
  nvarlevel<-length(saemix.model@ind.model)
  # Population parameters as a table with nb.modpar columns for each level of variability
  mean.phi<-NULL
  for(ivarlev in 1:nvarlevel) {
    mean.phi.var<-indiv.model@COV %*% indiv.model@MCOV
    mean.phi<-cbind(mean.phi,
                    repmatrix(mean.phi.var,saemix.model@nphirep[[ivarlev]]))
  }
  mean.phiM<-do.call(rbind,rep(list(mean.phi),nb.chains))
    
  # Find a valid set of parameters wrt to the structural.model
  # Any parameter set that does not generate NaN, inf or imaginary numbers will satisfy this criteria.
  kt<-0
  NM<-dim(mean.phiM)[1]
  itest.phi<-1:NM
  ltest.phi<-NM
  phiMc<-mean.phiM
  phiM<-matrix(data=0,nrow=NM,ncol=saemix.model@npar)
  while (ltest.phi>0) {
    kt<-kt+1
    if (kt==100) 
      stop("stats:fit.saemix:FailedInitialParameterGuess\nFailed to find a valid initial parameter guess\n")
#    end   
    phiMk<-matrix(data=0,nrow=NM,ncol=saemix.model@npar)
    for(ivarlev in 1:length(saemix.model@ind.model)) {
      i1<-1+(ivarlev-1)*saemix.model@npar
      ## Eco: faux... vérifier dans quel ordre on simule des eta et à quel niveau (pbt par chaine donc faire ça au bon niveau)
      etaMc<-0.5*matrix(rnorm(NM*saemix.model@ind.model[[ivarlev]]@nb.etas),ncol=saemix.model@ind.model[[ivarlev]]@nb.etas)%*%lchol.omega[[ivarlev]]
      etaMc<-repmatrix(etaMc,saemix.model@nphirep[[ivarlev]])
      phiMc[,(i1+saemix.model@ind.model[[ivarlev]]@index.omega.var)]<-mean.phiM[,(i1+saemix.model@ind.model[[ivarlev]]@index.omega.var)]+etaMc
      phiMk<-phiMk+phiMc[,i1:(i1+saemix.model@npar)]
    }
    
    phiM[itest.phi,]<-phiMk[itest.phi,]
    psiM<-transphi(phiM,saemix.model["transform.par"])
    fpred<-saemix.model@model(psiM, IdM, XM)
    inan<-(is.na(fpred)+is.infinite(fpred)+(Im(fpred)!=0))
    itest.phi<-unique(IdM[inan])
    ltest.phi<-length(itest.phi)
  }
  
  # removed from Dargs: N, error.model, etype.exp, modeltype
  ## error.model replaced by list of outcome (used to get the corresponding error models)
  ## modeltype replaced by is.lpdfM (1 if corresponding ytype is non-continuous, 0 for gaussian outcomes associated with an error model)
  ## check if N needed (in this case N would depend on the variability level, N for id, N*nocc for 2nd level, etc...) and etype.exp (probably not, need to transform using the error model for each relevant outcome) 
  # added to Dargs: nb.chains, index.phiM.novar (index.omega.novar for the different submodels, incremented by nb.modpar each time)
  # TODO: test with an exponential error model

  Dargs<-list(IdM=IdM, XM=XM, yM=yM, NM=NM, ind.ioM=ind.ioM, nobs=saemix.data["ntot.obs"], nb.chains=nb.chains,
              yobs=saemix.data["data"][,saemix.data["name.response"]],transform.par=saemix.model["transform.par"],
#              transform.par=rep(saemix.model["transform.par"],saemix.model@nvarlevel),
              model=saemix.model@model, outcome=saemix.model@outcome, is.lpdfM=is.lpdfM)

  return(list(phiM=phiM, saemix.model=saemix.model, Dargs=Dargs, DYF=DYF))
}

repmatrix <- function(mat, times) {
  matrep <- NULL
  for(icol in 1:ncol(mat))
    matrep<-cbind(matrep, rep(mat[,icol], times=times))
  return(matrep)
  #  return(as.data.frame(matrep))
}
