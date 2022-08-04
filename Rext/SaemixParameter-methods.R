
#' @name saemixPar
#' @aliases normalPar lognormalPar logitPar probitPar
#' @aliases convertArg2Parameter createSaemixParameter
#' @example 
#' # Defining a list of parameters for SaemixModel objects
#' # first parameter: only name is given, defaults to a parameter with a lognormal distribution
#' # second parameter: defined as a parameter with a normal distribution, setting options
#' # third parameter: given as name="distribution"
#' lpar <- list("ka",vd=normalPar(mu.start=20), cl="lognormal")
#' lpar 
#' @export saemixPar

########################################################################
# General parameters creation

saemixPar<-function(name="theta", distribution="lognormal", estimated=TRUE, prior=FALSE, mu.start=1, omega.start=1, omega.fix=FALSE, omega.level="id", covariate=list(), rho.param=list(), rho=list(), rho.fix=list(), omega=1, fixed=FALSE) {
  # Tests for transformation of vectors to list here
  if(missing(omega.level) & omega==0) omega.level<-c()
  if(missing(estimated) & fixed) estimated<-FALSE
  if(distribution=="lognormal")
    return(lognormalPar(name=name, estimated=estimated, prior=prior, mu.start=mu.start, omega.start=omega.start, omega.fix=omega.fix, omega.level=omega.level, covariate=covariate, rho.param=rho.param, rho=rho, rho.fix=rho.fix))
  if(distribution=="normal")
    return(normalPar(name=name, estimated=estimated, prior=prior, mu.start=mu.start, omega.start=omega.start, omega.fix=omega.fix, omega.level=omega.level, covariate=covariate, rho.param=rho.param, rho=rho, rho.fix=rho.fix))
  if(distribution=="logit")
    return(logitPar(name=name, estimated=estimated, prior=prior, mu.start=mu.start, omega.start=omega.start, omega.fix=omega.fix, omega.level=omega.level, covariate=covariate, rho.param=rho.param, rho=rho, rho.fix=rho.fix))
  if(distribution=="probit")
    return(probitPar(name=name, estimated=estimated, prior=prior, mu.start=mu.start, omega.start=omega.start, omega.fix=omega.fix, omega.level=omega.level, covariate=covariate, rho.param=rho.param, rho=rho, rho.fix=rho.fix))
  new(Class="SaemixParameter", name=name, distribution=distribution, estimated=estimated, prior=prior, mu.start=mu.start, omega.start=omega.start, omega.fix=omega.fix, omega.level=omega.level, covariate=covariate, rho.param=rho.param, rho=rho, rho.fix=rho.fix)
}

########################################################################
# Parameters with different distributions, default values

normalPar<-function(name="", mu.start=1, omega.start=0.5, estimated=TRUE, prior=FALSE, omega.fix=FALSE, omega.level="id", covariate=list(), rho.param=list(), rho=list(), rho.fix=list()) {
  if(missing(omega.start)) omega.start<-0.5*mu.start
  new(Class="SaemixParameter", name=name, distribution="normal", estimated=estimated, prior=prior, mu.start=mu.start, omega.start=omega.start, omega.fix=omega.fix, omega.level=omega.level, covariate=covariate, rho.param=rho.param, rho=rho, rho.fix=rho.fix)
}

lognormalPar<-function(name="", mu.start=1, omega.start=1, estimated=TRUE, prior=FALSE, omega.fix=FALSE, omega.level="id", covariate=list(), rho.param=list(), rho=list(), rho.fix=list()) {
  new(Class="SaemixParameter", name=name, distribution="lognormal",  estimated=estimated, prior=prior, mu.start=mu.start, omega.start=omega.start, omega.fix=omega.fix, omega.level=omega.level, covariate=covariate, rho.param=rho.param, rho=rho, rho.fix=rho.fix)
}

logitPar<-function(name="", mu.start=0.5, omega.start=1, estimated=TRUE, prior=FALSE, omega.fix=FALSE, omega.level="id", covariate=list(), rho.param=list(), rho=list(), rho.fix=list()) {
  new(Class="SaemixParameter", name=name, distribution="logit", estimated=estimated, prior=prior, mu.start=mu.start, omega.start=omega.start, omega.fix=omega.fix, omega.level=omega.level, covariate=covariate, rho.param=rho.param, rho=rho, rho.fix=rho.fix)
}
probitPar<-function(name="", mu.start=0.5, omega.start=1, estimated=TRUE, prior=FALSE, omega.fix=FALSE, omega.level="id", covariate=list(), rho.param=list(), rho=list(), rho.fix=list()) {
  new(Class="SaemixParameter", name=name, distribution="probit", estimated=estimated, prior=prior, mu.start=mu.start, omega.start=omega.start, omega.fix=omega.fix, omega.level=omega.level, covariate=covariate, rho.param=rho.param, rho=rho, rho.fix=rho.fix)
}

################################################
## Function to transform outcomes to a list of SaemixParameter objects
convertArg2Parameter <- function(parameter, verbose=FALSE) {
  if(!is(parameter,"list")) { # parameters are given as names only, names+type, or type only
    saemix.parameter<-vector(mode="list", length=length(parameter))
    if(!is.null(names(parameter))) names(saemix.parameter)<-names(parameter)
    for(i in 1:length(parameter)) {
      if(is(parameter[i],"numeric")) { # given as ka=1
        saemix.parameter[[i]]<-saemixPar(name=names(saemix.parameter)[i], mu.start=parameter[i])
      } else {
        if(parameter[i] %in% c("normal","lognormal","logit","probit")) {
          saemix.parameter[[i]]<-saemixPar(distribution=parameter[i])
          if(is.null(names(saemix.parameter)[i]) || names(saemix.parameter)[i]=="") names(saemix.parameter)[i]<-paste0("out",i)
        } else {
          saemix.parameter[[i]]<-saemixPar()
          names(saemix.parameter)[i]<-parameter[i]
        }
      }
    }
  } else {
    saemix.parameter<-vector(mode="list", length=length(parameter))
    nameout<-paste0("out",1:length(parameter))
    if(!is.null(names(parameter))) names(saemix.parameter)<-names(parameter) else names(saemix.parameter)<-nameout
    names(saemix.parameter)[names(saemix.parameter)==""]<-nameout[names(saemix.parameter)==""]
    for(i in 1:length(parameter)) {
      if(is(parameter[[i]],"SaemixParameter")) {
        saemix.parameter[[i]]<-parameter[[i]]
      } else {
        if(is(parameter[[i]],"numeric")) { # given as ka=1
          saemix.parameter[[i]]<-saemixPar(name=names(saemix.parameter)[i], mu.start=parameter[[i]])
        } else {
          
          if(is(parameter[[i]],"list")) {
            x1<-try(createSaemixParameter(parameter[[i]]))
            if(is(x1, "try-error")) {
              if(verbose) message("Parameter number",i,"of invalid type, assuming a continuous lognormal parameter")
              saemix.parameter[[i]]<-saemixPar(name=names(saemix.parameter)[i])
            } else saemix.parameter[[i]]<-x1
          }
        }
      }
    }
  }
  return(saemix.parameter)
}

createSaemixParameter<-function(paraslist, name="", verbose=FALSE) {
  deflist<-list(distribution="lognormal", estimated=TRUE, prior=FALSE, mu.start=1, omega.start=1, omega.fix=FALSE, omega.level="id")
  namlist<-names(paraslist)
  namdef<-names(deflist)
  for(i in namdef[!(namdef %in% namlist)]) {
    paraslist<-c(paraslist,deflist[[i]])
    names(paraslist)[length(paraslist)]<-i
  }
  if(is.null(paraslist$covariate)) paraslist$covariate <- list()
  if(is.null(paraslist$rho.param)) paraslist$rho.param <- list()
  if(is.null(paraslist$rho)) paraslist$rho <- list()
  if(is.null(paraslist$rho.fix)) paraslist$rho.fix <- list()
  new(Class="SaemixParameter", name=name, distribution=paraslist$distribution, estimated=paraslist$estimated, prior=paraslist$prior, 
      mu.start=paraslist$mu.start, covariate=paraslist$covariate, 
      omega.start=paraslist$omega.start, omega.fix=paraslist$omega.fix, omega.level=paraslist$omega.level, 
      rho.param=paraslist$rho.param, rho=rho, rho.fix=paraslist$rho.fix)
}

################################################
#' Functions for covariance structure
#' 
#' @aliases completeBlocks
#' @export getVarianceModel

getVarianceModel <- function(parameter, level="id", output="var.model") {
  # parameter=list of SaemixParameter objects
  nampar<-names(parameter)
  npar<-length(nampar)
  omega<-diag(x=0,nrow=npar, ncol=npar)
  colnames(omega)<-rownames(omega)<-nampar
  omega.model<-omega.model.fix<-omega
  variable<-level
  name.level<-""
  if(level=="id") name.level<-"iiv"
  if(level=="occ") name.level<-"iov"
  for(ipar in 1:npar) {
    idx<-match(level,parameter[[ipar]]@omega.level)
    if(!is.na(idx)) {
      if(length(parameter[[ipar]]@omega.model)>=idx) omega.model[ipar,ipar]<-parameter[[ipar]]@omega.model[idx] else omega.model[ipar,ipar]<-1
      if(omega.model[ipar,ipar]==1) {
        if(length(parameter[[ipar]]@omega.fix)>=idx) omega.model.fix[ipar,ipar]<-parameter[[ipar]]@omega.fix[idx]
        if(length(parameter[[ipar]]@omega)>=idx) omega[ipar,ipar]<-parameter[[ipar]]@omega[idx] else omega[ipar, ipar]<-ifelse(parameter[[ipar]]@distribution=="lognormal",1,0.5)
      }
      if(length(parameter[[ipar]]@rho.param)>0 && length(parameter[[ipar]]@rho.param[[idx]])>0) { # covariance parameters
        for(i in 1:length(parameter[[ipar]]@rho.param[[idx]])) {
          jpar<-match(parameter[[ipar]]@rho.param[[idx]][i], nampar)
          if(!is.na(jpar)) {
            omega.model[ipar,jpar]<-omega.model[jpar, ipar]<-1
            omega[ipar,jpar]<-omega[jpar, ipar] <- parameter[[ipar]]@rho[[idx]][i]
            omega.model.fix[ipar,jpar]<-omega.model.fix[jpar, ipar]<-parameter[[ipar]]@rho.fix[[idx]][i]
          }
        }
      }
    }
  }
  if(npar>1) {
    for(ipar in 1:(npar-1)) { # compute covariances from the correlations given in rho
      for(jpar in (ipar+1):npar) {
        if(omega[ipar, jpar]!=0) {
          x<-omega[ipar, jpar] * sqrt(omega[ipar,ipar])*sqrt(omega[jpar, jpar])
          omega[ipar, jpar]<-omega[jpar, ipar]<-x
        }
      }
    }
  }
  omega.model<-completeBlocks(omega.model)
  var.model<-saemixVarModel(size=npar, name.level=name.level, variable=variable, omega=omega, omega.model=omega.model, omega.model.fix=omega.model.fix)
  if(!is.null(colnames(var.model@omega))) {
    colnames(var.model@omega.model)<-rownames(var.model@omega.model)<-colnames(var.model@omega)
    colnames(var.model@omega.model.fix)<-rownames(var.model@omega.model.fix)<-colnames(var.model@omega)
  }
  if(output=="list") 
    return(list(omega=omega, omega.model=omega.model, omega.fix=omega.model.fix))
  var.model@omega.names <- paste0("var.",nampar)
  return(var.model)
}

#' Ensure a given covariance matrix has a block structure
#' @export completeBlocks

completeBlocks <- function(xmat) {
  # Check we have a square matrix of 0 and 1
  if(nrow(xmat)!=ncol(xmat) || sum(unique(!(c(xmat) %in% c(0,1)))>0)) return(NULL)
  # Remove covariances for parameters without variance
  for(ipar in 1:nrow(xmat)) {
    if(xmat[ipar,ipar]==0) xmat[,ipar]<-xmat[ipar,]<-0
  }
  # Check if any covariances are present
  if(sum(xmat[lower.tri(xmat)])==0) return(xmat)
  # Check covariance structure, add covariances to make blocks
  lcorr<-vector(mode="list", length=nrow(xmat))
  for(ipar in 1:nrow(xmat)) {
    lcorr[[ipar]]<-which(xmat[ipar,]==1)
  }
  lbloc<-vector(mode="list", length=nrow(xmat))
  lbloc[[1]]<-lcorr[[1]]
  nbloc<-1
  for(ibl in 2:length(lcorr)) {
    # check if we have another block related
    iother<-0
    for(ibl1 in 1:nbloc) {
      if(length(intersect(unlist(lbloc[[ibl1]]), unlist(lcorr[ibl])))>0) {
        lbloc[[ibl1]]<-union(unlist(lbloc[[ibl1]]), unlist(lcorr[ibl]))
        iother<-1
      }
    }
    if(iother==0) {
      nbloc<-nbloc+1
      lbloc[[nbloc]]<-lcorr[[ibl]]
    }
  }
  length(lbloc)<-nbloc
  for(ibl in 1:nbloc) {
    if(length(lbloc[[ibl]])>1) {
      for(ipar in lbloc[[ibl]]) {
        for(jpar in lbloc[[ibl]])
          xmat[ipar,jpar]<-1
      }
    }
  }
  return(xmat)
}

################################################
# Functions for covariate model

#' @aliases removeDuplicateCovDef checkSameCovDef
#' @export getCovariateModel

removeDuplicateCovDef<-function(parameter) {
  nampar<-names(parameter)
  npar<-length(nampar)
  # Sanitise covariate definitions for each parameter - remove duplicated covariate definitions for each parameter
  for(ipar in 1:npar) {
    if(length(parameter[[ipar]]@covariate)>0) {
      idx<-duplicated(parameter[[ipar]]@covariate)
      if(sum(idx)>0) {
        lcov1<-list()
        for(i in which(!idx)) lcov1<-append(lcov1, parameter[[ipar]]@covariate[[i]])
        names(lcov1)<-names(parameter[[ipar]]@covariate)[!(idx)]
        parameter[[ipar]]@covariate<-lcov1
      }
    }
  }
  lcov<-c()
  lcovdef<-list()
  # Sanitise covariates - check if similar definitions but different names, check if same names but different definitions
  for(ipar in 1:npar) {
    if(length(parameter[[ipar]]@covariate)>0) {
      if(length(lcov)==0) {
        lcov<-names(parameter[[ipar]]@covariate)
        lcovdef<-parameter[[ipar]]@covariate
      } else {
        for(icov in 1:length(parameter[[ipar]]@covariate)) {
          isamecov<-checkSameCovDef(parameter[[ipar]]@covariate[[icov]], lcovdef)
          if(isamecov>0) { # same definition, rename if different names
            if(names(parameter[[ipar]]@covariate)[icov] != lcov[isamecov])  names(parameter[[ipar]]@covariate)[icov] <- lcov[isamecov]
          } else {
            lcovdef<-append(lcovdef, parameter[[ipar]]@covariate[[icov]])
            namcov<-names(parameter[[ipar]]@covariate)[icov]
            idx<-match(namcov,lcov)
            if(!is.na(idx)) { # different definitions but same name
              namcov<-paste0(namcov,"2")
              names(parameter[[ipar]]@covariate)[icov]<-namcov
            }
            lcov<-c(lcov, namcov)
            names(lcovdef)[length(lcovdef)]<-namcov
          }
        }
      }
    }
  }
  return(parameter)
}
identicalCovariate<-function(covariate1, covariate2) {
  i0 <- identical(covariate1@covariate.transform, covariate2@covariate.transform, ignore.environment = TRUE) +
    identical(covariate1@name, covariate2@name)
  if(i0<2) return(FALSE) else return(TRUE)
}

checkSameCovDef<-function(covariate, listcov) {
  for(icov in 1:length(listcov)) {
    if(identicalCovariate(covariate, listcov[[icov]])) return(icov)
  }
  return(0)
}

getCovariateModel <- function(parameter) {
  # parameter: list of SaemixParameter objects, already sanitised (check for duplicate names and duplicate covariate definitions through removeDuplicateCovDef)
  nampar<-names(parameter)
  npar<-length(nampar)
  lcov<-c()
  for(ipar in 1:npar) {
    if(length(parameter[[ipar]]@covariate)>0) {
      lcov<-c(lcov, names(parameter[[ipar]]@covariate))
      lcov<-unique(lcov)
    }
  }
  if(length(lcov)>0) {
    covariate.model<-matrix(data=0, nrow=length(lcov), ncol=npar)
    rownames(covariate.model)<-lcov
    colnames(covariate.model)<-nampar
    covariate.model.fix<-covariate.model
    covariates<-vector(mode="list", length=length(lcov))
    names(covariates)<-lcov
    beta.start<-c()
    for(ipar in 1:npar) {
      if(length(parameter[[ipar]]@covariate)>0) {
        for(icov in 1:length(parameter[[ipar]]@covariate)) {
          idx<-match(names(parameter[[ipar]]@covariate)[icov], lcov)
          if(is.null(covariates[[idx]])) covariates[[idx]]<-parameter[[ipar]]@covariate[[icov]] else {
            # check compatibility between covariate definitions
            if(!identicalCovariate(covariates[[idx]],parameter[[ipar]]@covariate[[icov]]))
              message(paste0("Mismatch between covariate definitions for ",lcov))
          }
          covariate.model[idx, ipar]<-1
          if(length(parameter[[ipar]]@covariate[[icov]]@beta.fix)>0 && parameter[[ipar]]@covariate[[icov]]@beta.fix==1) covariate.model.fix[idx, ipar]<-1
          if(length(parameter[[ipar]]@covariate[[icov]]@beta)>0) beta.start<-c(beta.start, parameter[[ipar]]@covariate[[icov]]@beta) else beta.start<-c(beta.start, 0)
        }
      }
    }
  } else {
    covariate.model.fix<-covariate.model<-matrix(nrow=0,ncol=npar)
    colnames(covariate.model.fix)<-colnames(covariate.model)<-nampar
    beta.start<-c()
    covariates<-list()
  }
  return(list(covariates=covariates, covariate.model=covariate.model, covariate.model.fix=covariate.model.fix, beta.start=beta.start))
}


