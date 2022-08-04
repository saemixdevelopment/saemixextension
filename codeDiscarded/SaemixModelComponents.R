# Initial validated part is in:
source("/home/eco/work/saemix/saemixextension/Rext/SaemixModel2.R")

########################################################################
# Parameter model
#' @exportClass SaemixParameterModel

setClass(
  Class="SaemixParameterModel",
  contains = "SaemixStructuralModel",
  representation=representation(
    # model parameters
    nb.parameters="integer",	# nb of parameters in the model
    mu.start="numeric", # initial value for mean population parameters
    mu.fix="numeric",	# 1 when parameter is fixed to its initial value, 0 if estimated [=1-fixed.estim=1-mu.estimated]
    transform.par="numeric",	# distribution for model parameters
    # names
    name.modpar = "character", # name of parameters in the model
    name.betas="character",	# name of all fixed effects (mu + beta_cov)
#    name.sigma="character",	# name of residual parameters (maybe not necessary)
    # covariate model
    covariate.model="matrix",	# covariate model
    covariate.model.fix="matrix",	# fixed parameters in the covariate model (defaults to all 0's)
    transform.covariate="list",	# 
    beta.start="numeric", # initial value for covariate parameters
    betaest.model="matrix",	# 1st line=ones, next lines=covariate model [MAYBE REMOVE]
    # indices (used during the fit)
#    index.res="numeric",		# index of param of residual errors estimated (was indx.res)
    index.mu="numeric",		# index of mean param (mu) estimated (was indx.betaI then indx.fix)
    index.beta="numeric"		# index of covariate param (beta) estimated (was indx.betaC then indx.cov)
  ),
  validity=function(object){
    #    cat ("--- Checking SaemixModel object ---\n")
    return(TRUE)
  }
)

setMethod(
  f="initialize",
  signature="SaemixParameterModel",
  definition=function(.Object, outcome, model, description="", psi0, transform.par, mu.fix, covariate.model, covariate.model.fix, verbose=FALSE){
    # ECO how do we initialize a child using the parent ?
    if(missing(psi0)) {
      if(verbose) message("[ SaemixModel : Error ] Missing initial estimates of parameters (psi0)")
      return(NULL)
    }
    .Object <- callNextMethod(.Object, outcome, model, description, verbose)
    # Model parameters: names, initial values, distribution, fixed parameters
    nb.parameters<-dim(psi0)[2]
    mu.start<-psi0[1,]
    .Object@nb.parameters<-nb.parameters
    if(is.null(colnames(psi0))) name.modpar<-paste0("theta",.Object@nb.parameters) else name.modpar<-colnames(psi0)
    if(missing(transform.par)) transform.par<-c(0, nb.parameters)
    if(length(transform.par)!=nb.parameters) transform.par<-rep(transform.par, length.out=nb.parameters)
    .Object@transform.par<-transform.par
    if(!missing(mu.fix)) {
      if(length(mu.fix)==nb.parameters) .Object@mu.fix<-mu.fix else 
        if(verbose) message("Size mismatch between mu.fix and the number of parameters, ignoring")
    }
    # Covariate models: structure, fixed parameters, initial estimates
    if(!missing(covariate.model)) {
      if(dim(covariate.model)[2]!=nb.parameters) {
        if(verbose) message("Size mismatch between covariate.model and the number of parameters, ignoring")
        covariate.model<-NULL
      } else
        .Object@covariate.model<-covariate.model
    } else covariate.model<-NULL
    if(!missing(covariate.model.fix)) {
      if(dim(covariate.model.fix)[2]!=nb.parameters) {
        if(verbose) message("Size mismatch between covariate.model.fix and the number of parameters, ignoring")
      } else {
        if(is.null(covariate.model)) {
          if(verbose) message("Covariate model missing but fixed covariate parameters given, setting the same structure\n")
          .Object@covariate.model<-covariate.model.fix
        } else {
          .Object@covariate.model.fix<-covariate.model.fix
        }
      }
    }
    if(is.null(covariate.model)) .Object@betaest.model<-matrix(rep(1,dim(psi0)[2]), nrow=1) else {
      if(dim(psi0)[1]==1) beta.start<-rep(0,sum(covariate.model)) else {
        if(dim(psi0)[1]<(dim(covariate.model)[1]+1)) {
          l1<-rep(0,dim(psi0)[2])
          for(irow in (dim(psi0)[1]+1):(dim(covariate.model)[1]+1)) psi0<-rbind(psi0, l1)
        }
        psi1<-psi1[-c(1),,drop=FALSE]
        beta.start<-psi1[covariate.model]
      }
      .Object@beta.start<-beta.start
      .Object@betaest.model<-rbind(rep(1,dim(psi0)[2]),covariate.model)
    }
    .Object@name.betas<-paste0("mu.",.Object@name.modpar)
    validObject(.Object)
    return (.Object)
  }
)

# Getteur
setMethod(
  f ="[",
  signature = "SaemixParameterModel" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "nb.parameters"={return(x@nb.parameters)},
            "mu.start"={return(x@mu.start)},
            "mu.fix"={return(x@mu.fix)},
            "transform.par"={return(x@transform.par)},
            "name.modpar"={return(x@name.modpar)},
            "name.betas"={return(x@name.betas)},
            "covariate.model"={return(x@covariate.model)},
            "covariate.model.fix"={return(x@covariate.model.fix)},
            "beta.start"={return(x@beta.start)},
            "betaest.model"={return(x@betaest.model)},
            "index.mu"={return(x@index.mu)},
            "index.beta"={return(x@index.beta)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixParameterModel" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "nb.parameters"={x@nb.parameters<-value},
            "mu.start"={x@mu.start<-value},
            "mu.fix"={x@mu.fix<-value},
            "transform.par"={x@transform.par<-value},
            "name.modpar"={x@name.modpar<-value},
            "name.betas"={x@name.betas<-value},
            "covariate.model"={x@covariate.model<-value},
            "covariate.model.fix"={x@covariate.model.fix<-value},
            "beta.start"={x@beta.start<-value},
            "betaest.model"={x@betaest.model<-value},
            "index.mu"={x@index.mu<-value},
            "index.beta"={x@index.beta<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)


########################################################################
# Full model, adding variability levels
# Doesn't work, the arguments are not passed from SaemixModel to SaemixStructuralModel for some reason

#' @exportClass SaemixModel


setClass(
  Class="SaemixModel",
  contains = "SaemixParameterModel",
  representation=representation(
    # variability model
    nvarlevel = "numeric", # number of variability levels (currently only 1 is used in the model)
    var.model = "list", # list of variability levels associated to the model parameters
    #    index.iiv="numeric",	# index of random param estimated (was i1.omega2 then indx.omega) [now included in var.model]
    # names
    #    name.random="character",	# name of random parameters
    # elements which will be associated with the data (not-initialised)
    name.predictors="character",# name of predictors 
    name.X="character",	# name of variable to plot on X-axis for graphs
    name.cov="character",	# name of covariates
    Mcovariates="data.frame"	# matrix of individual covariates in the model
  ),
  validity=function(object){
    #    cat ("--- Checking SaemixModel object ---\n")
    # Check validity of var.model (should be a list of SaemixVarLevel, warning if more than 1 level, currently only one treated)
    # if(object@nb.outcome<1 | length(object@outcome)==0) {
    #   message("[ SaemixModel : Error ] Please specify at least one outcome")
    #   return("No outcome given")
    # }
    return(TRUE)
  }
)

#' @rdname initialize-methods
#' 
#' @param var.level can be a vector of strings giving the variability levels, or a list of SaemixVarLevel objects
#' 
#' @exportMethod initialize

# existsMethod("initialize", "SaemixParameterModel")

setMethod(
  f="initialize",
  signature="SaemixModel",
  definition=function(.Object, outcome, model, description="", psi0, beta.start, parameter, transform.par, mu.fix, covariate.model, covariate.model.fix, var.level, verbose=FALSE){
    # if(missing(parameter)) cat("Missing parameter in saemixModel\n")
    # if(missing(model)) cat("Missing model in saemixModel\n")
    # if(missing(outcome)) cat("Missing outcome in saemixModel\n")
    # if(missing(psi0)) cat("Missing psi0 in saemixModel\n")
    .Object <- callNextMethod(.Object, outcome=outcome, model=model, description=description, psi0=psi0, beta.start=beta.start, parameter=parameter, transform.par=transform.par, mu.fix=mu.fix, covariate.model=covariate.model, covariate.model.fix=covariate.model.fix, verbose=verbose)
    #    .Object <- callNextMethod()
    if(!missing(var.level)) {
      if(is.character(var.level)) {   # var.level given as a vector of variability levels
        nvarlevel<-length(var.level)
        var2<-vector(nvarlevel,mode=list)
        for(i in 1:nvarlevel) {
          x<-saemixVarModel(name.level=var.level[i],size=.Object@npar)
          var2[[i]]<-saemixVarNames(x, .Object@name.modpar)
        }
        var.level<-var2
      } else {
        if(is(var.level,"list")) { # var.level given as a list, check elements are all saemixVarLevel objects
          for(i in 1:length(var.level))
            i1<-i1+as.integer(is(var.level[[i]], "SaemixVarLevel"))
          if(i1!=length(var.level)) {
            if(verbose) message("If given, var.level must be either a vector of the variables associated with each variability level, or a list of SaemixVarLevel objects representing the variance structure for the different levels of random effects in the model.")
            var.level<-NULL}
        }
      }
    } else var.level<-NULL
    if(is.null(var.level)) {
      var.level<-saemixVarModel(size=.Object@npar)
      var.level<-list(saemixVarNames(var.level, .Object@name.modpar))
    }
    .Object@var.level<-var.level
    .Object@nvarlevel<-length(var.level)
    validObject(.Object)
    return (.Object)
  }
)


# Getteur
setMethod(
  f ="[",
  signature = "SaemixModel" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "var.model"={return(x@var.model)},
            "nvarlevel"={return(x@nvarlevel)},
            "name.predictors"={return(x@name.predictors)},
            "name.X"={return(x@name.X)},
            "name.cov"={return(x@name.cov)},
            "Mcovariates"={return(x@Mcovariates)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixModel" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "var.model"={x@var.model<-value},
            "nvarlevel"={x@nvarlevel<-value},
            "name.predictors"={x@name.predictors<-value},
            "name.X"={x@name.X<-value},
            "name.cov"={x@name.cov<-value},
            "Mcovariates"={x@Mcovariates<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)


####################################################################################
####  Auxiliary function

#' Matrix diagonal
#' 
#' Extract or replace the diagonal of a matrix, or construct a diagonal matrix (replace diag function from R-base)
#' 
#' @param x	a matrix, vector or 1D array, or missing.
#' @param nrow Optional number of rows for the result when x is not a matrix. 
#' @param ncol Optional number of columns for the result when x is not a matrix. 
#' 
#' @return If x is a matrix then diag(x) returns the diagonal of x. The resulting vector will have names if the matrix x has matching column and rownames.
#' @seealso \code{diag}
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @keywords models
#' @examples
#' 
#' mydiag(1)
#' mydiag(c(1,2))
#' 
#' @export mydiag

# Redefining diag function, too many problems with the R version
mydiag <- function (x = 1, nrow, ncol) {
  if (is.matrix(x)) {
    if (nargs() > 1L) 
      stop("'nrow' or 'ncol' cannot be specified when 'x' is a matrix")
    if ((m <- min(dim(x))) == 0L) 
      return(vector(typeof(x), 0L))
    y <- c(x)[1L + 0L:(m - 1L) * (dim(x)[1L] + 1L)]
    nms <- dimnames(x)
    if (is.list(nms) && !any(sapply(nms, is.null)) && identical((nm <- nms[[1L]][seq_len(m)]), 
                                                                nms[[2L]][seq_len(m)])) 
      names(y) <- nm
    return(y)
  }
  if (is.array(x) && length(dim(x)) != 1L) 
    stop("'x' is an array, but not 1D.")
  if (missing(x)) 
    n <- nrow
  else n <- length(x)
  if (!missing(nrow)) 
    n <- nrow
  if (missing(ncol)) 
    ncol <- n
  p <- ncol
  y <- array(0, c(n, p))
  if ((m <- min(n, p)) > 0L) 
    y[1L + 0L:(m - 1L) * (n + 1L)] <- x
  y
}

############ VALIDITY OF COVARIANCE MODEL 
#' Validate the structure of the covariance model
#' 
#' Check that a matrix corresponds to a structure defining a covariance model for a non-linear mixed effect model.
#' Such a matrix should be composed of only 0s and 1s, with at least one element set to 1, and should be square and symmetrical.
#' 1s on the diagonal indicate that the corresponding parameter has interindividual variability and that its variance will be estimated.
#' 1s as off-diagonal elements indicate that a covariance between the two corresponding parameters will be estimated.
#' 
#' @param x	a matrix
#' @param verbose	a boolean indicating whether warnings should be output if x is not a valid covariance model
#' 
#' @return a boolean, TRUE if x is an acceptable structure and FALSE if not. Messages will be output to describe why x isn't a valid covariance model if the argument verbose is TRUE.
#' @seealso \code{SaemixModel}
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Belhal Karimi
#' @keywords models
#' @examples
#' 
#' covarmodel<-diag(c(1,1,0))
#' validate.covariance.model(covarmodel) # should return TRUE
#' 
#' @export validate.covariance.model
#' 
validate.covariance.model <- function(x, verbose=TRUE){
  #non-square matrix
  if(dim(x)[1]!=dim(x)[2]) {
    if(verbose) message("Error initialising SaemixModel object:\n   The covariance model needs to be a square matrix, please check dimensions.\n")
    return(FALSE)
  }
  # only 0s
  s <- sum(abs(x))
  if(s==0) {
    if(verbose) message("At least one parameter should have IIV in the model, the covariance model may not be only 0s.")
    #   return(FALSE)
  }
  
  #values other than 1 or 0
  s <- sum(x[x!=1 & x!=0])
  if (s>0){
    if(verbose) message("Error initialising SaemixModel object:\n  Invalid covariance model, only 0 or 1 values accepted, please change covariance model.\n")
    return(FALSE)
  }
  
  #asymmetrical
  if (!all(t(x)==x)){
    if(verbose) message("Error initialising SaemixModel object:\n  The matrix defining the covariance model is not symmetrical, please change covariance model.\n")
    return(FALSE)
  }
  
  #values other than 0 when diagonal number is 0
  for (i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if (x[i,j]!=0){
        if(x[i,i]==0 |x[j,j]==0){
          if(verbose) message("Error initialising SaemixModel object:\n  Covariances can only be included between 2 parameters with variances.\n")
          return(FALSE)
        }
      }
    }
  }
  # Check that the matrix has a block structure - doesn't work, fails for simple block :-/
  # indx1<-which(diag(x)==0)
  # if(length(indx1)>0) x1<-x[-c(indx1),-c(indx1)] else x1<-x # removing the lines without variances
  # could maybe work by changing the off-diagonal elements to 0.5 instead of 1...
  # xchol<-try(chol(x1))
  # if(is(xchol, 'try-error')) {
  #   if(verbose) message("Error initialising SaemixModel object:\n  Covariance matrices should be block-diagonal.\n")
  #   return(FALSE)
  # }
  return(TRUE)
}

