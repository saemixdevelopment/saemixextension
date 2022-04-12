# Model structure
#' @exportClass SaemixStructuralModel

# initialize methods - only accept arguments with the proper format

setClass(
  Class="SaemixStructuralModel",
  representation=representation(
    # Outcome
    nb.outcome="integer", # number of outcome in the model (set to length(outcome))
    outcome="list", # list of outcomes in the model (class SaemixOutcome, either discrete SaemixDiscreteOutcome or continuous SaemixContinuousOutcome)
    name.outcome="character",	# name of responses
    # Model structure
    description="character",	# model description
    #    modeltype="character",     # type of model (structural, for continuous responses, likelihood, for discrete responses, combined, when both types are present in the model) => defaults to structural if nb.outcome not given or 1 [MAYBE REMOVE]
    model="function" 		# name of model function
  ),
  validity=function(object){
    #    cat ("--- Checking SaemixModel object ---\n")
    if(length(object@nb.outcome)==0) {
      message("[ SaemixModel : Error ] Please specify at least one outcome")
      return("No outcome given")
    }
    if(!is.function(object@model) || !identical(names(formals(object@model)),c("psi","id","xidep"))) {
      message("[ SaemixModel : Error ] Invalid type of model")
      return("Invalid model type")
    }
    return(TRUE)
  }
)

setMethod(
  f="initialize",
  signature="SaemixStructuralModel",
  definition=function(.Object, outcome, model, description="", verbose=FALSE){
    # if(missing(model)) cat("Missing model in structural model\n")  # else print(model)
    # if(missing(outcome)) cat("Missing outcome in structural model\n") # else print(outcome)
    # 
    if(missing(model)) {
      if(verbose) message("Error initialising SaemixModel object:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
      return(.Object)
    }
    if(missing(outcome)) {
      if(verbose) message("No outcome given, assuming a single continuous outcome with a constant error model")
      outcome<-list(y=new(Class="SaemixContinuousOutcome"))
    }
    if(!is(outcome, "list")) { # if a single outcome is given, transform to a list
      if(!is(outcome, "SaemixOutcome")) {
        if(verbose) message("Outcomes must be given as a list of SaemixOutcome objects")
        return(.Object)
      }
      outcome<-list(outcome)
    }
    .Object@nb.outcome<-length(outcome)
    if(is.null(names(outcome))) names(outcome)<-paste0("y",1:.Object@nb.outcome)
    .Object@outcome<-outcome
    .Object@model<-model
    .Object@name.outcome<-names(outcome)
    .Object@description<-description
    # Object validation
    validObject(.Object)
    return (.Object)
  }
)

#' Get/set methods for SaemixModel object
#' 
#' Access slots of a SaemixModel object using the object\["slot"\] format
#' 
#' @param x object
#' @param i element to be replaced
#' @param j element to replace with
#' @param drop whether to drop unused dimensions
#' @keywords methods
#' @exportMethod [
#' @exportMethod [<-


# Getteur
setMethod(
  f ="[",
  signature = "SaemixStructuralModel" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "outcome"={return(x@outcome)},
            "nb.outcome"={return(x@nb.outcome)},
            "name.outcome"={return(x@name.outcome)},
            "description"={return(x@description)},
            "model"={return(x@model)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixStructuralModel" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "outcome"={x@outcome<-value},
            "nb.outcome"={x@nb.outcome<-value},
            "name.outcome"={x@name.outcome<-value},
            "description"={x@description<-value},
            "model"={x@model<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

####################################################################################
# Parameter model
#' @exportClass SaemixParameterModel

setClass(
  Class="SaemixParameterModel",
  contains = "SaemixStructuralModel",
  representation=representation(
    # model parameters
    npar="integer",	# nb of parameters in the model
    mu.start="numeric", # initial value for mean population parameters
    mu.fix="numeric",	# 1 when parameter is fixed to its initial value, 0 if estimated [=1-fixed.estim=1-mu.estimated]
    transform.par="numeric",	# distribution for model parameters
    # names
    name.modpar = "character", # name of parameters in the model
    #    name.sigma="character",	# name of residual parameters (maybe not necessary)
    # covariate model
    covariate.model="matrix",	# covariate model
    covariate.model.fix="matrix",	# fixed parameters in the covariate model (defaults to all 0's)
    transform.covariate="list",	# [MAYBE REMOVE and associate instead to covariates]
    beta.start="numeric", # initial value for covariate parameters
    betaest.model="matrix"	# 1st line=ones, next lines=covariate model [MAYBE REMOVE]
  ),
  validity=function(object){
    #    cat ("--- Checking SaemixParameterModel object ---\n")
    return(TRUE)
  }
)

# existsMethod("initialize", "SaemixParameterModel")

# Initialize - all arguments with unique type, move ambiguous definitions to the creator function

setMethod(
  f="initialize",
  signature="SaemixParameterModel",
  definition=function(.Object, outcome, model, parameter, description="", mu.start, transform.par, mu.fix, covariate.model, covariate.model.fix, beta.start, verbose=FALSE){
    if(missing(parameter) || !is(parameter, "character")) {
      if(verbose) message("Please give the name of model parameters")
      return(.Object)
    }
    .Object <- callNextMethod(.Object, outcome=outcome, model=model, description=description, verbose=verbose)
    # Model parameters: names, initial values, distribution, fixed parameters
    npar<-length(parameter)
    if(missing(mu.start)) {
      mu.start<-rep(1,npar)
    } else {
      if(length(mu.start)!=npar & verbose) message("Size mismatch, resizing mu.start")
      length(mu.start)<-npar
      mu.start[is.na(mu.start)]<-1
    }
    if(missing(transform.par)) {
      transform.par<-rep(0,npar)
    } else {
      if(length(transform.par)!=npar & verbose) message("Size mismatch, resizing transform.par")
      length(transform.par)<-npar
      transform.par[is.na(transform.par)]<-0
    }
    if(missing(mu.fix)) {
      mu.fix<-rep(0,npar)
    } else {
      if(length(mu.fix)!=npar & verbose) message("Size mismatch, resizing mu.fix")
      length(mu.fix)<-npar
      mu.fix[is.na(mu.fix)]<-0
    }
    .Object@name.modpar <- parameter
    .Object@npar<-npar
    .Object@transform.par<-transform.par
    .Object@mu.start<-mu.start
    .Object@mu.fix<-mu.fix
    
    # Covariate model
    if(missing(covariate.model)) covariate.model<-NULL
    if(missing(covariate.model.fix)) covariate.model.fix<-NULL
    if(missing(beta.start)) beta.start<-NULL
    .Object <- setCovariateModel(.Object, covariate.model=covariate.model, covariate.model.fix=covariate.model.fix, beta.start=beta.start, verbose=verbose) 

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
            "npar"={return(x@npar)},
            "mu.start"={return(x@mu.start)},
            "mu.fix"={return(x@mu.fix)},
            "transform.par"={return(x@transform.par)},
            "name.modpar"={return(x@name.modpar)},
            "covariate.model"={return(x@covariate.model)},
            "covariate.model.fix"={return(x@covariate.model.fix)},
            "beta.start"={return(x@beta.start)},
            "betaest.model"={return(x@betaest.model)},
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
            "npar"={x@npar<-value},
            "mu.start"={x@mu.start<-value},
            "mu.fix"={x@mu.fix<-value},
            "transform.par"={x@transform.par<-value},
            "name.modpar"={x@name.modpar<-value},
            "covariate.model"={x@covariate.model<-value},
            "covariate.model.fix"={x@covariate.model.fix<-value},
            "beta.start"={x@beta.start<-value},
            "betaest.model"={x@betaest.model<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

# Auxiliary function to set the covariate model and initial estimates
setCovariateModel <- function(.Object, covariate.model=NULL, covariate.model.fix=NULL, beta.start=NULL, verbose=FALSE) {
  # Covariate models: structure, fixed parameters, initial estimates
  npar<-.Object@npar
  if(!is.null(covariate.model)) {
    if(dim(covariate.model)[2]!=npar) {
      if(verbose) message("Size mismatch between covariate.model and the number of parameters, ignoring")
      covariate.model<-NULL
    } else {
      colnames(covariate.model)<-.Object@name.modpar
      .Object@covariate.model<-covariate.model
    }
    
  } else covariate.model<-NULL
  if(!is.null(covariate.model.fix)) {
    if(dim(covariate.model.fix)[2]!=npar) {
      if(verbose) message("Size mismatch between covariate.model.fix and the number of parameters, ignoring")
    } else {
      colnames(covariate.model.fix)<-.Object@name.modpar
      if(is.null(covariate.model)) {
        if(verbose) message("Covariate model missing but fixed covariate parameters given, setting the same structure\n")
        .Object@covariate.model<-covariate.model.fix
      } else {
        .Object@covariate.model.fix<-covariate.model.fix
      }
    }
  }
  if(is.null(covariate.model)) {
    .Object@betaest.model<-matrix(rep(1,npar), nrow=1, dimnames=list(c("mu.psi"),.Object@name.modpar))
    beta.start<-numeric()
  } else {
    .Object@betaest.model<-rbind(rep(1,npar),covariate.model)
    colnames(.Object@betaest.model)<-.Object@name.modpar
    nb.beta<-sum(covariate.model)
    if(missing(beta.start)) beta.start<-rep(0,nb.beta)
    if(nb.beta>length(beta.start)) 
      beta.start<-c(beta.start, rep(0,nb.beta-length(beta.start)))
    beta.start<-beta.start[1:nb.beta] # if beta.start longer than the number of par-cov relationships
  }
  .Object@beta.start<-beta.start
  return(.Object)
}

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
    name.thetas="character",	# name of all fixed effects (mu + beta_cov) [filled later on when associated with data]
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
#' @param var.model can be a vector of strings giving the variability levels, or a list of SaemixVarLevel objects
#' 
#' @exportMethod initialize

# existsMethod("initialize", "SaemixParameterModel")

setMethod(
  f="initialize",
  signature="SaemixModel",
  definition=function(.Object, outcome, model, parameter, description="", mu.start, beta.start, transform.par, mu.fix, covariate.model, covariate.model.fix, var.model, verbose=FALSE){
    # if(missing(parameter)) cat("Missing parameter in saemixModel\n")
    # if(missing(model)) cat("Missing model in saemixModel\n")
    # if(missing(outcome)) cat("Missing outcome in saemixModel\n")
    # if(missing(mu.start)) cat("Missing mu.start in saemixModel\n")
    .Object <- callNextMethod(.Object, outcome=outcome, model=model, description=description, mu.start=mu.start, beta.start=beta.start, parameter=parameter, transform.par=transform.par, mu.fix=mu.fix, covariate.model=covariate.model, covariate.model.fix=covariate.model.fix, verbose=verbose)
    # Setting variability model(s) (for the moment only one used in the actual estimation)
    if(!missing(var.model)) {
      if(is(var.model,"SaemixVarLevel")) {
        var.model<-list(var.model) 
        } else {  # var.model given as a vector of variability levels
        ierr<-0
        if(is(var.model,"list")) { # check all items are valid SaemixVarLevel objects
          ierr<-0
          for(i in 1:length(var.model))
            ierr<-ierr+as.integer(!is(var.model[[i]], "SaemixVarLevel"))
        } else ierr<-1
        if(ierr>0) {
          if(verbose) message("If given, var.model must be either a vector of the variables associated with each variability level, or a list of SaemixVarLevel objects representing the variance structure for the different levels of random effects in the model.")
          var.model<-NULL
        } 
      }
      # if(is.character(var.model)) {   # var.model given as a vector of variability levels
      #   nvarlevel<-length(var.model)
      #   var2<-vector(nvarlevel,mode=list)
      #   for(i in 1:nvarlevel) {
      #     x<-saemixVarModel(name.level=var.model[i],size=.Object@npar)
      #     var2[[i]]<-saemixVarNames(x, .Object@name.modpar)
      #   }
      #   var.model<-var2
      # }
    } else var.model<-NULL
    if(is.null(var.model)) {
      var.model<-saemixVarModel(size=.Object@npar)
      var.model<-list(iiv=saemixVarNames(var.model, .Object@name.modpar))
    } else { # add names to omega parameters if missing
      namvar<-c()
      for(i in 1:length(var.model)) {
        if(length(var.model[[i]]@omega.names)==0) var.model[[i]]<-saemixVarNames(var.model[[i]])
        if(length(var.model[[i]]@name.level)==0) var.model[[i]]@name.level<-paste0("var",i)
        namvar<-c(namvar, var.model[[i]]@name.level)
      }
      names(var.model)<-namvar
    }
    nvarlevel<-length(var.model)
    if(verbose & nvarlevel>1) message("Currently only the first level of variability is taken into account in saemix")
    .Object@var.model<-var.model
    .Object@nvarlevel<-length(var.model)
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
            "name.thetas"={return(x@name.thetas)},
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
            "name.thetas"={x@name.thetas<-value},
            "Mcovariates"={x@Mcovariates<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

########################################################  Auxiliary functions
# Redefining diag function, too many problems with the R version

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

# Testing the Validity of covariance model

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

# TODO: a function to make sure levels of variability are nested and compatible with one another
