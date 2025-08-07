####################################################################################
####			SaemixModel class - definition				####
####################################################################################
#' Class "SaemixModel"
#' 
#' An object of the SaemixModel class, representing a nonlinear mixed-effect
#' model structure, used by the SAEM algorithm.
#' 
#' @name SaemixModel-class
#' @docType class
#' @aliases SaemixModel-class SaemixModel [<-,SaemixModel-method
#' @aliases print,SaemixModel showall,SaemixModel show,SaemixModel summary,SaemixModel 
#' @section Objects from the Class: 
#' An object of the SaemixModel class can be created by using the function \code{\link{saemixModel}} and contain the following slots:
#'   \describe{
#'     \item{\code{model}:}{Object of class \code{"function"}: name of the function used to get predictions from the model (see the User Guide and the online examples for the format and what this function should return).}
#'     \item{\code{description}:}{Object of class \code{"character"}: an optional text description of the model}
#'     \item{\code{psi0}:}{Object of class \code{"matrix"}: a matrix with named columns containing the initial estimates for the parameters in the model (first line) and for the covariate effects (second and subsequent lines, optional). The number of columns should be equal to the number of parameters in the model.}
#'     \item{\code{simulate.function}:}{Object of class \code{"function"}: for non-Gaussian data models, name of the function used to simulate from the model.}
#'     \item{\code{transform.par}:}{Object of class \code{"numeric"}: vector giving the distribution for each model parameter (0: normal, 1: log-normal, 2: logit, 3: probit). Its length should be equal to the number of parameters in the model.}
#'     \item{\code{fixed.estim}:}{Object of class \code{"numeric"}: for each parameter, 0 if the parameter is fixed and 1 if it should be estimated. Defaults to a vector of 1 (all parameters are estimated). Its length should be equal to the number of parameters in the model.}
#'     \item{\code{error.model}:}{Object of class \code{"character"}: name of the error model. Valid choices are "constant" (default), "proportional" and "combined" (see equations in User Guide, except for combined which was changed to y = f + sqrt(a^2+b^2*f^2)*e )}
#'     \item{\code{covariate.model}:}{Object of class \code{"matrix"}: a matrix of 0's and 1's, with a 1 indicating that a parameter-covariate relationship is included in the model (and an associated fixed effect will be estimated). The nmuber of columns should be equal to the number of parameters in the model and the number of rows to the number of covariates.}
#'     \item{\code{covariance.model}:}{Object of class \code{"matrix"}: a matrix f 0's and 1's giving the structure of the variance-covariance matrix. Defaults to the Identity matrix (diagonal IIV, no correlations between parameters)}
#'     \item{\code{omega.init}:}{Object of class \code{"matrix"}: a matrix giving the initial estimate for the variance-covariance matrix}
#'     \item{\code{error.init}:}{Object of class \code{"numeric"}: a vector giving the initial estimate for the parameters of the residual error}
#'   }
#'   Additional elements are added to the model object after a call to \code{saemix} and are used in the algorithm.
#' @section Methods:
#'   \describe{
#'     \item{[<-}{\code{signature(x = "SaemixModel")}: replace elements of object}
#'     \item{[}{\code{signature(x = "SaemixModel")}: access elements of object}
#'     \item{initialize}{\code{signature(.Object = "SaemixModel")}: internal function to initialise object, not to be used}
#'     \item{plot}{\code{signature(x = "SaemixModel")}: plot predictions from the model}
#'     \item{print}{\code{signature(x = "SaemixModel")}: prints details about the object (more extensive than show)}
#'     \item{showall}{\code{signature(object = "SaemixModel")}: shows all the elements in the object}
#'     \item{show}{\code{signature(object = "SaemixModel")}: prints details about the object}
#' 	 }
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' @author Audrey Lavenu
#' @author Marc Lavielle.
#' @seealso \code{\link{SaemixData}} \code{\link{SaemixObject}} \code{\link{saemixControl}} \code{\link{saemix}}
#' \code{\link{plot.saemix}}
#' @keywords classes
#' @exportClass SaemixModel
#' @examples
#' 
#' showClass("SaemixModel")
#' 

# Added:
## outcome
## nb.responses
## omega.fixed: fixed variance parameters (a matrix)

# Add (maybe)
## IOV and associated CI
## iov.fixed: fixed IOV variance parameters (a matrix)
## fixed => set, not estimated (notestim)

### Semantics
## fixed (mu,beta)/random (eta => omega, iov)
## parameters (=psi)/mu (=mu, beta)

# name change
## covariance.model: changed to omega.model
## fixed.estim changed to mu.fixed =1-mu.fixed
## name.fixed changed to name.mu
## name.random changed to name.omega
## name.fixed changed to name.mu

## name.response ? needed ? should it be the name of the outcomes or the name of the response column in data ?

# Restructure/change
## modeltype: maybe just keep to indicate whether multiple responses ?
## variability
## covariate model through formulas
## maybe (TBC) covariance model through formulas (would allow for nested IIV-IOV structure)

# Remove (TBC)
## error.model: now associated with outcome
## error.init: given in outcome
## betaest ? redundant with covariate.model ?

setClass(
  Class="SaemixModel",
  representation=representation(
    model="function", 		# name of model function
    simulate.function="function", 		# name of function used to simulate from data (used for non-Gaussian models)
    description="character",	# model description
    modeltype="character",     # type of model (structural, for continuous responses, likelihood, for discrete responses, combined, when both types are present in the model) => defaults to structural if nb.responses not given or 1
    nb.responses="integer", # number of responses in the model (defaults to 1)
    outcome="list", # list of outcomes in the model (class SaemixOutcome, either discrete SaemixDiscreteOutcome or continuous SaemixContinuousOutcome)
    psi0="matrix",		# CI for parameter estimates
    transform.par="numeric",	# distribution for model parameters
    mu.estimated="numeric",	# 1 when parameter is estimated, 0 if fixed to its initial value
    error.model="character",	# residual error model
    covariate.model="matrix",	# covariate model
    betaest.model="matrix",	# 1st line=ones, next lines=covariate model
    covariance.model="matrix",	# covariance model
    omega.init="matrix",	# CI for Omega
    error.init="numeric",	# CI for residual error
    nb.parameters="integer",	# nb of parameters in the model
    name.modpar="character",	# name of parameters in the model (columns of psi0)
    name.fixed="character",	# name of fixed effects
    name.random="character",	# name of random parameters
    name.sigma="character",	# name of residual parameters (maybe not necessary)
    name.predictors="character",# name of predictors 
    name.X="character",	# name of X 
    name.response="character",	# name of response
    name.cov="character",	# name of covariates
    indx.fix="numeric",		# index of mean param estimated (was indx.betaI)
    indx.cov="numeric",		# index of covariate param estimated (was indx.betaC)
    indx.omega="numeric",	# index of random param estimated (was i1.omega2)
    indx.res="numeric",		# index of param of residual errors estimated (was indx.res)
    Mcovariates="data.frame"	# matrix of individual covariates in the model
  ),
  validity=function(object){
#    cat ("--- Checking SaemixModel object ---\n")
    if (dim(object@psi0)[1]==0) {
      message("[ SaemixModel : Error ] Please provide initial estimates for the fixed effect (a matrix with columns named after the parameters in the model).")
      return("Missing psi0")
    }
    isize<-0
    npar<-dim(object@psi0)[2]
    if(npar!=length(object@transform.par)) isize<-1
    if(npar!=length(object@fixed.estim)) isize<-1
    if (npar!=dim(object@covariate.model)[2]) isize<-1
    if (npar!=dim(object@covariance.model)[1]) isize<-1
    if (npar!=dim(object@omega.init)[1]) isize<-1
#    cat("npar=",npar,length(object@transform.par),length(object@fixed.estim), dim(object@covariate.model)[2],dim(object@covariance.model)[1],dim(object@omega.init)[1],"\n")
    if(isize==1) {
      message("[ SaemixModel : Error ] The number of parameters should be the same in the following elements: psi0 (initial conditions), transform.par, fixed.estim, covariate.model, and the matrices covariance.model and omega.init should be square matrices of size equal to the number of parameters. Please check the input.")
      return("Size mismatch")
    }
    
    # if(npar<2) {
    #   message("[ SaemixModel : Error ] SAEM needs at least two parameters to work on.")
    #   return("Psi0 has size 1")
    # }

		if(sum(object@fixed.estim*mydiag(object@covariance.model))==0) {
			message("[ SaemixModel : Error ] ")
			if(sum(mydiag(object@covariance.model))==0) message("At least one parameter with IIV must be included in the model.") else message("At least one parameter with IIV must be estimated and not fixed in the model.")
			return("Invalid IIV structure")
		}
    if(sum(is.na((match(object@error.model,c('constant','proportional','combined', 'exponential', 'likelihood')))))>0) {
      message("[ SaemixModel : Error ] Invalid residual error model")
      return("Invalid residual error model")
    }
    if(sum(is.na(match(object@modeltype,c("structural","likelihood"))))>0) {
      message("[ SaemixModel : Error ] Invalid type of model (modeltypes should be either structural or likelihood)")
      return("Invalid model type")
    }
    if(!is.null(body(object@simulate.function))) { # Check the simulate.function is formally valid
      simulate.function <- object@simulate.function
      has.sim<-FALSE
      if(!is.function(simulate.function) || length(formals(simulate.function))!=3) {
        message("The simulate.function should have the same format as the model function, ignoring.\n")
        has.sim <- FALSE
      } else has.sim <- TRUE
      if(!has.sim) {
        return("Invalid simulate.function element")
      }
    }
    return(TRUE)
  }
)

### Eco: add vectors when dealing with multiple responses 

#' @rdname initialize-methods
#' 
#' @param model name of the function used to compute the structural model. The
#' function should return a vector of predicted values given a matrix of
#' individual parameters, a vector of indices specifying which records belong
#' to a given individual, and a matrix of dependent variables (see example
#' below).
#' @param description a character string, giving a brief description of the
#' model or the analysis
#' @param modeltype a character string, giving the type of the model for the analysis (one of "structural" or "likelihood", defaults to structural)
#' @param psi0 a matrix with a number of columns equal to the number of
#' parameters in the model, and one (when no covariates are available) or two
#' (when covariates enter the model) giving the initial estimates for the fixed
#' effects. The column names of the matrix should be the names of the
#' parameters in the model, and will be used in the plots and the summaries.
#' When only the estimates of the mean parameters are given, psi0 may be a
#' named vector.
## #' @param name.response a character string or a column number specifying which column of the data contains the dependent variable
#' @param name.sigma a vector of character string giving the names of the residual error parameters (defaults to "a" and "b")
#' @param transform.par the distribution for each parameter (0=normal,
#' 1=log-normal, 2=probit, 3=logit). Defaults to a vector of 1s (all parameters
#' have a log-normal distribution)
#' @param fixed.estim whether parameters should be estimated (1) or fixed to
#' their initial estimate (0). Defaults to a vector of 1s
#' @param error.model type of residual error model (valid types are constant,
#' proportional, combined and exponential). Defaults to constant
#' @param covariate.model a matrix giving the covariate model. Defaults to no
#' covariate in the model
#' @param covariance.model a square matrix of size equal to the number of parameters in the model, 
#' giving the variance-covariance matrix of the model: 1s correspond to estimated variances (in the diagonal) 
#' or covariances (off-diagonal elements). Defaults to the identity matrix
#' @param omega.init a square matrix of size equal to the number of parameters
#' in the model, giving the initial estimate for the variance-covariance matrix
#' of the model. The current default is a diagonal matrix with ones for all
#' transformed parameters as well as for all untransformed parameters with an
#' absolute value smaller than one.  For untransformed parameters greater or
#' equal to one, their squared value is used as the corresponding diagonal
#' element.
#' @param error.init a vector of size 2 giving the initial value of a and b in
#' the error model. Defaults to 1 for each estimated parameter in the error model
#' @param name.modpar names of the model parameters, if they are not given as
#' the column names (or names) of psi0
#' 
#' @exportMethod initialize

setMethod(
  f="initialize",
  signature="SaemixModel",
  definition=function(.Object, model, description, modeltype,psi0, name.response, name.sigma, transform.par,fixed.estim, error.model,covariate.model,covariance.model,omega.init,error.init, name.modpar, verbose=TRUE){
#    cat ("--- initialising SaemixModel Object --- \n")
    if(missing(name.response)) name.response<-""
    if(missing(modeltype)) modeltype<-rep("structural",length(name.response))
    .Object@modeltype<-modeltype
    if(missing(model)) {
#      cat("Error initialising SaemixModel object:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
      return(.Object)
    }
    .Object@model<-model
    if(missing(description)) description<-""
    .Object@description<-description
    if(missing(psi0) || length(psi0)==0) {
      if(verbose) message("Error initialising SaemixModel object:\n   Please provide initial estimates for the fixed effect (a matrix with columns named after the parameters in the model).\n")
      return(.Object)
    }
    npar<-dim(psi0)[2]
    if(missing(name.modpar) || length(name.modpar)==0) {
      y1<-try(name.modpar<-colnames(psi0))
      if(inherits(y1,"try-error")) {
        if(verbose) message("     Can't find parameter names.\n")
        name.modpar<-paste("theta",1:npar)
      }
    }
    if(is.null(colnames(psi0))) {
      y1<-try(colnames(psi0)<-name.modpar)
      if(inherits(y1,"try-error")) {
        if(verbose) message("Warning:\n   Problem with names of psi0\n")
        colnames(psi0)<-name.modpar<-paste("theta",1:npar)
      }
    }
    if(missing(name.response)) name.response<-""
    .Object@name.response<-name.response
    if(!missing(covariate.model)) {
    	if(dim(psi0)[1]<2 & sum(covariate.model)>0){
    		psi0<-rbind(psi0,rep(0,dim(psi0)[2]))
    	}
    }
    if(is.null(rownames(psi0))) {
      rownames(psi0)<-rep("",dim(psi0)[1])
      rownames(psi0)[1]<-"Pop.CondInit"
      if(dim(psi0)[1]>1) rownames(psi0)[2:dim(psi0)[1]]<-"Cov.CondInit"
    }
    .Object@psi0<-psi0    
    .Object@name.modpar<-name.modpar
    if(missing(error.model) || length(error.model)==0) error.model<-"constant"
    length(error.model)<-length(.Object@modeltype)
    error.model[.Object@modeltype=="likelihood"]<-"likelihood"
    error.model[is.na(error.model)]<-"constant"
    if(sum(!error.model %in% c('constant','proportional','combined', 'exponential', 'likelihood'))) {
      message("Invalid error model, switching to constant")
      error.model[!error.model %in% c('constant','proportional','combined', 'exponential','likelihood')] <- "constant"
    }
    # normally not needed now (already adjusted to size of modeltype)
#    if(length(error.model)<length(name.response)) error.model<-rep(error.model,length.out=length(name.response))
    .Object@error.model<-error.model
# Checking sizes
    .Object@nb.parameters<-npar
    if(missing(transform.par) || length(transform.par)==0) transform.par<-rep(0,npar)
    .Object@transform.par<-transform.par
    if(missing(fixed.estim) || length(fixed.estim)==0) fixed.estim<-rep(1,npar)
    .Object@fixed.estim<-fixed.estim
    if(missing(covariate.model) || length(covariate.model)==0 || sum(covariate.model)==0) covariate.model<-matrix(nrow=0,ncol=npar)
    if(is.null(dim(covariate.model)) & length(covariate.model)>0) covariate.model<-matrix(covariate.model,byrow=T,ncol=npar) # Covariate model given as a vector
    if(is.null(colnames(covariate.model))) colnames(covariate.model)<-colnames(psi0)
    .Object@covariate.model<-covariate.model
    if(missing(covariance.model) || length(covariance.model)==0) {
      covariance.model<-diag(nrow=npar,ncol=npar)
    } else {
      if(dim(covariance.model)[1]!=dim(covariance.model)[2]) {
        if(verbose) message("Error initialising SaemixModel object:\n   The covariance model needs to be a square matrix, please check dimensions.\n")
      return(.Object)
      }
    }
    nomg<-dim(covariance.model)[1]
    if(nomg!=npar) {
      if(verbose) message("Error initialising SaemixModel object:\n   The covariance model needs to have the same size as the number of parameters.\n")
      return(.Object)
    }
    if(!validate.covariance.model(covariance.model, verbose)) return(.Object)
    if(is.null(colnames(covariance.model))) colnames(covariance.model)<-rownames(covariance.model)<-colnames(psi0)
    .Object@covariance.model<-covariance.model
    indx.omega<-which(diag(covariance.model)>0)
    .Object@indx.omega<-indx.omega
    if(!missing(omega.init) && length(omega.init)>0) {
      if(dim(omega.init)[1]!=dim(omega.init)[2]) {
        if(verbose) message("Warning:   the matrix giving the initial conditions for the covariance model (omega.init) needs to be a square matrix. Changing it to the diagonal matrix\n")
        omega.init<-NULL
      }
    }
    if(missing(omega.init) || length(omega.init)==0) {
      omega.init<-diag(fixed.estim)
      diag.omegi<-rep(1,npar)
      j1<-which(transform.par==0)
      if(length(j1)>0) {
      	diag.omegi[j1]<-sapply(psi0[1,j1]**2,function(x) { x[x<1]<-1; return(x)})
#      for(i in j1) d[i]<-max(psi0[i]^2,1)
      }
      omega.init<-diag(diag.omegi,nrow=npar)
    }
    if(is.null(colnames(omega.init))) {
      if(dim(omega.init)[1]==length(colnames(psi0))) colnames(omega.init)<-rownames(omega.init)<-colnames(psi0) else message("The dimensions of omega.init don't agree with the number of parameters")
    }
    .Object@omega.init<-omega.init
		if(sum(.Object@fixed.estim*mydiag(.Object@covariance.model))==0) {
			message("Error initialising SaemixModel object:\n")
#			if(sum(mydiag(.Object@covariance.model))==0) cat("At least one parameter with IIV must be included in the model.\n") else cat("At least one parameter with IIV must be estimated and not fixed in the model.\n")
			return(.Object)
		}

## Residual Error model.
# error models are a + bf described by [a b]
# error models :
#   constant            y = f + a*e
#   proportional        y = f + b*f*e
#   combined            y = f + sqrt(a^2+b^2*f^2)*e
#   exponential         y = f*exp(a*e)    ( <=>  log(y) = log(f) + a*e )
    if(missing(error.init) || length(error.init)!=2*length(.Object@error.model)) {
      error.init<-c()
      for(i in 1:length(.Object@error.model)) {
        error.init<-c(error.init,switch(error.model[i],
        "constant"=c(1,0),
        "exponential"=c(1,0),
        "proportional"=c(0,1),
        "combined"=c(1,1),
        "likelihood"=c(0,0)))
     }
    }
    if (any(error.init < 0)) {
      error.init<-abs(error.init)
      if(verbose) message("Initial estimates for error model parameters should be non-negative, changing to absolute value")
    }
    xres<-c()
    if(missing(name.sigma)) mis.sig<-TRUE else mis.sig<-FALSE
    if(missing(name.sigma) || length(name.sigma)!=2) name.sigma<-c("a","b")
    if(!mis.sig) { # & .Object@name.response!="" # pb if response has more than 1 element
      for(i in 1:length(.Object@name.response)) {
        if(.Object@name.response[i]!="") xres<-c(xres,paste(name.sigma,.Object@name.response[i],sep=".")) else xres<-c(xres,paste(name.sigma,i,sep=".")) # xres<-c(xres,name.sigma) ?
      }
    } else xres<-rep(name.sigma,length(.Object@name.response))
    .Object@name.sigma<-xres
    .Object@error.init<-error.init
    indx.res<-c()
    indx.res1<-c()
    for(i in 1:length(.Object@error.model)) {
      if(.Object@error.model[i]=='constant') {
        indx.res1<-1
    } else {
        if(.Object@error.model[i]=='proportional') {
          indx.res1<-2
      } else {
          if(.Object@error.model[i]=='combined') {
            indx.res1<-c(1,2) 
        } else {
            if(.Object@error.model[i]=='exponential') {
              indx.res1<-1
           } else indx.res1<-c()
        }
      }
    }
      if(length(indx.res1)>0) indx.res<-c(indx.res,indx.res1+2*(i-1))
    }
    if(length(indx.res)>0) .Object@error.init[-indx.res]<-0 # if indx.res is c() then only likelihood type responses in the model
    if(length(indx.res)>0) .Object@indx.res<-indx.res
    .Object@betaest.model<-matrix(c(rep(1,.Object@nb.parameters), c(t(.Object@covariate.model))),ncol=.Object@nb.parameters,byrow=TRUE)
    colnames(.Object@betaest.model)<-colnames(.Object@covariate.model)
    if(!is.null(rownames(.Object@covariate.model))) {
      rownames(.Object@betaest.model)<-c("Fixed",rownames(.Object@covariate.model))
    } else {
      rownames(.Object@betaest.model)<-rep("",dim(.Object@betaest.model)[1])
      rownames(.Object@betaest.model)[1]<-"Fixed"
    }
# Object validation
    validObject(.Object)
    return (.Object)
  }
)

####################################################################################
####			saemixModel class - accesseur				####
####################################################################################

#' Get/set methods for SaemixModel object
#' 
#' Access slots of an SaemixModel object using the object\["slot"\] format
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
  signature = "SaemixModel" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "model"={return(x@model)},
    "simulate.function"={return(x@simulate.function)},
    "description"={return(x@description)},
    "modeltype"={return(x@modeltype)},
    "psi0"={return(x@psi0)},
    "transform.par"={return(x@transform.par)},
    "fixed.estim"={return(x@fixed.estim)},
    "error.model"={return(x@error.model)},
    "covariate.model"={return(x@covariate.model)},
    "betaest.model"={return(x@betaest.model)},
    "covariance.model"={return(x@covariance.model)},
    "omega.init"={return(x@omega.init)},
    "error.init"={return(x@error.init)},
    "nb.parameters"={return(x@nb.parameters)},
    "name.modpar"={return(x@name.modpar)},
    "name.fixed"={return(x@name.fixed)},
    "name.random"={return(x@name.random)},
    "name.sigma"={return(x@name.sigma)},
    "name.X"={return(x@name.X)},
    "name.response"={return(x@name.response)},
    "name.predictors"={return(x@name.predictors)},
    "name.cov"={return(x@name.cov)},
    "indx.fix"={return(x@indx.fix)},
    "indx.cov"={return(x@indx.cov)},
    "indx.omega"={return(x@indx.omega)},
    "indx.res"={return(x@indx.res)},
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
    "model"={x@model<-value},
    "simulate.function"={x@simulate.function<-value},
    "description"={return(x@description)},
    "modeltype"={return(x@modeltype)},
    "psi0"={x@psi0<-value},
    "transform.par"={x@transform.par<-value},
    "fixed.estim"={x@fixed.estim<-value},
    "error.model"={x@error.model<-value},
    "covariate.model"={x@covariate.model<-value},
    "betaest.model"={x@betaest.model<-value},
    "covariance.model"={x@covariance.model<-value},
    "omega.init"={x@omega.init<-value},
    "error.init"={x@error.init<-value},
    "nb.parameters"={x@nb.parameters<-value},
    "name.modpar"={x@name.modpar<-value},
    "name.fixed"={x@name.fixed<-value},
    "name.random"={x@name.random<-value},
    "name.sigma"={x@name.sigma<-value},
    "name.X"={x@name.X<-value},
    "name.response"={x@name.response<-value},
    "name.predictors"={x@name.predictors<-value},
    "name.cov"={x@name.cov<-value},
    "indx.fix"={x@indx.fix<-value},
    "indx.cov"={x@indx.cov<-value},
    "indx.omega"={x@indx.omega<-value},
    "indx.res"={x@indx.res<-value},
    "Mcovariates"={x@Mcovariates<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

