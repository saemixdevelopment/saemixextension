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
  contains = "SaemixIndivModel",
  representation=representation(
    description="character",	# model description
    model="function", 		# name of model function
    sim.model="function", 		# name of function used to simulate from data (used for non-Gaussian models)
    noutcome="integer", # number of responses in the model (defaults to 1 or length(outcome))
    outcome="list" # list of outcomes in the model (of class SaemixOutcome, either discrete SaemixDiscreteOutcome or continuous SaemixContinuousOutcome)
  ),
  validity=function(object){
    return(TRUE)
  }
)

#' @rdname initialize-methods
#' 
#' @param parameters a list of SaemixParameter objects used to create the statistical model
#' @param description a character string, giving a brief description of the model or the analysis
#' @param log a character string, giving the type of the model for the analysis (one of "structural" or "likelihood", defaults to structural)
#' @param model name of the function used to compute the structural model. The
#' function should return a vector of predicted values given a matrix of
#' individual parameters, a vector of indices specifying which records belong
#' to a given individual, and a matrix of dependent variables (see example
#' below).
#' @param indivmodel a list of objects of class SaemixIndivModel
#' @param sim.model if the model contains non-Gaussian outcome, the name of the function used to
#' simulate data from the model (optional, required for diagnotic plots)
#' @param outcome a list of objects of class SaemixOutcome. If missing, one outcome
#' 
#' @exportMethod initialize

setMethod(
  f="initialize",
  signature="SaemixModel",
  definition=function(.Object, parameters, model, description, outcome, verbose=TRUE){
#    cat ("--- initialising SaemixModel Object --- \n")
    if(missing(model)) {
      if(verbose) cat("Error initialising SaemixModel object:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
      return(.Object)
    }
    if(missing(description)) description<-""
    if(missing(outcome)) outcome<-list(y=new(Class="SaemixContinuousOutcome")) else {
      if(inherits(outcome,"character")) { # if only name is given assume continuous default outcome
        y1<-vector(mode="list", length=length(outcome))
        name(y1)<-outcome
        for(i in 1:length(outcome)) y1[[i]]<-list(new(Class="SaemixContinuousOutcome"))
        outcome<-y1
      }
      if(inherits(outcome,"SaemixOutcome")) {
        outcome<-list(new(Class="SaemixContinuousOutcome"))
        name(outcome)<-outcome[[1]]@name.outcome
      }
      if(is(outcome,"list")) { # check type
        is.ok<-0
        for(i in outcome) {
          if(!is(i,"SaemixOutcome")) is.ok<-1
        }
        if(is.ok==1) {
          if(verbose) cat("Error initialising SaemixModel object:\n   Valid outcome values are eithera list of SaemixOutcome objects or a vector of names. Please see the documentation for examples.\n")
          return(.Object)
        }
      }
      # outcome.type <- c()
      # for(i in 1:length(outcome)) {
      #   outcome.type <- c(outcome.type, outcome[[i]]@type.outcome)
      # }
      # if(sum(outcome.type=="continuous")!=length(outcome.type)) {
      #   if(missing(sim.model)) {
      #     logmsg<"There are non-Gaussian outcomes in the model. Diagnostic graphs will not be created for these outcomes unless a simulation function (sim.model) with the same arguments as the model function is added.\n"
      #     if(verbose) cat(logmsg)
      #     .Object@log <- paste0(.Object@log,logmsg)
      #   } else {
      #     if(!inherits(sim.model,"function")) {
      #       logmsg<"If given, the sim.model function must have the same arguments as the model function. Ignoring\n"
      #       if(verbose) cat(logmsg)
      #       .Object@log <- paste0(.Object@log,logmsg)
      #     } else .Object@sim.model<-sim.model
      #   }
      # }
    }
    .Object <- callNextMethod(.Object,parameters=parameters, verbose=verbose)
    
    # May need to add formal checks to model and sim.model (checking arguments)
    .Object@model<-model
    .Object@outcome<-outcome
    .Object@noutcome <- length(outcome)
    .Object@description<-description
			return(.Object)
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
  signature = "SaemixIndivModel" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "log"={return(x@log)},
            "nphi"={return(x@nphi)},
            "param.names"={return(x@param.names)},
            "distribution"={return(x@distribution)},
            "transform"={return(x@transform)},
            "invtransform"={return(x@invtransform)},
            "varlevel"={return(x@varlevel)},
            "covariate"={return(x@covariate)},
            "popmodel"={return(x@popmodel)},
            "varmodel"={return(x@varmodel)},
            "model"={return(x@model)},
            "sim.model"={return(x@sim.model)},
            "description"={return(x@description)},
            "outcome"={return(x@outcome)},
            "noutcome"={return(x@noutcome)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixIndivModel" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "log"={x@log<-value},
            "nphi"={x@nphi<-value},
            "param.names"={x@param.names<-value},
            "distribution"={x@distribution<-value},
            "transform"={x@transform<-value},
            "invtransform"={x@invtransform<-value},
            "varlevel"={x@varlevel<-value},
            "covariate"={x@covariate<-value},
            "popmodel"={x@popmodel<-value},
            "varmodel"={x@varmodel<-value},
            "description"={x@description<-value},
            "model"={x@model<-value},
            "sim.model"={x@sim.model<-value},
            "outcome"={x@outcome<-value},
            "noutcome"={x@noutcome<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)


