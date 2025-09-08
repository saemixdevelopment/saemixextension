#' @include aaa_generics.R
#' @include SaemixData.R
#' @include SaemixData-methods.R
#' @include SaemixData-methods_covariates.R
#' @include SaemixOutcome.R
#' @include SaemixParameter.R
#' @include SaemixVarModel.R
#' @include SaemixPopModel.R
#' @include SaemixIndivModel.R
NULL

# ToDo:
## finish description of SaemixIndivModel and copy the slots here
## check description for the slots of the child class
# showall class (currently defaults to print)

####################################################################################
####			SaemixModel class - definition				####
####################################################################################
#' Class "SaemixModel"
#' 
#' An object of the SaemixModel class, representing a nonlinear mixed-effect
#' model structure, used by the SAEM algorithm.
#' 
#' @name SaemixModel-class
#' 
#' @docType class
#' 
#' @aliases SaemixModel-class SaemixModel [<-,SaemixModel-method
#' @aliases print,SaemixModel showall,SaemixModel show,SaemixModel summary,SaemixModel 
#' 
#' @section Objects from the Class: 
#' An object of the SaemixModel class can be created by using the function \code{\link{saemixModel}} and contain the following slots:
#'   \describe{
#'     \item{\code{description}:}{Object of class \code{"character"}: an optional text description of the model}
#'     \item{\code{model}:}{Object of class \code{"function"}: name of the function used to get predictions from the model (see the User Guide and the online examples for the format and what this function should return).}
#'     \item{\code{sim.model}:}{Object of class \code{"function"}: for non-Gaussian data models, name of the function used to simulate from the model.}
#'     \item{\code{outcome}:}{List of objects of class \code{"SaemixOutcome"} giving the outcomes in the model}
#'     \item{\code{noutcome}:}{Number of outcomes}
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
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' @seealso \code{\link{SaemixData}} \code{\link{SaemixObject}} \code{\link{saemixControl}} \code{\link{saemix}}
#' \code{\link{plot.saemix}}

#' @keywords classes
#' @exportClass SaemixModel
#' @examples
#' 
#' showClass("SaemixModel")
#' 

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
  definition=function(.Object, model, parameters, outcome, description="", verbose=TRUE){
#    cat ("--- initialising SaemixModel Object --- \n")
    if(missing(model)) {
      if(verbose) cat("Error initialising SaemixModel object:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
      return("Creation of object SaemixModel failed\n")
    }
    if(missing(parameters)) {
      if(verbose) cat("Error initialising SaemixModel object:\n    Please provide a list of parameters.\n")
      return("Creation of object SaemixModel failed\n")
    }
    if(missing(outcome)) outcome<-list(y=new(Class="SaemixContinuousOutcome")) else {
      if(inherits(outcome,"character")) { # if only name is given assume continuous default outcome
        if(outcome=="") outcome<-"y"
        y1<-vector(mode="list", length=length(outcome))
        names(y1)<-outcome
        for(i in 1:length(outcome)) y1[[i]]<-new(Class="SaemixContinuousOutcome")
        outcome<-y1
      }
      if(inherits(outcome,"SaemixOutcome")) {
        outcome<-list(new(Class="SaemixContinuousOutcome"))
        name(outcome)<-outcome[[1]]@name.outcome
      }
      if(is(outcome,"list")) { # check type
        is.ok<-0
        for(i in 1:length(outcome)) {
          if(!inherits(outcome[[i]],"SaemixOutcome")) is.ok<-1 else {
            if(!is.null(names(outcome))) outcome[[i]]@name.outcome <- names(outcome)[i]
          }
        }
        if(is.ok==1) {
          if(verbose) cat("Error initialising SaemixModel object:\n   Valid outcome values are either a list of SaemixOutcome objects or a vector of names. Please see the documentation for examples.\n")
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
            "description"={return(x@description)},
            "model"={return(x@model)},
            "sim.model"={return(x@sim.model)},
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
            "description"={x@description<-value},
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

####################################################################################
####			SaemixModel class - method to print/show data		####
####################################################################################

#' @rdname print-methods
#' @exportMethod print

setMethod("print","SaemixModel",
          function(x,...) {
            cat("Nonlinear mixed-effects model\n")
            xout <- data.frame(Outcome=names(x@outcome), type= unlist(lapply(x@outcome, function(x) x@type.outcome)), distribution= unlist(lapply(x@outcome, function(x) x@distribution)), unit=unlist(lapply(x@outcome, function(x) x@unit)))
            cat("Outcomes\n")
            print(xout, row.names=FALSE)
            xpar <- data.frame(Parameters=x@param.names, distribution=x@distribution)
            cat("Model parameters\n")
            print(xpar, row.names=FALSE)
            if(length(x@covariate)>0) cat("Covariates",paste(x@covariate,collapse=", "),"\n") else cat("No covariates\n")
            cat("Variability levels:", x@varlevel,"\n")
            if( is.null(body(x@model))) {
              cat("No model function set yet\n")
              return()
            }
            cat("Model function:\n")
            print(x@model)
            if(!is.null(body(x@sim.model))) {
              print(x@sim.model)
            }
          }
)

#' @rdname show-methods
#' @exportMethod show

setMethod("show","SaemixModel",
          function(object) {
            cat("Nonlinear mixed-effects model\n")
            cat("     outcomes:",paste(names(object@outcome),collapse=", "),"\n")
            cat("     parameters:",paste(object@param.names,collapse=", "),"\n")
            if(length(object@covariate)>0) cat("     covariates:",paste(object@covariate,collapse=", "),"\n")
          }
)

#' @rdname showall-methods
#' @exportMethod showall

setMethod("showall","SaemixModel",
          function(object) {
            print(object)
          }
)

## ToDo: Move saemixModel to -methods

####################################################################################
####			SaemixModel class - User-level function			####
####################################################################################

#' Function to create an SaemixModel object
#' 
#' This function creates an SaemixModel object. The two mandatory arguments are
#' the name of a R function computing the model in the SAEMIX format (see
#' details and examples) and a matrix psi0 giving the initial estimates of the
#' fixed parameters in the model, with one row for the population mean
#' parameters and one row for the covariate effects (see documentation).
#' 
#' This function is the user-friendly constructor for the SaemixModel object
#' class.
#' 
#' @param model name of the function used to compute the structural model. The
#' function should return a vector of predicted values given a matrix of
#' individual parameters, a vector of indices specifying which records belong
#' to a given individual, and a matrix of dependent variables (see example
#' below).
#' @param parameters the list of parameters in the model. The preferred format is to use
#' the saemixParam() constructor, as in ka=saemixParam() with arguments describing
#' the distribution of the parameter, constraints and starting values (see \code{\link{saemixParam}} for examples).
#' Alternative forms with less flexibility are to use either a vector of names (eg c("ka","Vd"))
#' which will create parameters with corresponding names, or a named vector of numerical values
#' (eg c(ka=1, Vd=20)) which will create the same parameters and ascribe them starting values.
#' @param outcome the outcome in the model. The preferred format is to use
#' the saemixOutcome() constructors to specify outcomes through their distribution.
#' An alternative form in simple cases is to pass a vector of outcome names (eg c("parent","metabolite")) 
#' which will generate continuous outcomes with the default residual error model (constant variance). 
#' (see \code{\link{saemixOutcome}}, \code{\link{saemixContinuousOutcome}} 
#' and \code{\link{saemixDiscreteOutcome}} for details)
#' @param description a character string, giving a brief description of the
#' model or the analysis
#' @param simulate.function for non-Gaussian data models,
#' the name of the function used to simulate from the structural model. The
#' function should have the same header as the model function, and should return 
#' a vector of simulated values given a matrix of individual parameters, 
#' a vector of indices specifying which records belong to a given individual, 
#' and a matrix of dependent variables (see example in the documentation, section
#' discrete data examples)
#' @param verbose a boolean, controlling whether information about the created should be printed out. Defaults to TRUE
#' @return An SaemixModel object (see \code{\link{saemixModel}}).
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{saemixControl}},\code{\link{saemix}}
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
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' @export saemixModel

saemixModel<-function(model, parameters, outcome="", description="", simulate.function=NULL, verbose=FALSE) {
  # Creating model from class
  if(missing(model)) {
    if(verbose) cat("Error in saemixModel:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
    return("Creation of SaemixModel failed")  
  }
  xcal<-try(typeof(model))
  if(inherits(xcal,"try-error")) {
    if(verbose) cat("Error in saemixModel:\n   the model function does not exist.\n")
    return("Creation of SaemixModel failed")  
  }
  is.valid.model <- TRUE
  if(is(model,"character")) {
    if(exists(model)) model<-get(model) else is.valid.model<-FALSE
  }
  if(!is.function(model)) is.valid.model<-FALSE
  if(!is.valid.model) {
      if(verbose) cat("Error in saemixModel:\n   The argument model to saemixModel must be a valid function.\n")
      return("Creation of SaemixModel failed")
  }
  logmsg<-""
  is.valid.args <- TRUE
  if(length(formals(model))!=3) is.valid.args<-FALSE
  if(!identical(names(formals(model)),c("psi","id","xidep"))) is.valid.args<-FALSE
  if(!is.valid.args) {
    if(verbose) cat("Error in saemixModel:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
    return("Creation of SaemixModel failed")
  }
  has.sim<-FALSE
  if(!is.null(simulate.function)) {
    xcal<-try(typeof(simulate.function))
    if(inherits(xcal,"try-error")) {
        msg <- "The simulate.function does not exist, ignoring.\n"
        logmsg <- paste0(logmsg, msg)
        if(verbose) message(msg)
    } else {
      if(is(simulate.function,"character")) {
        if(exists(simulate.function)) {
          simulate.function<-get(simulate.function)
          if(inherits(simulate.function,"function")) has.sim <- TRUE else {
            msg <- "The simulate.function is not a valid function, ignoring.\n"
            logmsg <- paste0(logmsg, msg)
            if(verbose) message(msg)
            }
        }
      }
      if(!is.function(simulate.function) || length(formals(simulate.function))!=3) {
        msg<-"The simulate.function should have the same format as the model function, ignoring.\n"
        logmsg <- paste0(logmsg, msg)
        if(verbose) message(msg)
        has.sim <- FALSE
      } else has.sim <- TRUE
    }
  }

  if(missing(parameters)) {
    if(verbose) cat("Error in saemixModel:\n   please provide the parameters in the model as a list of saemixParameter objects.\n")
    return("Creation of SaemixModel failed")  
  }
  # test that parameters are valid
  if(inherits(parameters,"character")) { # only parameter names are given as a vector
    msg<-"Parameter given as names only\n"
    logmsg <- paste0(logmsg, msg)
    if(verbose) message(msg)
    xpar<-vector(mode='list',length=length(parameters))
    names(xpar)<-parameters
    for(i in 1:length(parameters)) xpar[[i]]<-saemixParam(name=parameters[i])
    parameters<-xpar
  }
  if(inherits(parameters,"numeric")) { # parameters given as c(ka=1, V=20, etc...)
    msg<-"Parameter given as name=starting value\n"
    logmsg <- paste0(logmsg, msg)
    if(verbose) message(msg)
    xpar<-vector(mode='list',length=length(parameters))
    names(xpar)<-parameters
    for(i in 1:length(parameters)) xpar[[i]]<-saemixParam(name=names(xpar)[i],mu.init=parameters[i])
    parameters<-xpar
  }
  if(inherits(parameters,"list")) { # parameters given as a list
    for(i in 1:length(parameters)) {
      if(inherits(parameters[[i]],"character")) {
        msg<-paste0("Parameter ",parameters[[i]]," given as name only, assuming log-normal distribution\n")
        logmsg <- paste0(logmsg, msg)
        if(verbose) message(msg)
        parameters[[i]]<-saemixParam(name=parameters[[i]]) # given as "CL"
      }
      if(inherits(parameters[[i]],"numeric")) {
        msg<-paste0("Parameter ",parameters[[i]]," given with its starting value, assuming log-normal distribution\n")
        logmsg <- paste0(logmsg, msg)
        if(verbose) message(msg)
        parameters[[i]]<-saemixParam(name=names(parameters)[i],mu.init=parameters[[i]]) # given as CL=0.5
      }
      if(!inherits(parameters[[i]],"SaemixParameter")) {
        if(verbose) cat("Error in saemixModel:\n   please provide the parameters in the model as a list of saemixParameter objects.\n")
        return("Creation of SaemixModel failed")
      }
    }
  }
  xmod <- try(new(Class="SaemixModel", model=model, parameters=parameters, outcome=outcome))
  if(is(xmod,"SaemixModel")) x1<-try(validObject(xmod),silent=FALSE) else x1<-xmod
  if(!inherits(x1,"try-error")) {
    if(verbose) cat("\n\nThe following SaemixModel object was successfully created:\n\n")
  } else xmod<-"Creation of SaemixModel failed"
  if(has.sim) xmod@sim.model <- simulate.function
  xmod@log <- logmsg
  if(verbose) print(xmod)
  return(xmod)
}


