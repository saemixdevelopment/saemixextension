####################################################################################
####			SaemixOutcome class - definition				####
####################################################################################

#' @include aaa_generics.R
NULL

#' Class "SaemixOutcome"
#' 
#' An object of the SaemixOutcome class, representing an observation model for an outcome. Outcomes can be discrete or survival-type data
#' (SaemixDiscreteOutcome) or continuous (SaemixContinuousOutcome), with both subclasses inheriting from the SaemixOutcome class
#' 
#' @name SaemixOutcome-class 
#' @docType class
#' @aliases SaemixOutcome SaemixOutcome-class 
#' @aliases print,SaemixOutcome showall,SaemixOutcome show,SaemixOutcome
#' @aliases SaemixDiscreteOutcome SaemixDiscreteOutcome-class 
#' @aliases print,SaemixDiscreteOutcome showall,SaemixDiscreteOutcome show,SaemixDiscreteOutcome
#' @aliases SaemixContinuousOutcome SaemixContinuousOutcome-class 
#' @aliases print,SaemixContinuousOutcome showall,SaemixContinuousOutcome show,SaemixContinuousOutcome
#' 
#' @section Objects from the Class: 
#' An object of the SaemixOutcome class contains the following slots:
#' @slot name.outcome Object of class \code{"character"}: name given to the outcome
#' @section Methods:
#'   \describe{
#'     \item{[<-}{\code{signature(x = "SaemixData")}: replace elements of object}
#'     \item{[}{\code{signature(x = "SaemixData")}: access elements of object}
#'     \item{initialize}{\code{signature(.Object = "SaemixData")}: internal function to initialise object, not to be used}
#'     \item{print}{\code{signature(x = "SaemixData")}: prints details about the object (more extensive than show)}
#'     \item{showall}{\code{signature(object = "SaemixData")}: shows all the elements in the object}
#'     \item{show}{\code{signature(object = "SaemixData")}: prints details about the object}
#' 	 }
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' 
#' @seealso \code{\link{saemixData}} \code{\link{SaemixModel}} \code{\link{saemixControl}} \code{\link{saemix}}
#' @examples
#' showClass("SaemixOutcome")
#' 
#' @keywords classes
#' @exportClass SaemixOutcome


# Observation class - generic
setClass(Class = "SaemixOutcome",
         representation=representation(
           name.outcome = "character", # outcome name
           type.outcome = "character", # Type: continuous, discrete or event
           distribution = "character" # Distribution 
         ),
         validity=function(object){
           return(TRUE)
         }
)


setMethod( 
  f="initialize",
  signature="SaemixOutcome",
  definition=function(.Object, name.outcome="y", type.outcome="continuous", distribution="normal"){
    .Object@name.outcome <- name.outcome
    .Object@type.outcome <- type.outcome
    .Object@distribution <- distribution
    validObject(.Object)
    return(.Object)
  }
)


# Getteur
setMethod(
  f ="[",
  signature = "SaemixOutcome" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.outcome"={return(x@name.outcome)},
            "type.outcome"={return(x@type.outcome)},
            "distribution"={return(x@distribution)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixOutcome" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.outcome"={x@name.outcome<-value},
            "type.outcome"={x@type.outcome<-value},
            "distribution"={x@distribution<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

########################################################################
# Observation class - Discrete outcome
#' @exportClass SaemixDiscreteOutcome

setClass(
  Class="SaemixDiscreteOutcome",
  contains = "SaemixOutcome",
  validity=function(object){
    # Check type.outcome is one of discrete or event
    if (!(object@type.outcome %in% c("discrete","event"))) {
      message("[ SaemixDiscreteOutcome : validation ] Please specify the type of the outcome (one of discrete, for count or categorical data, or event, for event-type data).")
      return("Outcome type not given")
    }
    return(TRUE)
  }
)

setMethod(
  f="initialize",
  signature="SaemixDiscreteOutcome",
  definition= function (.Object, name.outcome="y", type.outcome="categorical", distribution="binomial") {
    #    cat ("--- initialising SaemixDiscrete Object --- \n")
    .Object@name.outcome<-name.outcome
    .Object@type.outcome<-type.outcome
    .Object@distribution<-distribution
    return (.Object )
  }
)

########################################################################
# Observation class - Continuous outcome
#' @exportClass SaemixContinuousOutcome

setClass(
  Class="SaemixContinuousOutcome",
  contains = "SaemixOutcome",
  representation=representation(
    error.model = "character", # name of the error model associated
    error.npar = "numeric", # number of parameters in the error model
    error.parameters = "numeric", # value of parameters in the error model
    error.nameparameters = "character", # names of parameters in the error model
    error.function = "function"  # error model function
  ),
  validity=function(object){
    # Check error.model is one of constant, proportional, combined1, combined2, power or user
    # Check error.function exists if error.model=user
    #    cat ("--- Checking SaemixData object ---\n")
    return(TRUE)
  }
)

setMethod(
  f="initialize",
  signature="SaemixContinuousOutcome",
  definition= function (.Object, name.outcome="y", type.outcome="continuous", distribution="normal", error.model="constant", error.npar=NULL, error.function=NULL, error.parameters=NULL, error.nameparameters=NULL) {
    #    cat ("--- initialising SaemixDiscrete Object --- \n")
    .Object@name.outcome<-name.outcome
    .Object@type.outcome<-type.outcome
    .Object@distribution<-distribution
    if(!(error.model %in% c("constant","proportional", "combined1","combined2","power", "user"))) {
      message("[ SaemixContinuousOutcome : validation ] Please specify a valid error model for the continuous outcome (one of constant, proportional, combined1, combined2, power or user.")
      return("Error model not given")
    }
    .Object@error.model<-error.model
    if(error.model=="constant") {
      .Object@error.npar<-1
      .Object@error.function<-constantErrorModel
      if(is.null(error.parameters)) .Object@error.parameters<-c(1) else .Object@error.parameters<-error.parameters[1]
      if(is.null(error.nameparameters)) .Object@error.nameparameters<-paste("a",.Object@name.outcome,sep=".") else .Object@error.nameparameters<-error.nameparameters[1]
    }
    if(error.model=="proportional") {
      .Object@error.npar<-1
      .Object@error.function<-proportionalErrorModel
      if(is.null(error.parameters)) .Object@error.parameters<-c(1) else .Object@error.parameters<-error.parameters[1]
      if(is.null(error.nameparameters)) .Object@error.nameparameters<-paste("b",.Object@name.outcome,sep=".") else .Object@error.nameparameters<-error.nameparameters[1]
    }
    if(error.model=="combined1") {
      .Object@error.npar<-2
      .Object@error.function<-combined1ErrorModel
      if(is.null(error.parameters)) .Object@error.parameters<-c(1,1) else .Object@error.parameters<-error.parameters[1:2]
      if(is.null(error.nameparameters)) .Object@error.nameparameters<-paste(c("a","b"),.Object@name.outcome,sep=".") else .Object@error.nameparameters<-error.nameparameters[2]
    }
    if(error.model=="combined2") {
      .Object@error.npar<-2
      .Object@error.function<-combined2ErrorModel
      if(is.null(error.parameters)) .Object@error.parameters<-c(1,1) else .Object@error.parameters<-error.parameters[1:2]
      if(is.null(error.nameparameters)) .Object@error.nameparameters<-paste(c("a","b"),.Object@name.outcome,sep=".") else .Object@error.nameparameters<-error.nameparameters[2]
    }
    if(error.model=="power") {
      .Object@error.npar<-3
      .Object@error.function<-powerErrorModel
      if(is.null(error.parameters)) .Object@error.parameters<-c(1,1,2) else .Object@error.parameters<-error.parameters[1:3]
      if(is.null(error.nameparameters)) .Object@error.nameparameters<-paste(c("a","b","c"),.Object@name.outcome,sep=".") else .Object@error.nameparameters<-error.nameparameters[3]
    }
    if(error.model=="user") {
      if(is.null(error.npar) | is.null(error.function) | !(is.function(error.function))) {
        message("[ SaemixContinuousOutcome : validation ] When specifying a user-defined model, please give the number of parameters and a valid function with 2 arguments, f and ab, where the number of parameters is the size of the vector ab.")
        return("Error function or number of parameters not given")
      }
      if(!identical(names(formals(error.function)),c("f","ab"))) {
        message("[ SaemixContinuousOutcome : validation ] The error model should be a function with 2 arguments, f and ab, where the number of parameters is the size of the vector ab.")
        return("Error function arguments mismatch")
      }
      .Object@error.npar<-error.npar
      .Object@error.function<-error.function
      if(is.null(error.nameparameters)) error.nameparameters<-letters[1:error.npar]
      .Object@error.nameparameters<-paste(error.nameparameters,.Object@name.outcome,sep=".")
      if(is.null(error.parameters)) .Object@error.parameters<-rep(1,error.npar) else .Object@error.parameters<-error.parameters[1:error.npar]
    }
    return (.Object )
  }
)


# power model

user.error1<-function(f,ab) {
  g<-cutoff(sqrt((ab[1]+ab[2]*f)^ab[3]))
  return(g)
}

constantErrorModel<-function(f,ab) {
  g<-cutoff(ab[1])
  return(g)
}

proportionalErrorModel<-function(f,ab) {
  g<-cutoff(ab[1]*abs(f))
  return(g)
}

combined1ErrorModel<-function(f,ab) {
  g<-cutoff(ab[1]+ab[2]*abs(f))
  return(g)
}

# Johannes 02/21
combined2ErrorModel<-function(f,ab) {
  g<-cutoff(sqrt(ab[1]^2+ab[2]^2*f^2))  
  return(g)
}

powerErrorModel<-function(f,ab) {
  g<-cutoff(sqrt(ab[1]^2+ab[2]^2*f^ab[3]))
  return(g)
}


# Getteur
# setMethod(
#   f ="[",
#   signature = "SaemixDiscreteOutcome" ,
#   definition = function (x,i,j,drop ){
#     switch (EXPR=i,
#             "name.outcome"={return(x@name.outcome)},
#             "type.outcome"={return(x@type.outcome)},
#             "distribution"={return(x@distribution)},
#             stop("No such attribute\n")
#     )
#   }
# )

setMethod(
  f ="[",
  signature = "SaemixContinuousOutcome" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            # "name.outcome"={return(x@name.outcome)},
            # "type.outcome"={return(x@type.outcome)},
            # "distribution"={return(x@distribution)},
            "error.model"={return(x@error.model)},
            "error.npar"={return(x@error.npar)},
            "error.nameparameters"={return(x@error.nameparameters)},
            "error.parameters"={return(x@error.parameters)},
            "error.function"={return(x@error.function)},
            stop("No such attribute\n")
    )
  }
)


# Setteur
# setReplaceMethod(
#   f ="[",
#   signature = "SaemixDiscreteOutcome" ,
#   definition = function (x,i,j,value){
#     switch (EXPR=i,
#             "name.outcome"={x@name.outcome<-value},
#             "type.outcome"={x@type.outcome<-value},
#             "distribution"={x@distribution<-value},
#             stop("No such attribute\n")
#     )
#     validObject(x)
#     return(x)
#   }
# )


setReplaceMethod(
  f ="[",
  signature = "SaemixContinuousOutcome" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            # "name.outcome"={x@name.outcome<-value},
            # "type.outcome"={x@type.outcome<-value},
            # "distribution"={x@distribution<-value},
            "error.model"={x@error.model<-value},
            "error.npar"={x@error.npar<-value},
            "error.nameparameters"={x@error.nameparameters<-value},
            "error.parameters"={x@error.parameters<-value},
            "error.function"={x@error.function<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

########################################################################
# Show

setMethod("show","SaemixContinuousOutcome",
          function(object) {
            cat("Continuous outcome:",object@name.outcome,"\n")
            cat("   ",object@distribution,"distribution with ")
            if(object@error.model!="user") 
              cat(paste0(object@error.model," residual error model involving ",object@error.npar," parameter",ifelse(object@error.npar>1,"s",""),"\n")) else {
                cat(paste0("user-defined error model involving ",object@error.npar," parameter",ifelse(object@error.npar>1,"s",""),":\n"))
                print(object@error.function)
              }
          }
)

setMethod("showall","SaemixContinuousOutcome",
          function(object) {
            cat("Continuous outcome:",object@name.outcome,"\n")
            cat("   ",object@distribution,"distribution with ")
            if(object@error.model!="user") 
              cat(paste0(object@error.model," residual error model involving ",object@error.npar," parameter",ifelse(object@error.npar>1,"s",":")),"\n") else
                cat(paste0("user-defined error model involving ",object@error.npar," parameter",ifelse(object@error.npar>1,"s",""),":\n"))
            cat(object@error.nameparameters,"\n")
            print(object@error.function)
            cat("Initial parameter values:", paste(object@error.nameparameters,object@error.parameters,sep="="),"\n")
          }
)

setMethod("print","SaemixContinuousOutcome",
          function(x,nlines=10,...) {
            show(x)
          }
)

setMethod("show","SaemixDiscreteOutcome",
          function(object) {
            cat("Discrete outcome:",object@name.outcome,"\n")
            cat("    type:", object@type.outcome,"\n")
#            cat("    distribution:", object@distribution,"\n")
          }
)
setMethod("showall","SaemixDiscreteOutcome",
          function(object) {
            cat("Discrete outcome:",object@name.outcome,"\n")
            cat("    type:", object@type.outcome,"\n")
            cat("    distribution:", object@distribution,"\n")
          }
)

setMethod("print","SaemixDiscreteOutcome",
          function(x,nlines=10,...) {
            show(x)
          }
)
