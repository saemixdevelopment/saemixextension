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
           log = "character", # warning messages
           name.outcome = "character", # outcome name
           type.outcome = "character", # Type: continuous, discrete or event
           unit = "character", # unity of the outcome
           distribution = "character", # Distribution name
           density = "function", # density function for the outcome
           density.param = "numeric", # named vector of parameters for the density distribution
           simulate.function = "function" # function used to simulate outcome
         ),
         validity=function(object){
           return(TRUE)
         }
)


setMethod( 
  f="initialize",
  signature="SaemixOutcome",
  definition=function(.Object, name.outcome="y", type.outcome="continuous", distribution="normal", density=dnorm, density.param=c(mean=0, sd=1), simulate.function=rnorm, unit="", verbose=FALSE){
    .Object@log <- ""
    .Object@name.outcome <- name.outcome
    .Object@type.outcome <- type.outcome
    .Object@unit <- unit
    .Object@distribution <- tolower(distribution)
    .Object@density <- density
    .Object@density.param <- density.param
    .Object@simulate.function <- simulate.function
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
            "log"={return(x@log)},
            "name.outcome"={return(x@name.outcome)},
            "type.outcome"={return(x@type.outcome)},
            "unit"={return(x@unit)},
            "distribution"={return(x@distribution)},
            "density"={return(x@density)},
            "density.param"={return(x@density.param)},
            "simulate.function"={return(x@simulate.function)},
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
            "log"={x@log<-value},
            "name.outcome"={x@name.outcome<-value},
            "type.outcome"={x@type.outcome<-value},
            "unit"={x@unit<-value},
            "distribution"={x@distribution<-value},
            "density"={x@density<-value},
            "density.param"={x@density.param<-value},
            "simulate.function"={x@simulate.function<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

########################################################################
# Observation class - Discrete outcome (binary, categorical, count)
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
  definition= function (.Object, name.outcome="y", type.outcome="categorical", distribution="binomial", density=dbinom, density.param=c(size=1, prob=0.5), simulate.function=rbinom, unit="", verbose=FALSE) {
    #    cat ("--- initialising SaemixDiscrete Object --- \n")
    .Object@name.outcome<-name.outcome
    .Object@type.outcome<-type.outcome
    distribution<-tolower(distribution)
    .Object@distribution<-distribution
    .Object@unit<-unit
    logmsg<-""
    if(type.outcome=="categorical" & distribution %in% c("binomial","Poisson")) {
      msg<-paste("Setting density and simulation function to",distribution," distribution\n")
      if(verbose) message(msg)
      logmsg<-paste0(logmsg,msg)
    }
    if(type.outcome=="event" & distribution %in% c("exponential","weibull","logistic")) {
      msg<-paste("Setting density and simulation function to match the",distribution,"distribution\n")
      if(verbose) message(msg)
      logmsg<-paste0(logmsg,msg)
    }
    if(type.outcome=="categorical" & distribution %in% c("exponential","weibull","logistic")) {
      msg<-paste("Warning: ",distribution,"distribution is for event variables but type of outcome is set at categorical\n")
      if(verbose) message(msg)
      logmsg<-paste0(logmsg,msg)
    }
    if(type.outcome=="event" & distribution %in% c("binomial","poisson")) {
      msg<-paste("Warning: ",distribution,"distribution is for categorical variables but type of outcome is set at event\n")
      if(verbose) message(msg)
      logmsg<-paste0(logmsg,msg)
    }
    if(distribution=="poisson" & type.outcome=="categorical") {
      density <- dpois
      density.param <- c(lambda=1)
      simulate.function <- rpois
    }
    if(distribution=="binomial" & type.outcome=="categorical") {
      density <- dbinom
      density.param <- c(size=1, prob=0.5)
      simulate.function <- rbinom
    }
    if(distribution=="exponential" & type.outcome=="event") {
      density <- dexp
      density.param <- c(rate=1)
      simulate.function <- rexp
    }
    if(distribution=="weibull" & type.outcome=="event") {
      density <- dweibull
      density.param <- c(shape=1, scale=1)
      simulate.function <- rweibull
    }
    if(distribution=="logistic" & type.outcome=="event") {
      density <- dlogis
      density.param <- c(location=0, scale=1)
      simulate.function <- rlogis
    }
    .Object@density <- density
    .Object@density.param <- density.param
    .Object@simulate.function <- simulate.function
    .Object@log <- logmsg
    return (.Object )
  }
)

# Getteur
setMethod(
  f ="[",
  signature = "SaemixDiscreteOutcome" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "log"={return(x@log)},
            "name.outcome"={return(x@name.outcome)},
            "type.outcome"={return(x@type.outcome)},
            "unit"={return(x@unit)},
            "distribution"={return(x@distribution)},
            "density"={return(x@density)},
            "density.param"={return(x@density.param)},
            "simulate.function"={return(x@simulate.function)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixDiscreteOutcome" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "log"={x@log<-value},
            "name.outcome"={x@name.outcome<-value},
            "type.outcome"={x@type.outcome<-value},
            "unit"={x@unit<-value},
            "distribution"={x@distribution<-value},
            "density"={x@density<-value},
            "density.param"={x@density.param<-value},
            "simulate.function"={x@simulate.function<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

########################################################################
# Observation class - Continuous outcome
#' @exportClass SaemixContinuousOutcome

setClass(
  Class="SaemixContinuousOutcome",
  contains = "SaemixOutcome",
  representation=representation(
    error.model = "SaemixErrorModel" # an SaemixErrorModel
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
  definition= function (.Object, name.outcome="y", type.outcome="continuous", distribution="normal", density=dnorm, density.param=c(mean=0, sd=1), simulate.function=rnorm, error.model="constant", unit="", verbose=FALSE) {
    #    cat ("--- initialising SaemixDiscrete Object --- \n")
    logmsg<-""
    .Object@name.outcome<-name.outcome
    .Object@type.outcome<-type.outcome
    .Object@distribution<-distribution
    .Object@density <- density
    .Object@density.param <- density.param # distribution of epsilon, SD for Y is actually given by the error model
    .Object@simulate.function <- simulate.function 
    .Object@unit <- unit
    if(is(error.model,"character")) { # only name of error model is given
      if(!(error.model %in% c("constant","proportional", "combined1","combined2","power"))) {
        if(verbose) message("[ SaemixContinuousOutcome : validation ] Please specify a valid error model name for the continuous outcome (one of constant, proportional, combined1, combined2, power) or an SaemixErrorModel object.")
        return("Error model not recognised")
      }
      if(error.model=="user") {
        if(verbose) message("[ SaemixContinuousOutcome : validation ] For user-defined error models, please provide an SaemixErroModel object.")
        return("Error model incomplete")
      }
      error.model <- new(Class="SaemixErrorModel", name=error.model)
    }
    if(!is(error.model,"SaemixErrorModel")) {
        if(verbose) message("[ SaemixContinuousOutcome : validation ] An error model must be provided either through its name or through an SaemixErrorModel.")
        return("Invalid error model")
      }
    error.model@param.names <- paste0(error.model@param.names,".",.Object@name.outcome)
    .Object@error.model<-error.model
    return (.Object )
  }
)



setMethod(
  f ="[",
  signature = "SaemixContinuousOutcome" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "log"={return(x@log)},
            "name.outcome"={return(x@name.outcome)},
            "type.outcome"={return(x@type.outcome)},
            "unit"={return(x@unit)},
            "distribution"={return(x@distribution)},
            "density"={return(x@density)},
            "density.param"={return(x@density.param)},
            "simulate.function"={return(x@simulate.function)},
            "error.model"={return(x@error.model)},
            stop("No such attribute\n")
    )
  }
)

setReplaceMethod(
  f ="[",
  signature = "SaemixContinuousOutcome" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "log"={x@log<-value},
            "name.outcome"={x@name.outcome<-value},
            "type.outcome"={x@type.outcome<-value},
            "unit"={x@unit<-value},
            "distribution"={x@distribution<-value},
            "density"={x@density<-value},
            "density.param"={x@density.param<-value},
            "simulate.function"={x@simulate.function<-value},
            "error.model"={x@error.model<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

########################################################################
# Show/print

## Continuous Outcome
setMethod("print","SaemixContinuousOutcome",
          function(x,nlines=10,...) {
            show(x)
          }
)

setMethod("show","SaemixContinuousOutcome",
          function(object) {
            cat("Continuous outcome:",object@name.outcome,",",object@distribution,"distribution with",object@error.model@name,"error\n")
          }
)

setMethod("showall","SaemixContinuousOutcome",
          function(object) {
            cat("Continuous outcome:",object@name.outcome,"\n")
            cat("   ",object@distribution,"distribution with",object@error.model@name,"error\n")
            showall(object@error.model)
          }
)

## Discrete Outcome
setMethod("print","SaemixDiscreteOutcome",
          function(x,nlines=10,...) {
            show(x)
          }
)

setMethod("show","SaemixDiscreteOutcome",
          function(object) {
            cat("Discrete outcome:",object@name.outcome," of type", object@type.outcome,"with",object@distribution,"distribution\n")
          }
)
setMethod("showall","SaemixDiscreteOutcome",
          function(object) {
            cat("Discrete outcome:",object@name.outcome,"\n")
            cat("    type:", object@type.outcome,"\n")
            cat("    distribution:", object@distribution,"\n")
            cat("    density, with parameters",object@density.param,":")
            print(object@density)
            cat("    function to simulate from the distribution:")
            print(object@simulate.function)
          }
)

########################################################################
# Creator function

saemixOutcome <- function(name.outcome="y", type.outcome="continuous", distribution="normal", unit="",
                    density=dnorm, density.param=c(mean=0, sd=1), simulate.function=rnorm, 
                    error.model="constant", verbose=FALSE) {
  if(!(type.outcome %in% c("continuous","binary","categorical","count","event"))) {
    return("Creation of outcome failed: the type.outcome must be of one of continuous, categorical, count, event\n")
  }
  logmsg<-""
  if(type.outcome=="continuous") {
    if(distribution!="normal") {
      msg<-paste("Warning: ",distribution,"distribution not valid for continuous outcome, setting it to normal\n")
      if(verbose) message(msg)
      logmsg<-paste0(logmsg,msg)
      distribution<-"normal"
    }
    density <- dnorm
    density.param <- c(mean=0, sd=1)
    simulate.function <- rnorm
    y<- new(Class="SaemixContinuousOutcome", name.outcome=name.outcome, type.outcome="continuous", distribution=distribution, 
            unit=unit, density=density, density.param=density.param, simulate.function=simulate.function, 
            error.model=error.model, verbose=verbose)
  }
  
  if(type.outcome=="binary")  {
    if(distribution!="binomial") {
      msg<-paste("With binary outcome, the distribution should be binomial. Setting density and simulation function to match outcome type\n")
      if(verbose) message(msg)
      logmsg<-paste0(logmsg,msg)
      distribution<-"binomial"
    }
    density <- dbinom
    density.param <- c(size=1, prob=0.5)
    simulate.function <- rbinom
    y<- new(Class="SaemixDiscreteOutcome", name.outcome=name.outcome, type.outcome="categorical", unit=unit,
            distribution=distribution, density=density, density.param=density.param, simulate.function=simulate.function, verbose=verbose)
  }
  if(type.outcome %in% c("categorical","count","event")) {
    if(distribution=="normal" & type.outcome=="categorical") {
      distribution <- "binomial"
      msg<-paste("Distribution for type",type.outcome,"not given, assuming the outcome has a binomial distribution\n")
      if(verbose) message(msg)
      logmsg<-paste0(logmsg,msg)
    }
    if(distribution=="normal" & type.outcome=="event") {
      distribution <- "exponential"
      msg<-paste("Distribution for type",type.outcome,"not given, assuming an exponential hazard model\n")
      if(verbose) message(msg)
      logmsg<-paste0(logmsg,msg)
    }
    if(distribution=="normal" & type.outcome=="count") {
      distribution <- "poisson"
      msg<-paste("Distribution for type",type.outcome,"not given, assuming the outcome has a binomial distribution\n")
      if(verbose) message(msg)
      logmsg<-paste0(logmsg,msg)
    }
    if(distribution=="binomial") {
      density <- dbinom
      density.param <- c(size=1, prob=0.5)
      simulate.function <- rbinom
    }
    if(distribution=="exponential") {
      density <- dexp
      density.param <- c(rate=1)
      simulate.function <- rexp
    }
    if(distribution=="poisson") {
      density <- dpois
      density.param <- c(lambda=1)
      simulate.function <- rpois
    }
    y<- new(Class="SaemixDiscreteOutcome", name.outcome=name.outcome, type.outcome=type.outcome, 
            distribution=distribution, density=density, density.param=density.param, simulate.function=simulate.function, verbose=verbose)
  }
  y@log<-paste0(y@log, logmsg)
  return(y)
}


  

