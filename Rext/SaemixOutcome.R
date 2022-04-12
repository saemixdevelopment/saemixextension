####################################################################################
####			SaemixOutcome class - definition				####
####################################################################################

#' @include aaa_generics.R
NULL

#' Class "SaemixOutcome"
#' 
#' An object of the SaemixOutcome class, representing an observation model for an outcome. 
#' Outcomes can be discrete (SaemixDiscreteOutcome), survival-type (SaemixEventOutcome) 
#' or continuous (SaemixContinuousOutcome).
#' Both subclasses SaemixDiscreteOutcome and SaemixEventOutcome inherit 
#' from the SaemixOtherOutcome subclass,
#' and the SaemixContinuousOutcome and SaemixOtherOutcome subclasses inherit from
#' the SaemixOutcome class.
#' 
#' @name SaemixOutcome-class 
#' @docType class
#' @aliases SaemixOutcome SaemixOutcome-class 
#' @aliases print,SaemixOutcome showall,SaemixOutcome show,SaemixOutcome
#' @aliases SaemixOtherOutcome SaemixOtherOutcome-class 
#' @aliases SaemixDiscreteOutcome SaemixDiscreteOutcome-class 
#' @aliases SaemixEventOutcome SaemixEventOutcome-class 
#' @aliases print,SaemixDiscreteOutcome showall,SaemixDiscreteOutcome show,SaemixDiscreteOutcome
#' @aliases print,SaemixEventOutcome showall,SaemixEventOutcome show,SaemixEventOutcome
#' @aliases SaemixContinuousOutcome SaemixContinuousOutcome-class 
#' @aliases print,SaemixContinuousOutcome showall,SaemixContinuousOutcome show,SaemixContinuousOutcome
#' 
#' @section Objects from the Class: 
#' An object of the SaemixOutcome class contains the following slots:
#' @slot name Object of class \code{"character"}: name given to the outcome
#' @section Methods:
#'   \describe{
#'     \item{[<-}{\code{signature(x = "SaemixOutcome")}: replace elements of object}
#'     \item{[}{\code{signature(x = "SaemixOutcome")}: access elements of object}
#'     \item{initialize}{\code{signature(.Object = "SaemixOutcome")}: internal function to initialise object, not to be used}
#'     \item{print}{\code{signature(x = "SaemixOutcome")}: prints details about the object (more extensive than show)}
#'     \item{showall}{\code{signature(object = "SaemixOutcome")}: shows all the elements in the object}
#'     \item{show}{\code{signature(object = "SaemixOutcome")}: prints details about the object}
#' 	 }
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' 
#' @seealso \code{\link{saemixData}} \code{\link{SaemixModel}} \code{\link{saemixControl}} \code{\link{saemix}}
#' @examples
#' showClass("SaemixOutcome")
#' concentration<-continuousOutcome()
#' 
#' @keywords classes
#' @exportClass SaemixOutcome


# Observation class - generic
setClass(Class = "SaemixOutcome",
         representation=representation(
           name = "character", # outcome name
           type = "character", # Type: continuous, discrete or event
           unit = "character" # Unit (defaults to empty="")
         ),
         validity=function(object){
           return(TRUE)
         }
)


setMethod( 
  f="initialize",
  signature="SaemixOutcome",
  definition=function(.Object, name="y", type="continuous", distribution="normal", unit=""){
    .Object@name <- name
    .Object@type <- type
    .Object@unit <- unit
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
            "name"={return(x@name)},
            "type"={return(x@type)},
            "unit"={return(x@unit)},
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
            "name"={x@name<-value},
            "type"={x@type<-value},
            "unit"={x@unit<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

########################################################################
# Observation class - Discrete outcome
setClass(
  Class="SaemixOtherOutcome",
  contains = "SaemixOutcome",
  validity=function(object){
    # Check type is discrete
    if (!(object@type %in% c("binary","categorical","count","event"))) {
      message("[ SaemixOtherOutcome : validation ] Please specify the type of the outcome (binary, count, categorical, for discrete data, or event, for event-type data).")
      return("Outcome type not given")
    }
    return(TRUE)
  }
)

# Check if I need to initialise or if I can skip directly to initialising the children
# setMethod(
#   f="initialize",
#   signature="SaemixDiscreteOutcome",
#   definition= function (.Object, name="y", type="categorical", distribution="binomial") {
#     #    cat ("--- initialising SaemixDiscrete Object --- \n")
#     .Object@name<-name
#     .Object@type<-type
#     .Object@distribution<-distribution
#     return (.Object )
#   }
# )


########################################################################
# Observation class - Discrete outcome
#' @exportClass SaemixDiscreteOutcome

setClass(
  Class="SaemixDiscreteOutcome",
  contains = "SaemixOtherOutcome",
  representation=representation(
    distribution = "character", # Distribution 
    levels = "factor" # categories
  ),  
  validity=function(object){
    # Check type is discrete
    if (!(object@type %in% c("binary","categorical","count"))) {
      message("[ SaemixDiscreteOutcome : validation ] Please specify the type of the outcome (binary, count or categorical).")
      return("Outcome type not given")
    }
    if(object@type %in% c("categorical","count") & length(object@levels)<3) {
      message("[ SaemixDiscreteOutcome : validation ] Outcome is of type categorical or count but levels have not been set")
    }
    return(TRUE)
  }
)

setMethod(
  f="initialize",
  signature="SaemixDiscreteOutcome",
  definition= function (.Object, name="y", type, distribution="", levels=NULL) {
    #    cat ("--- initialising SaemixDiscrete Object --- \n")
    .Object@name<-name
    if(missing(type) & distribution=="") {
      if(is.null(levels)) {
        message("[ SaemixDiscreteOutcome : validation ] Please specify the type of the outcome (binary, count or categorical) or the distribution of the outcome (one of bernouilli (binary data), categorical (categorical data), poisson, binomial, negativebinomial, poissonzip (count data), or user (please also specify type in this case)).")
        return("Missing distribution and/or type")
      } # else use nb of levels to toggle between binary or categorical
      if(length(levels)==2) type<-"binary"
      if(length(levels)>2) {
        message("Type of outcome and distribution not given, inferring a categorical distribution")
        type<-"categorical"
      }
    }
    if(distribution=="") {
      if(missing(type)) {
        message("[ SaemixDiscreteOutcome : validation ] Please specify the type of the outcome (binary, count or categorical) or the distribution of the outcome (one of bernouilli (binary data), categorical (categorical data), poisson, binomial, negativebinomial, poissonzip (count data), or user (please also specify type in this case)).")
        return("Missing distribution and/or type")
      }
      if(type=="binary") distribution<-"bernouilli"
      if(type=="categorical") distribution<-"categorical"
      if(type=="count") distribution<-"poisson"
    }
    if(distribution=="categorical") type<-"categorical"
    if(distribution=="bernouilli") type<-"binary"
    if(distribution %in% c("poisson", "binomial", "negativebinomial", "poissonzip")) type<-"count"
    .Object@type<-type
    .Object@distribution<-distribution
    if(is.null(levels)) levels<-c(0,1)
    if(!is.factor(levels)) levels<-factor(levels, levels=levels)
    .Object@levels<-levels
    return (.Object )
  }
)

# Getteur
setMethod(
  f ="[",
  signature = "SaemixDiscreteOutcome" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "levels"={return(x@levels)},
            "distribution"={return(x@distribution)},
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
            "levels"={x@levels<-value},
            "distribution"={x@distribution<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

########################################################################
# Observation class - Survival outcome
#' @exportClass SaemixEventOutcome

setClass(
  Class="SaemixEventOutcome",
  contains = "SaemixOtherOutcome",
  representation=representation(
    distribution = "character", # Distribution 
    maxEvents = "numeric", # maximum number of events (1 for single TTE)
    intervalCensored = "logical", # TRUE if observation of event is interval censored, defaults to FALSE
    intervalLength = "numeric"  # if intervalCensored is TRUE, length of observation interval
  ),  
  validity=function(object){
    # Check type is one of discrete or event
    if (!(object@type %in% c("event"))) {
      message("[SaemixEventOutcome : validation ] Please specify the type of the outcome (event, for event-type data).")
      return("Outcome type not given")
    }
    if (!(is.infinite(object@maxEvents))) {
      if(object@maxEvents<=0) {
        message("[SaemixEventOutcome : validation ] Maximum number of events should be a positive integer.")
        return("maxEvents negative")
      }
    }
    if (object@intervalCensored && intervalLength<=0) {
      message("[SaemixEventOutcome : validation ] Please give the length of the observation interval")
      return("Interval censored event but length of observation interval missing")
    }
    return(TRUE)
  }
)

setMethod(
  f="initialize",
  signature="SaemixEventOutcome",
  definition= function (.Object, name="y", distribution="constantHazard", type="event", maxEvents=1, intervalCensored=FALSE, intervalLength=0) {
    #    cat ("--- initialising SaemixEventOutcome Object --- \n")
    .Object@name<-name
    .Object@type<-"event"
    .Object@distribution<-distribution
    .Object@maxEvents<-maxEvents
    .Object@intervalCensored<-intervalCensored
    .Object@intervalLength<-intervalLength
    return (.Object )
  }
)

# Getteur
setMethod(
  f ="[",
  signature = "SaemixEventOutcome" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "maxEvents"={return(x@maxEvents)},
            "intervalCensored"={return(x@intervalCensored)},
            "intervalLength"={return(x@intervalLength)},
            "distribution"={return(x@distribution)},
            stop("No such attribute\n")
    )
  }
)
# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixEventOutcome" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "maxEvents"={x@maxEvents<-value},
            "intervalCensored"={x@intervalCensored<-value},
            "intervalLength"={x@intervalLength<-value},
            "distribution"={x@distribution<-value},
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
    distribution = "character", # Distribution 
    error.model = "character", # name of the error model associated
    error.npar = "numeric", # number of parameters in the error model
    error.parameters = "numeric", # value of parameters in the error model
    error.nameparameters = "character", # names of parameters in the error model
    error.function = "function",  # error model function
    error.fix = "logical" # fixed parameters (not estimated)
  ),
  validity=function(object){
    # Check error.model is one of constant, proportional, combined1, combined2, power, exponential or user
    # Check error.function exists if error.model=user
    # Check distribution is normal (message otherwise)
    #    cat ("--- Checking SaemixData object ---\n")
    return(TRUE)
  }
)

setMethod(
  f="initialize",
  signature="SaemixContinuousOutcome",
  definition= function (.Object, name="y", type="continuous", distribution="normal", error.model="constant", error.npar=NULL, error.function=NULL, error.parameters=NULL, error.nameparameters=NULL, error.fix=NULL) {
    #    cat ("--- initialising SaemixDiscrete Object --- \n")
    .Object@name<-name
    .Object@type<-type
    .Object@distribution<-distribution
    if(!(error.model %in% c("constant","proportional", "combined1","combined2","power", "exponential", "user"))) {
      message("[ SaemixContinuousOutcome : validation ] Please specify a valid error model for the continuous outcome (one of constant, proportional, combined1, combined2, power or user.")
      return("Error model not given")
    }
    .Object@error.model<-error.model
    if(error.model=="constant" | error.model=="exponential") {
      .Object@error.npar<-1
      .Object@error.function<-constantErrorModel
      if(is.null(error.parameters)) .Object@error.parameters<-c(1) else .Object@error.parameters<-error.parameters[1]
      if(is.null(error.nameparameters)) .Object@error.nameparameters<-paste("a",.Object@name,sep=".") else .Object@error.nameparameters<-error.nameparameters[1]
    }
    if(error.model=="proportional") {
      .Object@error.npar<-1
      .Object@error.function<-proportionalErrorModel
      if(is.null(error.parameters)) .Object@error.parameters<-c(1) else .Object@error.parameters<-error.parameters[1]
      if(is.null(error.nameparameters)) .Object@error.nameparameters<-paste("b",.Object@name,sep=".") else .Object@error.nameparameters<-error.nameparameters[1]
    }
    if(error.model=="combined1") {
      .Object@error.npar<-2
      .Object@error.function<-combined1ErrorModel
      if(is.null(error.parameters)) .Object@error.parameters<-c(1,1) else .Object@error.parameters<-error.parameters[1:2]
      if(is.null(error.nameparameters)) .Object@error.nameparameters<-paste(c("a","b"),.Object@name,sep=".") else .Object@error.nameparameters<-error.nameparameters[2]
    }
    if(error.model=="combined2") {
      .Object@error.npar<-2
      .Object@error.function<-combined2ErrorModel
      if(is.null(error.parameters)) .Object@error.parameters<-c(1,1) else .Object@error.parameters<-error.parameters[1:2]
      if(is.null(error.nameparameters)) .Object@error.nameparameters<-paste(c("a","b"),.Object@name,sep=".") else .Object@error.nameparameters<-error.nameparameters[2]
    }
    if(error.model=="power") {
      .Object@error.npar<-3
      .Object@error.function<-powerErrorModel
      if(is.null(error.parameters)) .Object@error.parameters<-c(1,1,2) else .Object@error.parameters<-error.parameters[1:3]
      if(is.null(error.nameparameters)) .Object@error.nameparameters<-paste(c("a","b","c"),.Object@name,sep=".") else .Object@error.nameparameters<-error.nameparameters[3]
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
      .Object@error.nameparameters<-paste(error.nameparameters,.Object@name,sep=".")
      if(is.null(error.parameters)) .Object@error.parameters<-rep(1,error.npar) else .Object@error.parameters<-error.parameters[1:error.npar]
    }
    if(length(error.fix)==0) error.fix<-rep(FALSE,.Object@error.npar) else {
      if(length(error.fix)!=.Object@error.npar) error.fix<-rep(error.fix, length.out=.Object@error.npar)
    }
    .Object@error.fix<-error.fix
    return (.Object )
  }
)

############ Create outcomes from a list

createSaemixOutcome<-function(response) {
  if(response$type=="continuous") saemix.outcome<-new(Class="SaemixContinuousOutcome", error.model=response$error.model, error.npar=outcome$error.npar, error.function=response$error.function, error.parameters=response$start, error.fix=response$error.fix)
  if(response$type %in% c("binary","categorical")) saemix.outcome<-new(Class="SaemixDiscreteOutcome", levels=response$levels, distribution=response$distribution, type=response$type)
  if(response$type=="event") saemix.outcome<-new(Class="SaemixEventOutcome", distribution=response$distribution, maxEvents=response$maxEvents, intervalCensored=response$intervalCensored, intervalLength=response$intervalLength)
  return(saemix.outcome)
}

################################################################################

# Set/Get Methods
setMethod(
  f ="[",
  signature = "SaemixContinuousOutcome" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            # "name"={return(x@name)},
            # "type"={return(x@type)},
            "distribution"={return(x@distribution)},
            "error.model"={return(x@error.model)},
            "error.npar"={return(x@error.npar)},
            "error.nameparameters"={return(x@error.nameparameters)},
            "error.parameters"={return(x@error.parameters)},
            "error.fix"={return(x@error.fix)},
            "error.function"={return(x@error.function)},
            stop("No such attribute\n")
    )
  }
)

setReplaceMethod(
  f ="[",
  signature = "SaemixContinuousOutcome" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            # "name"={x@name<-value},
            # "type"={x@type<-value},
            "distribution"={x@distribution<-value},
            "error.model"={x@error.model<-value},
            "error.npar"={x@error.npar<-value},
            "error.nameparameters"={x@error.nameparameters<-value},
            "error.parameters"={x@error.parameters<-value},
            "error.fix"={x@error.fix<-value},
            "error.function"={x@error.function<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

################################################################################
# Error models

#' @rdname saemix.internal
#' 
#' @aliases cutoff cutoff.eps cutoff.max cutoff.res
#' 
#' @keywords internal

cutoff<-function(x,seuil=.Machine$double.xmin) {x[x<seuil]<-seuil; return(x)}
cutoff.max<-function(x) max(x,.Machine$double.xmin)
cutoff.eps<-function(x) max(x,.Machine$double.eps)
cutoff.res<-function(x,ares,bres) max(ares+bres*abs(x),.Machine$double.xmin)

user.error1<-function(f,ab) {
  g<-cutoff(sqrt((ab[1]+ab[2]*f)^ab[3]))
  return(g)
}

constantErrorModel<-function(f,ab) {
  g<-rep(cutoff(ab[1]), length(f))
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


################################################################################
# Show and print methods

# Continuous outcome
setMethod("show","SaemixContinuousOutcome",
          function(object) {
            cat("Continuous outcome:",object@name,"\n")
            cat("   ",object@distribution,"distribution with ")
            if(object@error.model!="user") 
              cat(paste0(object@error.model," residual error model involving ",object@error.npar," parameter",ifelse(object@error.npar>1,"s",""),"\n")) else {
                cat(paste0("user-defined error model involving ",object@error.npar," parameter",ifelse(object@error.npar>1,"s",""),":\n"))
                print(object@error.function)
              }
            if(sum(object@error.fix)>0) cat("Parameters fixed:",paste0(object@error.nameparameters[object@error.fix],collapse=", "),"\n")
          }
)

setMethod("showall","SaemixContinuousOutcome",
          function(object) {
            cat("Continuous outcome:",object@name,"\n")
            cat("   ",object@distribution,"distribution with ")
            if(object@error.model!="user") 
              cat(paste0(object@error.model," residual error model involving ",object@error.npar," parameter",ifelse(object@error.npar>1,"s",":")),"\n") else
                cat(paste0("user-defined error model involving ",object@error.npar," parameter",ifelse(object@error.npar>1,"s",""),":\n"))
            cat(object@error.nameparameters,"\n")
            print(object@error.function)
            cat("Initial parameter values:", paste(object@error.nameparameters,object@error.parameters,sep="="),"\n")
            if(sum(object@error.fix)>0) cat("Parameters fixed:",paste0(object@error.nameparameters[object@error.fix],collapse=", "),"\n")
          }
)

setMethod("print","SaemixContinuousOutcome",
          function(x,nlines=10,...) {
            show(x)
          }
)

# Discrete Outcome
setMethod("show","SaemixDiscreteOutcome",
          function(object) {
            cat("Discrete outcome:",object@name,"\n")
            cat("    type:", object@type,"\n")
#            cat("    distribution:", object@distribution,"\n")
          }
)
setMethod("showall","SaemixDiscreteOutcome",
          function(object) {
            cat("Discrete outcome:",object@name)
            if(object@type!="event") cat(" with",length(object@levels),"categories:",levels(object@levels),"\n") else cat("\n")
            if(object@distribution!="") cat("    distribution:", object@distribution,"\n")
          }
)

setMethod("print","SaemixDiscreteOutcome",
          function(x,nlines=10,...) {
            showall(x)
          }
)

# Survival Outcome
setMethod("show","SaemixEventOutcome",
          function(object) {
            cat("Time-to-event outcome:",object@name,"\n")
            #            cat("    distribution:", object@distribution,"\n")
          }
)
setMethod("showall","SaemixEventOutcome",
          function(object) {
            cat("Time-to-event outcome:",object@name,"\n")
            if(object@intervalCensored) 
              cat("    interval censored with length",object@intervalLength,"\n")
            cat("    distribution:", object@distribution,"\n")
          }
)

setMethod("print","SaemixEventOutcome",
          function(x,nlines=10,...) {
            showall(x)
          }
)

setMethod("show","SaemixOutcome",
          function(object) {
            cat("Outcome ")
            if(object@name!="") cat(object@name,"- ") else cat("of ")
            cat("type:",object@type)
            if(object@unit!="") cat(paste0(" (units: ",object@unit,")\n"))
          }
)

setMethod("print","SaemixOutcome",
          function(x,nlines=10,...) {
            show(x)
          }
)

################################################################
# functions to create continuous and discrete outcomes as lists

#' @examples
#' conc <- saemixOutcome("concentration", unit="mg/L")
#' conc
#' pain <- saemixOutcome("pain score",type="categorical", unit="(-)")
#' pain
#' resp <- saemixOutcome("response",type="binary")
#' resp

saemixOutcome<-function(name="",type="continuous", unit="") {
  new(Class = "SaemixOutcome", name=name, type=type, unit=unit)
}


continuousOutcome<-function(model="constant", start=c(1), error.function=constantErrorModel, npar=1, fix=c(FALSE)) {
  if(model=="user") {
    ierr<-0
    if(identical(error.function, constantErrorModel)) {
      model<-"constant"
      ierr<-1
    }
    if(identical(error.function, proportionalErrorModel)) {
      model<-"proportional"
      ierr<-1
    }
    if(identical(error.function, combined1ErrorModel)) {
      model<-"combined1"
      ierr<-1
    }
    if(identical(error.function, combined2ErrorModel)) {
      model<-"combined2"
      ierr<-1
    }
    if(identical(error.function, powerErrorModel)) {
      model<-"power"
      ierr<-1
    }
    if(ierr==1) 
      message("User error model set to ",model," for consistency with the error.function given\n")
  }
  outlist<-list(type="continuous", error.model=model, start=start, error.function=error.function, error.npar=npar, error.fix=fix)
  if(outlist$error.model=="constant" | outlist$error.model=="exponential") {
    outlist$error.function<-constantErrorModel
    outlist$error.npar<-1
    outlist$error.fix<-outlist$error.fix[1]
    outlist$start<-outlist$start[1]
    if(length(outlist$start)==0) outlist$start<-c(1)
    if(outlist$start<=0) {
      outlist$start<-1
      message("Starting value for constant error model must be positive, changed to 1\n")
    }
  }
  if(outlist$error.model=="proportional") {
    outlist$error.npar<-1
    outlist$error.function<-proportionalErrorModel
    outlist$start<-outlist$start[1]
    outlist$error.fix<-outlist$error.fix[1]
    if(length(outlist$start)==0) outlist$start<-c(0.5)
    if(outlist$start<=0) {
      outlist$start<-0.5
      message("Starting value for proportional error model must be positive, changed to 0.5\n")
    }
  }
  if(outlist$error.model=="combined1") {
    outlist$error.function<-combined1ErrorModel
    outlist$error.npar<-2
    if(length(outlist$start)>2) outlist$start<-abs(outlist$start[1:2])
    if(length(outlist$start)==1) outlist$start<-c(abs(outlist$start),0.5)
    if(length(outlist$start)==0) outlist$start<-c(1,0.5)
    if(length(outlist$error.fix)<2) outlist$error.fix<-rep(outlist$error.fix, length.out=2) else outlist$error.fix<-outlist$error.fix[1:2]
    x<-try(outlist$error.function(c(1), outlist$start))
    if(is.na(x)) message("Check starting value for error function, returns NA for x=1\n")
  }
  if(outlist$error.model=="combined2") {
    outlist$error.function<-combined2ErrorModel
    outlist$error.npar<-2
    if(length(outlist$start)>2) outlist$start<-abs(outlist$start[1:2])
    if(length(outlist$start)==1) outlist$start<-c(abs(outlist$start),0.5)
    if(length(outlist$start)==0) outlist$start<-c(1,0.5)
    if(length(outlist$error.fix)<2) outlist$error.fix<-rep(outlist$error.fix, length.out=2) else outlist$error.fix<-outlist$error.fix[1:2]
    x<-try(outlist$error.function(c(1), outlist$start))
    if(is.na(x)) message("Check starting value for error function, returns NA for x=1\n")
  }
  if(outlist$error.model=="power") {
    outlist$error.function<-powerErrorModel
    outlist$error.npar<-3
    if(length(outlist$start)>3) outlist$start<-outlist$start[1:3]
    if(length(outlist$start)==0) outlist$start<-c(1,0.5,2)
    if(length(outlist$start)==1) outlist$start<-c(outlist$start,0.5,2)
    if(length(outlist$start)==2) outlist$start<-c(outlist$start,2)
    outlist$start[1:2]<-abs(outlist$start[1:2]) # a and b must be positive
    if(length(outlist$error.fix)<3) outlist$error.fix<-rep(outlist$error.fix, length.out=3) else outlist$error.fix<-outlist$error.fix[1:3]
    x<-try(outlist$error.function(c(1), outlist$start))
    if(is.na(x)) message("Check starting value for error function, returns NA for x=1\n")
    }
  if(outlist$error.model=="user") {
    if(length(start)!=npar) outlist$error.npar<-length(start)
    x<-try(error.function(c(1), start))
    if(length(outlist$error.fix)<outlist$error.npar) outlist$error.fix<-rep(outlist$error.fix, length.out=outlist$error.npar) else outlist$error.fix<-outlist$error.fix[1:outlist$error.npar]
    if(is.na(x)) {
      message("Check user-defined error function, returns NA for x=1\n")
      return()
    }
  }
  return(outlist)
}

discreteOutcome<-function(type="binary", distribution="", ncat=NULL, levels=NULL, maxEvents=NULL, intervalCensored=FALSE, intervalLength=NULL) {
  if(is.null(levels) & !is.null(ncat)) levels<-c(1:ncat)
  outlist<-list(type=type, distribution=distribution, levels=levels, maxEvents=maxEvents, intervalCensored=intervalCensored, intervalLength=intervalLength)
  if(outlist$distribution=="") {
    if(type=="categorical") outlist$distribution<-"categorical"
    if(type=="binary") outlist$distribution<-"bernouilli"
  }
  if(is.null(levels) & type=="binary")
    outlist$levels<-c(0,1)
  if(outlist$type=="event") {
    if(is.null(outlist$maxEvents)) outlist$maxEvents<-1
  }
  if(is.null(outlist$intervalLength)) outlist$intervalLength<-0
  return(outlist)
}
