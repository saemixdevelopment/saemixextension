# ToDo
# add documentation + aliases and formulas for the error models

setClass(
  Class="SaemixErrorModel",
  representation=representation(
    name = "character", # name of the error model
    log = "character", # warning messages
    npar = "numeric", # number of parameters in the error model
    param = "numeric", # value of parameters in the error model
    param.names = "character", # names of parameters in the error model
    model = "function"  # error model function
  ),
  validity=function(object){
    # Check model is one of constant, proportional, combined1, combined2, power or user
    # Check function exists if model=user
    #    cat ("--- Checking SaemixData object ---\n")
    return(TRUE)
  }
)

# ToDo: add checks for size of objects param and param.names

setMethod(
  f="initialize",
  signature="SaemixErrorModel",
  definition= function (.Object, name="constant", npar=NULL, model=NULL, param=NULL, param.names=NULL, unit="", verbose=FALSE) {
    #    cat ("--- initialising SaemixDiscrete Object --- \n")
    logmsg<-""
    if(!(name %in% c("constant","proportional", "combined1","combined2","power", "user"))) {
      if(verbose) message("[ SaemixContinuousOutcome : validation ] Please specify a valid error model for the continuous outcome (one of constant, proportional, combined1, combined2, power or user.")
      return("Invalid error model name")
    }
    .Object@name<-name
    if(name=="constant") {
      .Object@npar<-1
      .Object@model<-constantErrorModel
      if(is.null(param)) .Object@param<-c(1) else .Object@param<-param[1]
      if(!is.null(param.names)) .Object@param.names<-param.names[1] else .Object@param.names<-c("a")
    }
    if(name=="proportional") {
      .Object@npar<-1
      .Object@model<-proportionalErrorModel
      if(is.null(param)) .Object@param<-c(1) else .Object@param<-param[1]
      if(!is.null(param.names)) .Object@param.names<-param.names[1] else .Object@param.names<-c("b")
    }
    if(name=="combined1") {
      .Object@npar<-2
      .Object@model<-combined1ErrorModel
      if(is.null(param)) .Object@param<-c(1,1) else .Object@param<-param[1:2]
      if(!is.null(param.names)) .Object@param.names<-param.names[1:2] else .Object@param.names<-c("a","b")
    }
    if(name=="combined2") {
      .Object@npar<-2
      .Object@model<-combined2ErrorModel
      if(is.null(param)) .Object@param<-c(1,1) else .Object@param<-param[1:2]
      if(!is.null(param.names)) .Object@param.names<-param.names[1:2] else .Object@param.names<-c("a","b")
    }
    if(name=="power") {
      .Object@npar<-3
      .Object@model<-powerErrorModel
      if(is.null(param)) .Object@param<-c(1,1,2) else .Object@param<-param[1:3]
      if(!is.null(param.names)) .Object@param.names<-param.names[1:3] else .Object@param.names<-c("a","b","c")
    }
    if(name=="user") {
      if(is.null(npar) | is.null(model) | !(is.function(model))) {
        if(verbose) message("[ SaemixContinuousOutcome : validation ] When specifying a user-defined model, please give the number of param and a valid model function with 2 arguments, f and ab, where the number of param is the size of the vector ab.")
        return("Error model or number of param not given")
      }
      if(!identical(names(formals(model)),c("f","ab"))) {
        if(verbose) message("[ SaemixContinuousOutcome : validation ] The error model should be a function with 2 arguments, f and ab, where the number of parameters is the size of the vector ab.")
        return("Error model function arguments mismatch")
      }
      .Object@npar<-npar
      .Object@model<-model
      if(is.null(param.names)) param.names<-letters[1:npar]
      .Object@param.names<param.names
      if(is.null(param)) .Object@param<-rep(1,npar) else .Object@param<-param[1:npar]
    }
    return (.Object )
  }
)

setMethod(
  f ="[",
  signature = "SaemixErrorModel" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name"={return(x@name)},
            "npar"={return(x@npar)},
            "param.names"={return(x@param.names)},
            "param"={return(x@param)},
            "model"={return(x@model)},
            stop("No such attribute\n")
    )
  }
)

setReplaceMethod(
  f ="[",
  signature = "SaemixErrorModel" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name"={x@name<-value},
            "npar"={x@npar<-value},
            "param.names"={x@param.names<-value},
            "param"={x@param<-value},
            "model"={x@model<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

#################################################################################

#' @rdname saemix.internal
#' 
#' @aliases cutoff cutoff.eps cutoff.max cutoff.res
#' @aliases normcdf norminv
#' 
#' @keywords internal

cutoff<-function(x,seuil=.Machine$double.xmin) {x[x<seuil]<-seuil; return(x)}
cutoff.max<-function(x) max(x,.Machine$double.xmin)
cutoff.eps<-function(x) max(x,.Machine$double.eps)
cutoff.res<-function(x,ares,bres) max(ares+bres*abs(x),.Machine$double.xmin)


# Predefined error model functions
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

#################################################################################
# print/show functions

setMethod("print","SaemixContinuousOutcome",
          function(x,nlines=10,...) {
            show(x)
          }
)

setMethod("show","SaemixErrorModel",
          function(object) {
            cat("Residual error model with variance model",object@name,"\n")
          }
)

setMethod("showall","SaemixErrorModel",
          function(object) {
            cat("Residual error model with variance model",object@name,"\n")
            xout<-data.frame(Parameter=object@param.names, CI=object@param)
            print(xout, row.names=FALSE)
            cat("Error model function:\n")
            print(object@model)
          }
)
