# Observation class

# types: continuous, discrete (count/categorical), event

setClass(
  Class="SaemixContinuousOutcome",
  representation=representation(
    name.outcome = "character", # outcome name
    type.outcome = "character", # Type, continuous (for compatibility with discrete outcome
    error.model = "character", # name of the error model associated
    error.npar = "numeric", # number of parameters in the error model
    error.function = "function"  # error model function
  ),
  validity=function(object){
    # Check error.model is one of constant, proportional, combined1, combined2, power or user
    # Check error.function exists if error.model=user
    #    cat ("--- Checking SaemixData object ---\n")
    return(TRUE)
  }
)

setClass(
  Class="SaemixDiscreteOutcome",
  representation=representation(
    name.outcome = "character", # outcome name
    type.outcome = "character" # Type, one of discrete, event
  ),
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
  signature="SaemixContinuousOutcome",
  definition= function (.Object, name.outcome, error.model="constant", error.npar=NULL, error.function=NULL) {
    #    cat ("--- initialising SaemixDiscrete Object --- \n")
    .Object@name.outcome<-name.outcome
    .Object@type.outcome<-"continuous"
    if(!(error.model %in% c("constant","proportional", "combined1","combined2", "user"))) {
      message("[ SaemixContinuousOutcome : validation ] Please specify a valid error model for the continuous outcome (one of constant, proportional, combined1, combined2 or user.")
      return("Error model not given")
    }
    .Object@error.model<-error.model
    if(error.model=="constant") {
      .Object@error.npar<-1
      .Object@error.function<-constantErrorModel
    }
    if(error.model=="proportional") {
      .Object@error.npar<-1
      .Object@error.function<-proportionalErrorModel
    }
    if(error.model=="combined1") {
      .Object@error.npar<-2
      .Object@error.function<-combined1ErrorModel
    }
    if(error.model=="combined2") {
      .Object@error.npar<-2
      .Object@error.function<-combined2ErrorModel
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
    }
    return (.Object )
  }
)

setMethod(
  f="initialize",
  signature="SaemixDiscreteOutcome",
  definition= function (.Object, name.outcome, type.outcome) {
    #    cat ("--- initialising SaemixDiscrete Object --- \n")
        .Object@name.outcome<-name.outcome
        .Object@type.outcome<-type.outcome
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
combined2ErrorModel<-function(f,ab) {
  #  g<-cutoff(ab[1]+ab[2]*abs(f))
  g<-cutoff(sqrt(ab[1]^2+ab[2]^2*f^2))  # Johannes 02/21
  return(g)
}


# Getteur
setMethod(
  f ="[",
  signature = "SaemixDiscreteOutcome" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.outcome"={return(x@name.outcome)},
            "type.outcome"={return(x@type.outcome)},
            stop("No such attribute\n")
    )
  }
)

setMethod(
  f ="[",
  signature = "SaemixContinuousOutcome" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.outcome"={return(x@name.outcome)},
            "type.outcome"={return(x@type.outcome)},
            "error.model"={return(x@error.model)},
            "error.npar"={return(x@error.npar)},
            "error.function"={return(x@error.function)},
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
            "name.outcome"={x@name.outcome<-value},
            "type.outcome"={x@type.outcome<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

setReplaceMethod(
  f ="[",
  signature = "SaemixContinuousOutcome" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.outcome"={x@name.outcome<-value},
            "type.outcome"={x@type.outcome<-value},
            "error.model"={x@error.model<-value},
            "error.npar"={x@error.npar<-value},
            "error.function"={x@error.function<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

# Show

setMethod("show","SaemixContinuousOutcome",
          function(object) {
            cat("Continuous outcome:",object@name.outcome,"\n")
            cat("    gaussian distribution with ")
            if(object@error.model!="user") 
              cat(object@error.model,"residual error model involving",object@error.npar,"parameter",ifelse(object@error.npar>1,"s",""),"\n") else {
                cat("user-defined error model involving",object@error.npar,"parameters:\n")
                print(object@error.function)
              }
          }
)

setMethod("showall","SaemixContinuousOutcome",
          function(object) {
            cat("Continuous outcome:",object@name.outcome,"\n")
            cat("    gaussian distribution with ")
            if(object@error.model!="user") 
              cat(object@error.model,"residual error model involving",object@error.npar,"parameter",ifelse(object@error.npar>1,"s",""),":\n") else
                cat("user-defined error model involving",object@error.npar,"parameters:\n")
            print(object@error.function)
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
          }
)

setMethod("print","SaemixDiscreteOutcome",
          function(x,nlines=10,...) {
            show(x)
          }
)
