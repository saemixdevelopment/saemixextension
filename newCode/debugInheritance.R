# From Genolini

setClass(
  Class="Trajectories",
  representation(times="numeric",traj="matrix"),
  validity=function(object){
    cat("~~~ Trajectories: inspector ~~~ \n")
    if(length(object@times)!=ncol(object@traj)){
      stop ("[Trajectories: validation] the number of temporal measurements does not correspond")
    }else{}
    return(TRUE)
  }
)
new(Class="Trajectories",times=1:2,traj=matrix(1:2,ncol=2))

setMethod (
  f="initialize",
  signature="Trajectories",
  definition=function(.Object,times,traj){
    cat ("~~~~~ Trajectories: initializator ~~~~~ \n")
    if(!missing(traj)){
      colnames(traj) <- paste("T",times,sep="")
      rownames(traj) <- paste("I",1:nrow(traj),sep="")
      .Object@times <- times
      .Object@traj <- traj
      validObject(.Object) # call of the inspector
    }
    return(.Object)
  }
)

setClass(
  Class="Partition",
  representation=representation (
    nbGroups="numeric",
    part="factor"
  )
)
setClass(
  Class="TrajPartitioned",
  representation=representation(listPartitions="list"),
  contains= "Trajectories"
)

setMethod("initialize","TrajPartitioned",
          function(.Object,times,traj,listPartitions){
            cat("~~~~TrajPartitioned: initializator ~~~~ \n")
            if(!missing(traj)){
              .Object@times <- times
              .Object@traj <- traj # Assignment of attributes
              .Object@listPartitions <- listPartitions
            }
            return(.Object) # return of the object
          }
)

trajPitie <- new(Class="Trajectories")
trajCochin <- new(
  Class= "Trajectories",
  times=c(1,3,4,5),
  traj=rbind (
    c(15,15.1, 15.2, 15.2),
    c(16,15.9, 16,16.4),
    c(15.2, NA, 15.3, 15.3),
    c(15.7, 15.6, 15.8, 16)
    )
  )
partCochin <- new(Class="Partition",nbGroups=2,part=factor(c("A","B","A","B")))
tdCochin <- new(
  Class="TrajPartitioned",
  traj=trajCochin@traj,
  times=c(1,3,4,5),
  listPartitions=list(partCochin,partCochin))

######## need to change return(NULL) into return(.Object)

setMethod(
  f="initialize",
  signature="SaemixStructuralModel",
  definition=function(.Object, outcome, model, description="", verbose=FALSE){
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


setMethod(
  f="initialize",
  signature="SaemixModel",
  definition=function(.Object, outcome, model, description="", psi0, beta.start, parameter, transform.par, mu.fix, covariate.model, covariate.model.fix, var.level, verbose=FALSE){
    # ECO how do we initialize a child using the parent ?
    .Object <- callNextMethod(.Object, outcome=outcome, model=model, description=description, verbose=verbose)
    if(missing(psi0) & missing(parameter)) {
      if(verbose) message("[ SaemixModel : Error ] Missing initial estimates of parameters (psi0) or parameter names (parameter)")
      return (.Object)
    }
    if(missing(outcome)) {
      if(verbose) message("No outcome given, assuming a single continuous outcome with a constant error model")
      outcome<-list(y=new(Class="SaemixContinuousOutcome"))
    }
    if(!is(outcome, "list")) { # if a single outcome is given, transform to a list
      if(!is(outcome, "SaemixOutcome")) {
        if(verbose) message("Outcomes must be given as a list of SaemixOutcome objects")
        return (.Object)
      }
      outcome<-list(outcome)
    }
    .Object@nb.outcome<-length(outcome)
    if(is.null(names(outcome))) names(outcome)<-paste0("y",1:.Object@nb.outcome)
    .Object@outcome<-outcome
    .Object@model<-model
    .Object@name.outcome<-names(outcome)
    .Object@description<-description
    print(.Object@model)
    #    print(.Object)
    
    validObject(.Object)
    return (.Object)
  }
)



