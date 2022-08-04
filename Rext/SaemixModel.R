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
    covariates="list",	# the list of covariates defined in the parameters
    covariate.model="matrix",	# covariate model
    covariate.model.fix="matrix",	# fixed parameters in the covariate model (defaults to all 0's)
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
    if(missing(mu.start) || length(mu.start)==0 || !is.numeric(mu.start)) {
      mu.start<-rep(1,npar)
    } else {
      if(length(mu.start)!=npar & verbose) message("Size mismatch, resizing mu.start")
      length(mu.start)<-npar
      mu.start[is.na(mu.start)]<-1
    }
    if(missing(transform.par) || length(transform.par)==0 || !is.numeric(transform.par)) {
      transform.par<-rep(0,npar)
    } else {
      if(length(transform.par)!=npar & verbose) message("Size mismatch, resizing transform.par")
      length(transform.par)<-npar
      transform.par[is.na(transform.par)]<-0
    }
    if(missing(mu.fix) || length(mu.fix)==0) {
      mu.fix<-rep(0,npar)
    } else {
      if(is.logical(mu.fix)) mu.fix<-as.integer(mu.fix)
      mu.fix<-as.integer(mu.fix>0)
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
            "covariates"={return(x@covariates)},
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
            "covariates"={x@covariates<-value},
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
  if(!is.null(covariate.model)) 
    covariate.model<-validate.covariate.model(covariate.model, npar=npar, verbose=verbose)
  if(!is.null(covariate.model)) { 
      colnames(covariate.model)<-.Object@name.modpar
      .Object@covariate.model<-covariate.model
  }
  if(!is.null(covariate.model.fix)) 
    covariate.model.fix<-validate.covariate.model(covariate.model.fix, npar=npar, verbose=verbose)
  if(!is.null(covariate.model.fix)) {
    colnames(covariate.model.fix)<-.Object@name.modpar
    if(is.null(covariate.model)) {
      if(verbose) message("covariate.model missing but covariate.model.fix given, setting the same structure\n")
      covariate.model<-.Object@covariate.model<-covariate.model.fix
    }
    .Object@covariate.model.fix<-covariate.model.fix
    }
  if(is.null(covariate.model)) {
    .Object@betaest.model<-matrix(rep(1,npar), nrow=1, dimnames=list(c("mu.psi"),.Object@name.modpar))
    beta.start<-numeric()
  } else {
    .Object@betaest.model<-rbind(rep(1,npar),covariate.model)
    colnames(.Object@betaest.model)<-.Object@name.modpar
    nb.beta<-sum(.Object@covariate.model)
    if(missing(beta.start)) beta.start<-rep(0,nb.beta)
    if(nb.beta>length(beta.start)) 
      beta.start<-c(beta.start, rep(0,nb.beta-length(beta.start)))
    beta.start<-beta.start[1:nb.beta] # if beta.start longer than the number of par-cov relationships
  }
  .Object@beta.start<-beta.start
  return(.Object)
}

validate.covariate.model <- function(x, npar, verbose=TRUE){
  if(dim(x)[2]!=npar || !(is.numeric(x))) {
    if(verbose) message("Size mismatch between covariate model and the number of parameters, ignoring")
    x<-NULL
  } else {
    if(sum(!(c(x) %in% c(0,1)))>0) {
      if(verbose) message("covariate model must contain only 0/1, ignoring")
      x<-NULL
    }
  }
  return(x)
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
    ind.model = "list", # list of individual variability models (built when associating covariate.model, covariance.model and data) corresponding the the variability levels
    nphirep = "list", # a list, each element is a vector giving how many times do we need to repeat the phi of each variability level
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
            "nvarlevel"={return(x@nvarlevel)},
            "var.model"={return(x@var.model)},
            "ind.model"={return(x@ind.model)},
            "nphirep"={return(x@nphirep)},
            "name.predictors"={return(x@name.predictors)},
            "name.X"={return(x@name.X)},
            "name.cov"={return(x@name.cov)},
            "name.thetas"={return(x@name.thetas)},
            "Mcovariates"={return(x@Mcovariates)},
            "npar"={return(x@npar)},
            "mu.start"={return(x@mu.start)},
            "mu.fix"={return(x@mu.fix)},
            "transform.par"={return(x@transform.par)},
            "name.modpar"={return(x@name.modpar)},
            "covariates"={return(x@covariates)},
            "covariate.model"={return(x@covariate.model)},
            "covariate.model.fix"={return(x@covariate.model.fix)},
            "beta.start"={return(x@beta.start)},
            "betaest.model"={return(x@betaest.model)},
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
  signature = "SaemixModel" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "nvarlevel"={x@nvarlevel<-value},
            "var.model"={x@var.model<-value},
            "ind.model"={x@ind.model<-value},
            "nphirep"={x@nphirep<-value},
            "name.predictors"={x@name.predictors<-value},
            "name.X"={x@name.X<-value},
            "name.cov"={x@name.cov<-value},
            "name.thetas"={x@name.thetas<-value},
            "Mcovariates"={x@Mcovariates<-value},
            "npar"={x@npar<-value},
            "mu.start"={x@mu.start<-value},
            "mu.fix"={x@mu.fix<-value},
            "transform.par"={x@transform.par<-value},
            "name.modpar"={x@name.modpar<-value},
            "covariates"={x@covariates<-value},
            "covariate.model"={x@covariate.model<-value},
            "covariate.model.fix"={x@covariate.model.fix<-value},
            "beta.start"={x@beta.start<-value},
            "betaest.model"={x@betaest.model<-value},
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

######################################################################################
####			S4 methods: show, print, showall	
####			SaemixModel class - method to print/show data		
####################################################################################

#' @rdname print-methods
#' @exportMethod print

setMethod("print","SaemixModel",
          function(x,...) {
            cat("Nonlinear mixed-effects model\n")
            if( is.null(body(x@model))) {
              cat("No model function set yet\n")
              return()
            }
            cat("  Model function")
            if(length(x@description)>0 && nchar(x@description)>0) cat(": ",x@description)
            cat("\n")
            # Responses
            nout<-length(x@outcome)
            cat(paste0("  ","with ",nout," outcome",ifelse(nout>1,"s: ",": ")), x@name.outcome,"\n")
            # Parameters
            distrib<-c("normal","log-normal","probit","logit")
            cat("  Nb of parameters:",x@npar,"\n")
            cat("      parameter names: ",x@name.modpar,"\n")
            cat("      distribution:\n")
            tab<-cbind(Parameter=x@name.modpar,Distribution=distrib[x@transform.par+1], Estimated=ifelse(x@mu.fix==0,"estim","fixed"), mu.CI=x@mu.start, omega.CI=diag(x@var.model[[1]]@omega))
            if(x@nvarlevel>1) {
              for(i in 2:x@nvarlevel) {
                tab<-cbind(tab, diag(x@var.model[[i]]@omega)) 
                colnames(tab)[ncol(tab)]<-paste0("varlevel",i,".CI")
              }
            }
            print(tab,quote=FALSE)
            if(length(x@covariate.model)>0) {
              cat("Covariate model:\n")
              print(x@covariate.model)
              if(length(x@covariate.model.fix)>0) {
                cat("Fixed parameters in covariate model:\n")
                print(x@covariate.model.fix)
              }
            }
          }
)
# distrib<-c("normal","log-normal","probit","logit")
# cat("  Model type")
# if(length(x@modeltype)>0 && nchar(x@modeltype)>0) cat(": ",x@modeltype)
# cat("\n")
# print(x@model)
# cat("  Variance-covariance matrix:\n")
# tab<-x@covariance.model
# #    try(colnames(tab)<-rownames(tab)<-x@name.modpar)
# print(tab,quote=FALSE)
# st1<-paste(x@name.sigma,x@error.init,sep="=")
# if (x@modeltype=="structural"){
#   cat("  Error model:",x@error.model,", initial values:",st1[x@indx.res],"\n")
# }
# if(dim(x@covariate.model)[1]>0) {
#   cat("  Covariate model:")
#   if(sum(x@covariate.model)==0) cat(" none\n") else {
#     cat("\n")
#     print(x@covariate.model)
#   }
# } else cat("    No covariate in the model.\n")

#' @rdname show-methods
#' @exportMethod show

setMethod("show","SaemixModel",
          function(object) {
            cat("Nonlinear mixed-effects model\n")
            if(length(object@description)>0 && nchar(object@description)>0)
              cat(": ",object@description,"\n")
            if( is.null(body(object@model))) {
              cat("No model function set yet\n")
              return()
            } 
            # Parameters
            cat("  ",object@npar,"parameter(s): ",object@name.modpar,"\n")
            if(dim(object@covariate.model)[1]>0) {
              cat("     covariate model with",sum(object@covariate.model),"covariate effects\n")
            }
            # Responses
            nout<-length(object@outcome)
            cat(paste0("   ","with ",nout," outcome",ifelse(nout>1,"s: ",": ")), object@name.outcome,"\n")
          }
)

#' @rdname showall-methods
#' @exportMethod showall

setMethod("showall","SaemixModel",
          function(object) {
            cat("Nonlinear mixed-effects model\n")
            if(length(object@description)>0 && nchar(object@description)>0)
              cat(": ",object@description,"\n")
            if( is.null(body(object@model))) {
              cat("No model function set yet\n")
              return()
            } 
            # Parameters
            distrib<-c("normal","log-normal","probit","logit")
            cat("  ",object@npar,"parameter(s): ",object@name.modpar,"\n")
            if(dim(object@covariate.model)[1]>0) {
              cat("     covariate model with",sum(object@covariate.model),"covariate effects\n")
            }
            # Responses
            nout<-length(object@outcome)
            cat(paste0("   ","with ",nout," outcome",ifelse(nout>1,"s: ",":")),"\n")
            for(iout in 1:nout) showall(x@outcome[[iout]])
            if(length(object@name.predictors)>0) cat("      predictors: ",object@name.predictors,"\n")
            if(length(object@name.X)>0) cat("      predictor used for X-axis on graphs: ",object@name.X,"\n")
            # Items in object
            if(length(object@name.thetas)>0) cat("      fixed parameters: ",object@name.thetas,"\n")
            if(length(object@name.cov)>0) cat("      covariates: ",object@name.cov,"\n")
            cat("      distribution:\n")
            # Parameters
            tab<-cbind(Parameter=object@name.modpar,Distribution=distrib[object@transform.par+1], Estimated=ifelse(object@mu.fix==0,"estim","fixed"), mu.CI=object@mu.start, omega.CI=diag(object@var.model[[1]]@omega))
            if(object@nvarlevel>1) {
              for(i in 2:object@nvarlevel) {
               tab<-cbind(tab, diag(object@var.model[[i]]@omega)) 
               colnames(tab)[ncol(tab)]<-paste0("varlevel",i,".CI")
              }
            }
            print(tab,quote=FALSE)
            # Covariates TODO
            # if(length(object@ind.model)>0) {}
            if(dim(object@covariate.model)[1]>0) {
              cat("  Covariate model:")
              if(sum(object@covariate.model)==0) cat(" none\n") else {
                cat("\n")
                print(object@covariate.model)
              }
            } else cat("  No covariate in the model.\n")
          }
)


