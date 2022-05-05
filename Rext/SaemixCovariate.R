# Observation class - generic
#' Class "SaemixCovariateType"
#' 
#' An object of the SaemixCovariateType class, representing a continuous or categorical(/binary) covariate
#' 
#' The SaemixCovariateType is the parent class, containing only the name, type and unit of the covariate, for use 
#' in the creation of a data object in saemix.
#' Two child classes have also been defined, for use in the model object in saemix:
#' - SaemixContinuousCovariate for continuous covariates
#' - SaemixDiscreteCovariate for discrete (categorical/binary) covariates
#' These classes contain transformation models 
#' - for the continuous covariates, the user can specify a function to apply to the vector of values 
#' and either a centering function (eg median) or a value (eg 60)
#' (example of use: transform a covariate to the logarithm of the covariate, centering with respect to the median)
#' - for categorical covariates, the grouping and reference values 
#' (example of use: regroup CT and TT alleles and define CC as reference category)
#' Covariate transformations will be applied to the covariate in the dataset once combined with a data object
#' 
#' @name SaemixCovariateType-class 
#' @docType class
#' @aliases SaemixCovariateType SaemixCovariateType-class 
#' @aliases print,SaemixCovariateType showall,SaemixCovariateType show,SaemixCovariateType
#' @aliases SaemixDiscreteCovariate SaemixDiscreteCovariate-class 
#' @aliases print,SaemixDiscreteCovariate showall,SaemixDiscreteCovariate show,SaemixDiscreteCovariate
#' @aliases SaemixContinuousCovariate SaemixContinuousCovariate-class 
#' @aliases print,SaemixContinuousCovariate showall,SaemixContinuousCovariate show,SaemixContinuousCovariate
#' 
#' @examples
#' showClass("SaemixCovariateType")
#' 
#' @keywords classes
#' @exportClass SaemixCovariateType

setClass(Class = "SaemixCovariateType",
         representation=representation(
           name = "character", # outcome name
           type = "character", # Type: continuous, binary, categorical
           unit = "character" # Covariate unit
         ),
         validity=function(object){
           if(!(object@type %in% c("continuous","binary","categorical"))) {
             message("Covariate type should be one of continuous, binary, categorical")
             return(NULL)
         }
           return(TRUE)
         }
)

setMethod( 
  f="initialize",
  signature="SaemixCovariateType",
  definition=function(.Object, name="", type="continuous", unit=""){
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
  signature = "SaemixCovariateType" ,
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
  signature = "SaemixCovariateType" ,
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
setClass(Class = "SaemixCovariate",
         contains="SaemixCovariateType",
         representation=representation(
           covariate.transform = "SaemixCovariateTransform", # Transformation applied to the covariate (a function)
           beta="numeric",  # covariate effect parameter
           beta.fix="numeric" # 1 if corresponding element of beta is fixed, 0 if estimated (default)
         ),
         validity=function(object){
           return(TRUE)
         }
)



# Getteur
setMethod(
  f ="[",
  signature = "SaemixCovariate" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "covariate.transform"={return(x@covariate.transform)},
            "beta"={return(x@beta)},
            "beta.fix"={return(x@beta.fix)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixCovariate" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "covariate.transform"={x@covariate.transform<-value},
            "beta"={x@beta<-value},
            "beta.fix"={x@beta.fix<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)


setMethod( 
  f="initialize",
  signature="SaemixCovariate",
  definition=function(.Object, name="", type="continuous", unit="", beta=1, beta.fix=0, covariate.transform){
    .Object@name <- name
    .Object@type <- type
    .Object@unit <- unit
    if(missing(covariate.transform) || !is(covariate.transform,"SaemixCovariateTransform")) {
      if(type=="continuous") covariate.transform<-new(Class="covmodelCont2Cont",name=name) else covariate.transform<-new(Class="covmodelCat2Cat",name=name)
    }
    covariate.transform@name<-name
    .Object@covariate.transform<-covariate.transform
    .Object@type<-covariate.transform@type
    if(covariate.transform@type=="categorical" && length(covariate.transform@ncat)>0) {
      length(beta)<-length(beta.fix)<-(covariate.transform@ncat-1)
      beta[is.na(beta)]<-1
      beta.fix[is.na(beta.fix)]<-0
    } 
    .Object@beta<-beta
    .Object@beta.fix<-beta.fix

    validObject(.Object)
    return(.Object)
  }
)

########################################################################
# Print and show methods

setMethod("show","SaemixCovariateType",
          function(object) {
            cat(object@type,"covariate")
            if(object@name!="") cat(": ",object@name)
            if(object@unit!="") cat(paste0(" (unit: ",object@unit,")"))
            cat("\n")
          }
)

setMethod("print","SaemixCovariateType",
          function(x,nlines=10,...) {
            show(x)
          }
)


# SaemixCovariate
setMethod("show","SaemixCovariate",
          function(object) {
            cat(object@type,"covariate")
            if(object@name!="") cat(": ",object@name)
            if(object@unit!="") cat(paste0(" (unit: ",object@unit,")"))
            cat("\n")
          }
)


setMethod("print","SaemixCovariate",
          function(x,nlines=10,...) {
            show(x)
          }
)

setMethod("showall","SaemixCovariate",
          function(object) {
            # cat(object@type,"covariate")
            # if(object@name!="") cat(": ",object@name)
            # if(object@unit!="") cat(paste0(" (unit: ",object@unit,")"))
            # cat("\n")
            showall(object@covariate.transform)
          }
)

########################################################################
# Creator functions
contCov <- function(name, transform, centering, beta=0, beta.fix=0, verbose=FALSE) {
  if(missing(name)) name<-"cov"
  if(!missing(centering)) {
    #    if(is(centering, "character") && is(as.function(centering),"function")) centering<-as.function(centering)
    #    if(is(centering, "character") && is(as.double(centering),"numeric")) centering<-as.double(centering)
    if(is(centering, "function")) {
      covmodel <- new(Class="covmodelCont2Cont", name=name, transform.function=transform, centering.function=centering, verbose=verbose)
    } else {
      if(is(centering,"numeric")) {
        covmodel <- new(Class="covmodelCont2Cont", name=name, transform.function=transform, centering.value=centering, verbose=verbose)
      } else {
        if(verbose) message("The argument centering must be either a function or a numeric value")
        covmodel <- new(Class="covmodelCont2Cont", name=name, transform.function=transform, verbose=verbose)
      }
    }
  } else 
    covmodel <- new(Class="covmodelCont2Cont", name=name, transform.function=transform, verbose=verbose)
  xcov <- new(Class="SaemixCovariate", name=name, type="continuous", covariate.transform=covmodel)
  if(!is.na(beta)) xcov@beta<-beta
  xcov@beta.fix<-as.integer(beta.fix!=0)
  return(xcov)
}

# special kind of catCov with ncat=2
binCov <- function(name, breaks, groups, name.cat=character(), reference=character(),beta=0, beta.fix=0, verbose=FALSE) {
  xcov<-catCov(name=name, ncat=2, breaks=breaks, groups=groups, name.cat=name.cat, reference=reference,beta=beta, beta.fix=beta.fix, verbose=verbose)
  xcov@covariate.transform@name.cat<-name.cat # remove the names if none given, otherwise catCov sets them to G1, G2, ...
  xcov@covariate.transform@reference<-reference # remove the reference if none given, otherwise catCov sets it to G1
  return(xcov)
}

catCov <- function(name, breaks, groups, name.cat=character(), reference=character(), ncat=0, beta=0, beta.fix=0, verbose=FALSE) {
  if(missing(name)) name<-"cov"
  if(!missing(breaks) & !missing(groups)) { # If both groups and breaks are given, we can't decide, return empty
    if(verbose) message("Please specify only one of either 'groups' if the original covariate is categorical or 'breaks' if the original covariate is continuous and should be categorised")
    return(NULL)
  }
  if(ncat==0 & length(name.cat)==0 & missing(breaks) & missing(groups)) { # nothing given, we assume a categorical covariate unchanged
    if(verbose) message("No information given, assuming a categorical covariate unchanged")
    covmodel<-new(Class="covmodelCat2Cat", name=name, reference=reference, verbose=verbose)
  #  if(missing(name.cat)) name.cat<-covmodel@name.cat
  } else {
    if(missing(breaks) & missing(groups)) { # Again nothing given, we assume a categorical covariate
      if(verbose) message("No information given, assuming a categorical covariate unchanged")
#      if(!missing(name.cat) & length(name.cat)>0) ncat<-length(name.cat)
      if(length(name.cat)==0) name.cat<-paste0("G",1:ncat)
      covmodel<-new(Class="covmodelCat2Cat", name=name, name.cat=name.cat, reference=reference, verbose=verbose)
    }
    if(!missing(groups)) { # original covariate is categorical
      covmodel<-new(Class="covmodelCat2Cat", name=name, name.cat=name.cat, groups=groups, reference=reference, verbose=verbose)
    }
    if(!missing(breaks)) { # original covariate is continuous
      covmodel<-new(Class="covmodelCont2Cat", name=name, name.cat=name.cat, breaks=breaks, reference=reference, verbose=verbose)
    }
  }
  xcov <- new(Class="SaemixCovariate", name=name, type="categorical", covariate.transform=covmodel)
  if(!is.na(beta[1]) & length(xcov@covariate.transform@ncat)>0) {
    length(beta)<-xcov@covariate.transform@ncat-1
    beta[is.na(beta)]<-1
  }
  xcov@beta<-beta
  if(!is.na(beta.fix[1]) & length(xcov@covariate.transform@ncat)>0) {
    length(beta.fix)<-xcov@covariate.transform@ncat-1
    beta.fix[is.na(beta.fix)]<-0
  }
  xcov@beta.fix<-as.integer(beta.fix!=0)
  return(xcov)
}
########################################################################
