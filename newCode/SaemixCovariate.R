# Generic, move to aaa_generics.R
setGeneric(name="transformCovariate",
           def=function(covariateModel, x=NULL) standardGeneric("transformCovariate")
)

# Observation class - generic
#' Class "SaemixCovariate"
#' 
#' An object of the SaemixCovariate class, representing a continuous or categorical(/binary) covariate
#' 
#' The SaemixCovariate is the parent class, containing only the name, type and unit of the covariate, for use 
#' in the creation of a data object in saemix.
#' Two child classes have also been defined, for use in the model object in saemix:
#' - SaemixContinuousCovariate for continuous covariates
#' - SaemixDiscreteCovariate for discrete (categorical/binary) covariates
#' These classes contain transformations for the continuous covariates, with the option to center them
#' (for instance, to transform a covariate to the logarithm of the covariate, centering with respect to the median
#' or to a given value), for categorical covariates, the grouping and reference values.
#' 
#' @name SaemixCovariate-class 
#' @docType class
#' @aliases SaemixCovariate SaemixCovariate-class 
#' @aliases print,SaemixCovariate showall,SaemixCovariate show,SaemixCovariate
#' @aliases SaemixDiscreteCovariate SaemixDiscreteCovariate-class 
#' @aliases print,SaemixDiscreteCovariate showall,SaemixDiscreteCovariate show,SaemixDiscreteCovariate
#' @aliases SaemixContinuousCovariate SaemixContinuousCovariate-class 
#' @aliases print,SaemixContinuousCovariate showall,SaemixContinuousCovariate show,SaemixContinuousCovariate
#' 
#' @examples
#' showClass("SaemixCovariate")
#' 
#' @keywords classes
#' @exportClass SaemixCovariate

setClass(Class = "SaemixCovariate",
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
  signature="SaemixCovariate",
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
  signature = "SaemixCovariate" ,
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
  signature = "SaemixCovariate" ,
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
# Continuous covariate
setClass(Class = "SaemixContinuousCovariate",
         contains="SaemixCovariate",
         representation=representation(
           transform.function = "function", # Transformation applied to the covariate (a function)
           centering.function = "function", # if given, the function is applied (eg: median, mean)
           centering.value = "numeric" # if given, 
         ),
         validity=function(object){
           return(TRUE)
         }
)

setMethod( 
  f="initialize",
  signature="SaemixContinuousCovariate",
  definition=function(.Object, name="", transform.function="log", centering.function, centering.value){
    .Object <- callNextMethod(.Object, name, type="continuous")
    if(!missing(transform.function)) .Object@transform.function <- transform.function else .Object@transform.function<-function(x) x
    if(!missing(centering.function)) .Object@centering.function<-centering.function else .Object@centering.function<-function(x) 1
    if(!missing(centering.value)) .Object@centering.value<-centering.value
    validObject(.Object)
    return(.Object)
  }
)

# Getteur
setMethod(
  f ="[",
  signature = "SaemixContinuousCovariate" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "transform.function"={return(x@transform.function)},
            "centering.function"={return(x@centering.function)},
            "centering.value"={return(x@centering.value)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixContinuousCovariate" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "transform.function"={x@transform.function<-value},
            "centering.function"={x@centering.function<-value},
            "centering.value"={x@centering.value<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

# Discrete covariate
setClass(Class = "SaemixDiscreteCovariate",
         contains="SaemixCovariate",
         representation=representation(
           groups = "list", # groups for categorical category
           reference = "character" # reference category
         ),
         validity=function(object){
           return(TRUE)
         }
)

setMethod( 
  f="initialize",
  signature="SaemixDiscreteCovariate",
  definition=function(.Object, name="", type="binary", reference, groups){
    .Object <- callNextMethod(.Object, name, type)
    if(!missing(reference)) .Object@reference <- reference else .Object@reference <- ""
    if(!missing(groups)) .Object@groups <- groups else .Object@groups <- list()
    validObject(.Object)
    return(.Object)
  }
)

# Getteur
setMethod(
  f ="[",
  signature = "SaemixDiscreteCovariate" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "groups"={return(x@groups)},
            "reference"={return(x@reference)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixDiscreteCovariate" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "groups"={x@groups<-value},
            "reference"={x@reference<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

########################################################################

#' Function to create a SaemixCovariate object
#' 
#' This function creates a continuous, categorical or binary covariate structure
#' to be passed on to a SaemixData object.
#' 
#' @name saemixCov
#' 
#' @param name name of the covariate (defaults to empty (""))
#' @param type possible types are continuous, categorical or binary. Defaults to continuous
#' @param unit unit of the covariate (defaults to empty (""))
#'  
#' @examples
#' age<-saemixCov(name="age", unit="yr")
#' age
#' age<-saemixCov("age", unit="yr") # same result
#' age
#' comorb1<-saemixCov(name="diabetes", type="binary")
#' comorb1
#' pgp<-saemixCov(name="PgP", type="categorical")
#' pgp
#' @export 

saemixCov<-function(name="",type="continuous", unit="") {
  new(Class="SaemixCovariate", name=name, type=type, unit=unit)
}


continousCov<-function(name="",type="continuous", unit="") {
  new(Class="SaemixCovariate", name=name, type=type, unit=unit)
}

# Associated testthat test
# 
# context("Creating SaemixCovariate objects")
# 
# test_that("Creating covariates", {
#   cov1<-saemixCov()
#   cov2<-saemixCov(name="cov2")
#   age<-saemixCov("age", unit="yr")
#   comorb1<-saemixCov(name="diabetes", type="binary")
#   expect_equal(cov1@name, "")
#   expect_equal(cov2@name, "cov2")
#   expect_equal(age@type, "continuous")
#   expect_equal(comorb1@type, "binary")
# })
# 

########################################################################
# Print and show methods

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


setMethod("show","SaemixDiscreteCovariate",
          function(object) {
            cat(object@type,"covariate")
            if(object@name!="") cat(": ",object@name)
            if(object@unit!="") cat(paste0(" (unit: ",object@unit,")"))
            cat("\n")
          }
)

setMethod("showall","SaemixDiscreteCovariate",
          function(object) {
            cat(object@type,"covariate")
            if(object@name!="") cat(": ",object@name)
            if(object@unit!="") cat(paste0(" (unit: ",object@unit,")"))
            cat("\n")
            if(length(object@groups)>0) {
              cat("     groups: ")
              for(i in 1:length(object@groups)) cat(paste0(i,"=(", paste(object@groups[[i]],collapse=","),")  "))
              cat("\n")
            }
            cat("     reference class:", object@reference,"\n")
          }
)


setMethod("print","SaemixDiscreteCovariate",
          function(x,nlines=10,...) {
            show(x)
          }
)

setMethod("show","SaemixContinuousCovariate",
          function(object) {
            cat(object@type,"covariate")
            if(object@name!="") cat(": ",object@name)
            if(object@unit!="") cat(paste0(" (unit: ",object@unit,")"))
            cat("\n")
            if(!identical(object@transform.function, function(x) x)) {
              cat("    transformation applied\n")
            }
            if(!identical(object@centering.function, function(x) 1)) {
              cat("    centering function applied\n")
            } else {
              if(length(object@centering.value)>0) cat("    centering on value:",object@centering.value,"\n")
            }
          }
)

setMethod("showall","SaemixContinuousCovariate",
          function(object) {
            cat(object@type,"covariate")
            if(object@name!="") cat(": ",object@name)
            if(object@unit!="") cat(paste0(" (unit: ",object@unit,")"))
            cat("\n")
            if(!identical(object@transform.function, function(x) x)) {
              cat("    transformation:\n")
              print(object@transform.function)
            }
            if(!identical(object@centering.function, function(x) 1)) {
              cat("    centering function:\n")
              print(object@centering.function)
            } else {
              if(length(object@centering.value)>0) cat("    centering on value:",object@centering.value,"\n")
            }
          }
)

setMethod("print","SaemixContinuousCovariate",
          function(x,nlines=10,...) {
            show(x)
          }
)


########################################################################
#' Transforming covariates
#' 
#' This method allows to apply the transformation contained in a SaemixContinuousCovariate or SaemixDiscreteCovariate
#' object to a vector of covariate values in order to include the covariate in a saemix model.
#' 
#' @name transformCovariate-methods
#' @aliases transformCovariate 
#' @aliases transformCovariate,SaemixContinuousCovariate transformCovariate,SaemixContinuousCovariate-method
#' @aliases transformCovariate,SaemixDiscreteCovariate transformCovariate,SaemixDiscreteCovariate-method
#' 
#' @param object an object of class SaemixDiscreteCovariate or SaemixContinuousCovariate
#' @param x the values of the covariate to transform
#'     
#' @return a vector or dataframe with the transformed values
#' 
#' @details 
#' For continuous covariates, the transformation contained in object is applied to the values in x.
#' For binary covariates, x is transformed to 0 for the values corresponding to the reference value and 1 otherwise.
#' For categorical covariates with ncat (=3 or more) categories, a dataframe is created with (ncat-1) dummy variables in columns. Each column is 1 if the corresponding value of x is in the corresponding group and 0 otherwise.
#' 
#' @examples 
#' # Transforming a vector of weight to log(weight/mean(weight))
#' weightCov<-new(Class="SaemixContinuousCovariate", name="Weight", transform.function=log, centering.function=median)
#' print(c(60, 70, 80))
#' transformCovariate(weightCov, c(60, 70, 80))
#' 
#' # Transforming a gender covariate given as Female/Male to 0/1 with Female as reference (0)
#' sexCov<-new(Class="SaemixDiscreteCovariate", name="gender", reference="Female")
#' transformCovariate(sexCov, c("Female","Male","Male"))
#' # Also works with factors
#' transformCovariate(sexCov, as.factor(c("Female","Male","Male")))
#' 
#' # Regrouping a covariate with 5 categories in 3 categories
#' scoreCov<-new(Class="SaemixDiscreteCovariate", name="score", type="categorical",groups=list(c(1,2),c(3,4),c(5)))
#' transformCovariate(scoreCov, c(1,2,3,4,5))
#' 
#' @exportMethod transformCovariate


setMethod(
  f ="transformCovariate",
  signature = "SaemixContinuousCovariate" ,
  definition = function (covariateModel, x=NULL){
    if(is.null(x)) return(NULL)
    if(length(covariateModel@centering.value)==0) 
      covariateModel@transform.function(x/covariateModel@centering.function(x)) else 
        covariateModel@transform.function(x/covariateModel@centering.value)
  }
)

setMethod(
  f ="transformCovariate",
  signature = "SaemixDiscreteCovariate" ,
  definition = function (covariateModel, x=NULL){
    if(is.null(x)) return(NULL)
    if(covariateModel@reference=="") reference<-sort(unique(as.character(x)))[1] else reference<-as.character(covariateModel@reference)
    if(covariateModel@type=="binary" & length(unique(x))<=2) {
      return(1-as.integer(as.character(x)==reference))
    } else { # type="categorical" OR x has more than 2 unique values
      ncat<-length(unique(x))
      if(ncat>length(x)/2) {
        message("x may not be a categorical covariate, it seems to have too many different values")
      }
      # Create (ncat-1) categories
      if(is.factor(x)) xcat<-levels(x) else xcat<-sort(unique(as.character(x)))
      if(length(covariateModel@groups)==0) {
        groups<-as.vector(xcat, mode="list")
        names(groups)<-xcat
      } else groups<-covariateModel@groups
      if(is.null(names(groups))) {
        l1<-c()
        for(i in 1:length(groups)) l1<-c(l1,as.character(groups[[i]][1]))
        names(groups)<-l1
      }
      idx1<-grep(reference,groups)
      if(length(idx1)==0) {
        message(paste("Can't find reference category ", reference,", using first category instead"))
        reference<-groups[[1]][1]
        idx1<-1
      }
      ncat<-length(groups)
      tabcat<-NULL
      for(i in 1:ncat) {
        if(i!=idx1) {
          icat<-as.integer(as.character(x) %in% groups[[i]])
          tabcat<-cbind(tabcat, icat)
        }
      }
      colnames(tabcat)<-paste0(covariateModel@name,".",names(groups)[-idx1])
      if(dim(tabcat)[2]==1) tabcat<-c(tabcat) # regrouped to binary covariate
      return(tabcat)
    }
  }
)

########################################################################
