#' @include aaa_generics.R
#' @include SaemixData.R
#' @include SaemixData-methods.R
#' @include SaemixData-methods_covariates.R
#' @include SaemixParameter.R
#' @include SaemixVarLevel.R
NULL

setClass(Class = "SaemixIndivModel",
         representation=representation(
           log="character",		# A record of the warnings and messages during the creation of the object
           nphi = "numeric", # number of model parameters (size of param.names)
           param.names="character", # names of the model parameters
           distribution = "character", # a vector specifiying the distribution of each parameter (currently one of normal, lognormal, logit, probit)
           transform="list", # list of functions to transform parameters
           invtransform = "list", # list of inverse transformation for parameters
           dtransform="list", # list of derivatives of transform functions
           varlevel = "character", # variability levels (vector of grouping levels)
           covariate = "character", # a vector giving the names of the covariates in the model
           popmodel = "list", # a list with for each level in varlevel, the fixed effect model as a SaemixPopModelHat object
           varmodel = "list" # a list with for each level in varlevel, the variability model as a SaemixVarModelHat object
         ),
         validity=function(object){
           # Check sizes, check types of lists
           return(TRUE)
         }
)

# More structured definition of an individual model matching covariate model and data for a given variability level
# defining design variables and setting indices to access items
# in the end, should be created by a function associating a variability model, a covariate model (par-cov relationship and cov transformation) and data (find cov, find level of variability)
# Contains most of what used to be Uargs


setMethod( 
  f="initialize",
  signature="SaemixIndivModel",
  definition=function(.Object, parameters, verbose=FALSE){
    if(missing(parameters)) {
      if(verbose) cat("Please supply a list of parameters\n")
      return("Creation of SaemixIndivModel object failed \n")
    }
    for(i in parameters) {
      if(!inherits(i,"SaemixParameter")) {
        if(verbose) cat("Please supply a list of parameters created via the SaemixParam() function\n")
        return("Creation of SaemixIndivModel object failed \n")
      }
    }
    xcheck<-checkParameters(parameters)
    logmsg<-xcheck$logmsg
    parameters<-xcheck$parameters
    .Object@varlevel <- xcheck$varlevel.order
    nphi<-length(parameters)
    .Object@nphi <- nphi
    .Object@param.names <- names(parameters)
    distribution<-c()
    transform<-invtransform<-vector(mode="list",length=nphi)
    for(ipar in 1:nphi) {
      distribution<-c(distribution, parameters[[ipar]]@distribution)
      transform[[ipar]]<-parameters[[ipar]]@transform
      invtransform[[ipar]]<-parameters[[ipar]]@invtransform
    }
    covariate<-c()
    for(ipar in 1:nphi) 
      if(length(parameters[[ipar]]@covariate)>0) covariate<-c(covariate, parameters[[ipar]]@covariate)
    if(length(covariate)>0) covariate<-unique(covariate) else covariate<-c()
    .Object@covariate <- covariate
    .Object@distribution <- distribution
    .Object@transform <- transform
    .Object@invtransform <- invtransform
    .Object@varmodel <- extractVarModel(parameters)
    .Object@popmodel <- extractFixedEffectModel(parameters)
    .Object@log <- logmsg
    
    validObject(.Object)
    return(.Object)
  }
)


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
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)


################################################################################
# Print/show functions for SaemixVarModel  (statistical model)

setMethod("print","SaemixIndivModel",
          function(x,nlines=10,...) {
            cat("Individual model object\n")
            cat("      parameters",x@param.names,"\n      variability levels",x@varlevel,"\n")
            if(length(x@covariate)>0) cat("      covariates",x@covariate,"\n")
          }
)

setMethod("show","SaemixIndivModel",
          function(object) {
            cat("Individual model object\n")
            cat("      parameters",object@param.names,"\n      variability levels",object@varlevel,"\n")
            if(length(object@covariate)>0) cat("      covariates",object@covariate,"\n")
          }
)
setMethod("showall","SaemixIndivModel",
          function(object) {
            cat("Individual model object\n")
            cat("      parameters",object@param.names,"\n      variability levels",object@varlevel,"\n")
            if(length(object@covariate)>0) cat("      covariates",object@covariate,"\n")
            cat("Fixed effect model(s)\n")
            showall(object@popmodel)
            cat("Variability model(s)\n")
            showall(object@varmodel)
          }
)

################################################################################
# Computational functions

# Replaces 
## transphi (with param=phi and transform=model@transform)
## transpsi (with param=psi and transform=model@invtransform)
## dtransphi (with param=phi and transform=model@dtransform) [check format and maybe create different function]

transformPar <- function(param, transform) {
  # param needs to be a matrix (not a vector)
  tparam<-param
  # check if this is needed (then change phi to param)
  # if(is.null(dim(phi))) {
  #    dpsi<-as.matrix(t(rep(1,length(phi))))
  #    psi<-as.matrix(t(phi),nrow=1)
  # } else 
  #   dpsi<-matrix(1,dim(phi)[1],dim(phi)[2])
  for(i in 1:dim(param)[2]) {
    tparam[,i]<-transform[[i]](param[,i])
  }
#  if(is.null(dim(phi))) dpsi<-c(dpsi)
  return(tparam)
}


