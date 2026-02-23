#' @include aaa_generics.R
NULL

#' Class "SaemixVarPhi"
#' 
#' An object of the SaemixVarPhi class, containing the individual parameters for a given level of variability, on the eta scale
#' 
#' @name SaemixVarPhi-class 
#' @docType class
#' @aliases SaemixVarPhi SaemixVarPhi-class 
#' @aliases print,SaemixVarPhi showall,SaemixVarPhi show,SaemixVarPhi
#' 
#' @section Objects from the Class: 
#' An object of the SaemixVarPhi class contains the following slots:
#' @slot name name of the grouping level (needed ? TBD)
#' @slot log A string recording the warnings and messages during the creation of the object
#' @slot phi individual parameters (mu+beta.Cov + eta)
#' @slot eta random effects associated to this level of variability (TBD, maybe here only phi)
#' @slot Mcovar design matrix of the covariates associated to this level of variability
#' @slot mean.phi mean (population) individual parameters (mu+beta.Cov)
#' @slot eta.samp etas sampled at this level of variability (TBD)
#' @slot phi.samp etas sampled at this level of variability (TBD)
#' @slot index indicator to match the lines in the dataframe to the corresponding lines in the full dataset (ex: for IIV, this would be 1 for all lines corresponding to subject 1, 2 for all lines of subject 2, etc...; for IOV, 1.1 would be occasion 1 in subject 1, 1.2 occasion 2 in subject 1, 2.1 occasion 1 in subject 2, etc...)
#' @section Methods:
#'   \describe{
#'     \item{[<-}{\code{signature(x = "SaemixVarPhi")}: replace elements of object}
#'     \item{[}{\code{signature(x = "SaemixVarPhi")}: access elements of object}
#'     \item{initialize}{\code{signature(.Object = "SaemixVarPhi")}: internal function to initialise object, not to be used}
#'     \item{print}{\code{signature(x = "SaemixVarPhi")}: prints details about the object (more extensive than show)}
#'     \item{showall}{\code{signature(object = "SaemixVarPhi")}: shows all the elements in the object}
#'     \item{show}{\code{signature(object = "SaemixVarPhi")}: prints details about the object}
#' 	 }
#' 
#' @details
#' Intermediate object holding the levels 
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' 
#' @examples
#' showClass("SaemixVarPhi")
#' 
#' @keywords classes
#' @exportClass SaemixVarPhi

# Variability level class - generic
setClass(Class = "SaemixVarPhi",
         representation=representation(
           name = "character", # name of the variability level
           log="character",		# A record of the warnings and messages during the creation of the object
           phi = "matrix", # individual parameters
           eta = "matrix", # individual random effects
           Mcovar = "matrix", # matrix of covariates
           mean.phi = "matrix", # mean individual parameter
           eta.samp = "matrix", # sampled random effects
           phi.samp = "matrix", # sampled individual parameters
           index = "numeric" # index to match to the dataset
         ),
         validity=function(object){
           return(TRUE)
         }
)

# Getteur
setMethod(
  f ="[",
  signature = "SaemixVarPhi" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name"={return(x@name)},
            "log"={return(x@log)},
            "phi"={return(x@phi)},
            "eta"={return(x@eta)},
            "Mcovar"={return(x@Mcovar)},
            "mean.phi"={return(x@mean.phi)},
            "eta.samp"={return(x@eta.samp)},
            "index"={return(x@index)},
            "phi.samp"={return(x@phi.samp)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixVarPhi" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name"={x@name<-value},
            "log"={x@log<-value},
            "phi"={x@phi<-value},
            "eta"={x@eta<-value},
            "Mcovar"={x@Mcovar<-value},
            "mean.phi"={x@mean.phi<-value},
            "eta.samp"={x@eta.samp<-value},
            "index"={x@index<-value},
            "phi.samp"={x@phi.samp<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

# Initialising ab initio
setMethod( 
  f="initialize",
  signature="SaemixVarPhi",
  definition=function(.Object, name="iiv", phi, index=c(), Mcovar=matrix(ncol=0,nrow=0)) {
    .Object@name <- name
    .Object@phi <- phi
    if(length(index)>0) .Object@index <- index
    .Object@Mcovar <- Mcovar
    validObject(.Object)
    return(.Object)
  }
)

