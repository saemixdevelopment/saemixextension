
#' Class "SaemixParameter"
#' 
#' An object of the SaemixParameter class, representing an observation model for an outcome. Outcomes can be discrete or survival-type data
#' (SaemixDiscreteOutcome) or continuous (SaemixContinuousOutcome), with both subclasses inheriting from the SaemixParameter class
#' 
#' @name SaemixParameter-class 
#' @docType class
#' @aliases SaemixParameter SaemixParameter-class 
#' @aliases print,SaemixParameter showall,SaemixParameter show,SaemixParameter
#' @aliases SaemixDiscreteOutcome SaemixDiscreteOutcome-class 
#' @aliases print,SaemixDiscreteOutcome showall,SaemixDiscreteOutcome show,SaemixDiscreteOutcome
#' @aliases SaemixContinuousOutcome SaemixContinuousOutcome-class 
#' @aliases print,SaemixContinuousOutcome showall,SaemixContinuousOutcome show,SaemixContinuousOutcome
#' 
#' @section Objects from the Class: 
#' An object of the SaemixParameter class contains the following slots:
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
#' showClass("SaemixParameter")
#' 
#' @keywords classes
#' @exportClass SaemixParameter


# Variability level class - generic
setClass(Class = "SaemixParameter",
         representation=representation(
           name.par = "character", # parameter name
           distribution = "character", # distribution (currently one of normal, log-normal, logit, probit)
           estim = "character", # one of estimated, fixed or prior (currently only estimated and fixed are valid)
           initial = "numeric", # initial value (psi0)
           varlevel = "list", # variability levels
           covariate = "character" # a list of covariates
         ),
         validity=function(object){
           if(!(object@distribution %in% c("normal","lognormal","logit","probit"))) {
             message(paste("SaemixParameter: distribution",object@distribution,"cannot be handled, please use one of normal, lognormal, logit, or probit \n"))
             return("Creation of SaemixParameter failed")
           }
           return(TRUE)
         }
)


setMethod( 
  f="initialize",
  signature="SaemixParameter",
  definition=function(.Object, name.par="theta", distribution="normal", estim='estimated', initial=1, varlevel="id", covariate=""){
    .Object@name.par <- name.par
    .Object@distribution <- distribution
    .Object@estim <- estim
    .Object@initial <- initial
    if(!is.list(varlevel)) varlevel<-list(varlevel)
    .Object@varlevel <- varlevel
    .Object@covariate <- covariate
    validObject(.Object)
    return(.Object)
  }
)


# Getteur
setMethod(
  f ="[",
  signature = "SaemixParameter" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.par"={return(x@name.par)},
            "distribution"={return(x@distribution)},
            "estim"={return(x@estim)},
            "initial"={return(x@initial)},
            "varlevel"={return(x@varlevel)},
            "covariates"={return(x@covariates)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixParameter" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.par"={x@name.par<-value},
            "distribution"={x@distribution<-value},
            "estim"={x@estim<-value},
            "initial"={x@initial<-value},
            "varlevel"={x@varlevel<-value},
            "covariates"={x@covariates<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)
