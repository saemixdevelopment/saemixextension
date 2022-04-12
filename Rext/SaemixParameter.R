
#' Class "SaemixParameter"
#' 
#' An object of the SaemixParameter class, representing a parameter in the model
#' 
#' @name SaemixParameter-class 
#' @docType class
#' @aliases SaemixParameter SaemixParameter-class 
#' @aliases print,SaemixParameter showall,SaemixParameter show,SaemixParameter
#' 
#' @section Objects from the Class: 
#' An object of the SaemixParameter class contains the following slots:
#' @slot name.par Object of class \code{"character"}: parameter name in the model (psi)
#' @slot distribution distribution (currently one of normal, lognormal, logit, probit)
#' @slot transform associated function, transforming a linear form of fixed effects, covariates and random effects into the parameter transformation h where psi=h(phi)
#' @slot inversetransform inverse transformation h-1 (phi=h-1(psi))
#' @slot mu.fix "estimated" if the parameter is estimated, "fixed" otherwise (defaults to estimated)
#' @slot mu typical population value for the parameter
#' @slot var.level associated variability level (currently only IIV (associated to the subject/group) is possible)
#' @slot var.fixed "estimated" if the variance is estimated, "fixed" otherwise (a vector if var.level has more than one level, currently not available)
#' @slot var the value of the variance
#' @slot covariate the vector of covariate names which have an effect on the parameter
#' @slot beta.fix "estimated" if the associated covariate parameter beta is estimated, "fixed" otherwise
#' @slot beta the value of the parameter(s) beta
#' 
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


# Parameter class - generic
setClass(Class = "SaemixParameter",
         representation=representation(
           name.par = "character", # parameter name
           distribution = "character", # distribution (currently one of normal, lognormal, logit, probit)
           transform = "function", # transformation h where psi=h(phi) [transphi]
           inversetransform = "function", # inverse transformation h-1 (phi=h-1(psi)) [transpsi]
           mu.fix = "character", # one of estimated, fixed or prior (currently only estimated and fixed are valid)
           mu = "numeric", # initial value (psi0)
           var.level = "list", # variability levels
           var.fix = "character", # one of estimated, fixed (a vector if var.level has more than one level)
           var = "numeric", # initial value for variances
           covariate = "character", # a list of covariates
           beta.fix = "character", # one of estimated, fixed
           beta = "numeric" # initial value for the covariate effects
         ),
         validity=function(object){
           return(TRUE)
         }
)

# Inverse of the normal cumulative distribution fct: using erfcinv from ?pnorm
norminv<-function(x,mu=0,sigma=1)  mu-sigma*qnorm(x,lower.tail=FALSE)

# Truncated gaussian distribution (verifie par rapport a definition de erf/matlab)
normcdf<-function(x,mu=0,sigma=1)
  cutoff(pnorm(-(x-mu)/sigma,lower.tail=FALSE),1e-30)

setMethod( 
  f="initialize",
  signature="SaemixParameter",
  definition=function(.Object, name.par="theta", distribution="normal", estim='estimated', initial=1, var.level="id", covariate=""){
    .Object@name.par <- name.par
    if(!(distribution %in% c("normal","lognormal","logit","probit"))) {
      message(paste("SaemixParameter: distribution",object@distribution,"cannot be handled, please use one of normal, lognormal, logit, or probit \n"))
      return("Creation of SaemixParameter failed")
    }
    .Object@distribution <- distribution
    if(distribution=="normal") {
      .Object@transform <-function(x) x
      .Object@inversetransform <-function(x) x
    }
    if(distribution=="lognormal") {
      .Object@transform <-exp
      .Object@inversetransform <-log
    }
    if(distribution=="logit") {
      .Object@transform <- function(x) 1/(1+exp(-x))
      .Object@inversetransform <-function(x) log(x/(1-x))
    }
    if(distribution=="probit") {
      .Object@transform <-function(x) normcdf(x)
      .Object@inversetransform <-function(x) norminv(x)
    }
    if(!(estim %in% c("estimated","fixed"))) estim<-"estimated"
#    if(!(estim %in% c("estimated","fixed","prior"))) estim<-"estimated"
    .Object@estim <- estim
    .Object@initial <- initial
    if(!is.list(var.level)) var.level<-list(var.level)
    .Object@var.level <- var.level
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
            "transform"={return(x@transform)},
            "inversetransform"={return(x@inversetransform)},
            "mu.fix"={return(x@mu.fix)},
            "mu"={return(x@mu)},
            "var"={return(x@var)},
            "var.fix"={return(x@var.fix)},
            "var.level"={return(x@var.level)},
            "covariate"={return(x@covariate)},
            "beta"={return(x@beta)},
            "beta.fix"={return(x@beta.fix)},
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
            "transform"={x@transform<-value},
            "inversetransform"={x@inversetransform<-value},
            "mu"={x@mu<-value},
            "mu.fix"={x@mu.fix<-value},
            "var.level"={x@var.level<-value},
            "var.fix"={x@var.fix<-value},
            "var"={x@var<-value},
            "covariate"={x@covariate<-value},
            "beta.fix"={x@beta.fix<-value},
            "beta"={x@beta<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

########################################################################
# could be transphi

# setMethod(
#   f ="transformParameter",
#   signature = "SaemixParameter" ,
#   definition = function (parameterModel, x=NULL){
#     if(is.null(x)) return(NULL)
#     parameterModel@transform(x)
#   }
# )


