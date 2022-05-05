
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
#' @slot name Object of class \code{"character"}: parameter name in the model (psi)
#' @slot distribution distribution (currently one of normal, lognormal, logit, probit)
#' @slot transform associated function, transforming a linear form of fixed effects, covariates and random effects into the parameter transformation h where psi=h(phi)
#' @slot inversetransform inverse transformation h-1 (phi=h-1(psi))
#' @slot mu.fix "estimated" if the parameter is estimated, "fixed" otherwise (defaults to estimated)
#' @slot mu typical population value for the parameter
#' @slot omega.level associated variability level (currently only IIV (associated to the subject/group) is possible)
#' @slot omega.fixed 1 if the SD is fixed to its starting value
#' @slot omega the value of the SD
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
           name = "character", # parameter name
           distribution = "character", # distribution (currently one of normal, lognormal, logit, probit)
           transform = "function", # transformation h where psi=h(phi) [transphi]
           inversetransform = "function", # inverse transformation h-1 (phi=h-1(psi)) [transpsi]
           mu = "numeric", # initial value (psi0)
           mu.fix = "numeric", # 1 if mu is fixed, default is 0 (estimated)
           prior = "logical", # ignored and set to FALSE for the moment; can be extended to a SaemixPrior object (to be defined, containing the distribution and the parameters of the prior distribution, could have a default to eg identity)
           omega.level = "character", # vector of variables associated with the variability levels
           omega.model = "numeric", # for each variability level, 1 if the variance is present in the model, 0 if no IIV
           omega.fix = "numeric", # for each variability level, 1 if the variance is fixed (a vector if omega.level has more than one level)
           omega = "numeric", # for each variability level, initial value for SD
           rho.param = "list", # for each level of variability, a list with which parameters is the parameter correlated
           rho = "list", # a list, for each level of variability, starting values for correlation terms
           rho.fix = "list", # a list, for each level of variability, which correlation parameters are fixed (same order)
           covariate = "list" # for each variability level, a list of covariate models of type SaemixCovariate
           # beta.fix = "character", # one of estimated, fixed
           # beta = "numeric" # initial value for the covariate effects
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
  definition=function(.Object, name="theta", distribution="lognormal", estimated=TRUE, prior=FALSE, mu.start=1, omega.start=1, omega.fix=FALSE, omega.level="id", covariate=list(), rho.param=list(), rho=list(), rho.fix=list()){
    .Object@name <- name
    if(tolower(distribution) %in% c("normal","norm","n"))  distribution<-"normal"
    if(tolower(distribution) %in% c("lognormal","ln"))  distribution<-"lognormal"
    if(!(distribution %in% c("normal","lognormal","logit","probit"))) {
      message(paste("SaemixParameter: distribution",object@distribution,"cannot be handled, please use one of normal, lognormal, logit, or probit \n"))
      return("Creation of SaemixParameter failed")
    }
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
    .Object@distribution <- distribution
    #    if(!(estim %in% c("estimated","fixed","prior"))) estim<-"estimated"
    .Object@mu.fix <- 1-as.integer(estimated)
    .Object@mu <- mu.start
    # Choose: either omega.level is a SaemixVarLevel object or it's just the list of variables defining variability levels
#    if(!is.list(omega.level)) omega.level<-list(omega.level)
    if(length(omega.level)>0) {
      if(length(omega.level)>1) message("Currently only one level of variability is handled by saemix, ignoring nested variability levels")
      .Object@omega.level <- omega.level[1]
      .Object@omega.model <- 1
      .Object@omega <- omega.start[1]
      .Object@omega.fix <- as.integer(omega.fix[1])
      if(length(rho.param)>0) {
        if(is(rho.param,"character")) rho.param<-list(rho.param) # given as c() instead of list() (valid for one level of variability)
        if(length(rho)>0 && is(rho,"numeric")) rho<-list(rho)
        if(is(rho.fix,"numeric")) rho.fix<-list(rho.fix)
        if(length(rho.param)>1) rho.param<-rho.param[[1]]
        if(length(rho)>1) rho<-rho[[1]]
        if(length(rho.fix)>1) rho.fix<-rho.fix[[1]]
        .Object@rho.param<-rho.param
        if(length(rho)==0  | length(unlist(rho))!=length(unlist(rho.param)) | !is(unlist(rho),"numeric")) { # starting values for correlations not given or wrong size
          rho<-vector(mode="list", length=length(rho.param))
          for(i in 1:length(rho)) rho[[i]]<-rep(0.5, length(rho.param[[i]])) # initialise all correlations to 0.5
        }
        .Object@rho<-rho
        if(length(rho.fix)==0 & length(unlist(rho.fix))!=length(unlist(rho.param))) { # fixed correlations
          rho.fix<-vector(mode="list", length=length(rho.param))
          for(i in 1:length(rho.fix)) rho.fix[[i]]<-rep(0, length(rho.param[[i]])) # by default, all correlations to be estimated
        }
        .Object@rho.fix<-rho.fix
      }
    } else {
      .Object@omega.model <- 0
    }
    if(length(covariate)>0) {
      ipb<-0
      covariate<-as.list(covariate)
      for(i in 1:length(covariate)) {
        if(!is(covariate[[i]],"SaemixCovariate")) ipb<-ipb+1 else covariate[[i]]@name<-names(covariate)[i]
      }
      if(ipb==0) .Object@covariate <- covariate else  
        message("Element covariate should be a list of SaemixCovariate elements created through the contCov(), binCov() or CatCov() function, wrong format ignored")
    }

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
            "name"={return(x@name)},
            "distribution"={return(x@distribution)},
            "transform"={return(x@transform)},
            "inversetransform"={return(x@inversetransform)},
            "mu"={return(x@mu)},
            "mu.fix"={return(x@mu.fix)},
            "prior"={return(x@prior)},
            "omega"={return(x@omega)},
            "omega.model"={return(x@omega.model)},
            "omega.fix"={return(x@omega.fix)},
            "omega.level"={return(x@omega.level)},
            "rho.param"={return(x@rho.param)},
            "rho"={return(x@rho)},
            "rho.fix"={return(x@rho.fix)},
            "covariate"={return(x@covariate)},
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
            "name"={x@name<-value},
            "distribution"={x@distribution<-value},
            "transform"={x@transform<-value},
            "inversetransform"={x@inversetransform<-value},
            "mu"={x@mu<-value},
            "mu.fix"={x@mu.fix<-value},
            "prior"={x@prior<-value},
            "omega"={x@omega<-value},
            "omega.model"={x@omega.model<-value},
            "omega.fix"={x@omega.fix<-value},
            "omega.level"={x@omega.level<-value},
            "rho.param"={x@rho.param<-value},
            "rho"={x@rho<-value},
            "rho.fix"={x@rho.fix<-value},
            "covariate"={x@covariate<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

########################################################################
# Print and show methods

setMethod("show","SaemixParameter",
          function(object) {
            cat("Model parameter", object@name,"\n")
            cat("     distribution:",object@distribution,"\n")
            cat(paste0("     mu.",object@name,":"),object@mu)
            if(object@mu.fix==1) cat(" (fixed)")
            cat("\n")
            if(length(object@omega.level)>0) {
              cat(paste0("     var.",object@name,":"),object@omega[1])
              if(object@omega.fix[1]==1) cat(" (fixed)")
              cat("\n")
              if(length(object@rho)>0 && length(object@rho[[1]])) {
                cat("     covariances:")
                wrho<-c()
                for(i in 1:length(object@rho[[1]]))
                  wrho<-c(wrho,paste0(object@rho.param[[1]][i]," (",object@rho[[1]][i],ifelse(object@rho.fix[[1]][i]==1," fixed",""),")"))
                cat(paste(object@rho,collapse=","),"\n")
              }
            }
            if(length(object@covariate)>0) {
              cat("     covariates:",names(object@covariate),"\n")
            }
          }
)


setMethod("print","SaemixParameter",
          function(x,nlines=10,...) {
            show(x)
          }
)

setMethod("showall","SaemixParameter",
          function(object) {
            show(object)
            cat("Transformation h(phi):\n")
            print(object@transform)
            cat("Inverse transformation h-1(psi):\n")
            print(object@inversetransform)
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


