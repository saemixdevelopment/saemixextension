#' @include aaa_generics.R
NULL

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
#' @slot name character string giving the name of the parameter (defaults to theta)
#' @slot log A string recording the warnings and messages during the creation of the object
#' @slot distribution character string, with the name of the distribution associated with the parameter; currently one of normal, lognormal, logit, probit (defaults to lognormal)
#' @slot varlevel vector giving the grouping levels with variability associated to the parameter (defaults to "iiv" usually indicating subject-level variability)
#' @slot corr a list of vector giving for each level in varlevel which parameters are correlated
#' @slot transform the name of the transformation function which transforms the random effects (which have a normal distribution)
#' to the distribution of the parameter. For the built-in distributions, the transform function is respectively the identity, 
#' logarithmic, logit, and probit transformations.
#' @slot invtransform the name of the inverse transformation function which transforms the parameters back to the scale of the random effects.
#' For the built-in distributions, the inverse transform function is respectively the identity, exponential, logistic and inverse cdf.
#' @slot mu.init numeric, initial estimate for the population value (mu) of the parameter (defaults to 1)
#' @slot mu.estim character string, one of estimated (default), fixed or prior (currently only estimated and fixed are valid)
#' @slot sd.init numeric, value of the initial estimate of the standard deviation of the parameter for each level of variability (defaults to 1 for each level of variability)
#' @slot sd.estim  vector of character strings, with values estimated (default), fixed or prior (currently only estimated and fixed are valid)
#' @slot corr.init a list of vectors of numeric values between -1 and 1, giving the initial estimate for the correlations associated with the parameter for each level of variability (defaults to 0)
#' @slot corr.estim list of vector of character strings, with values estimated (default), fixed or prior (currently only estimated and fixed are valid) to specify which correlation terms are estimated
#' @slot covariate vector giving the covariate model to be associated with the parameter as a list of covariate names
#' @slot ncovariate number of covariates associated to the parameter (length of covariate vector)
#' @slot covariate.estim  vector of character strings, with values estimated (default), fixed or prior (currently only estimated and fixed are valid)
#' @slot covariate.init numeric, vector of the initial values for the covariate effects (beta) associated to the parameter (default values are 0)
#' @slot covariate.varlevel Vector of strings of the same size as covariate giving the associated covariate level. If given, the names will be matched to the names given in varlevel. Will default to the first level of variability if not specified.
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
#' @details
#' Currently four built-in distributions are available for saemix parameters. Giving the name of the distribution 
#' (normal, lognormal, logit, probit) will automatically fill in the transform and invtransform slots of the object, 
#' overriding other functions which may be given by the user.
#' In future versions, 
#' 
#' Covariates names given in the covariate argument will be matched to the variables in the dataset when creating the model.
#' The covariate model can also be changed on the fly later when running the model. By default initial values for the covariate
#' effects are set to 0, but can be set using the covariate.init argument.
#' 
#' Correlation terms will be collated from all the parameters given in the model and only need to be specified once. 
#' The variance-covariance matrix will then be checked for consistency, and non-block diagonal matrices will be completed by adding extra covariance terms.
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' 
#' @seealso \code{\link{saemixData}} \code{\link{SaemixModel}} \code{\link{saemixControl}} \code{\link{saemix}}
#' @examples
#' showClass("SaemixParameter")
#' 
#' # default definition (log-normal parameter with initial value 1)
#' ka <- saemixParam()
#' # specifying initial values and covariate influencing the parameter
#' cl<-saemixParam(name="cl",mu.init=2, covariate=c("wt","sex","age"), covariate.init=c(0.75), covariate.mu.estim=c("fixed"))
#' 
#' @keywords classes
#' @exportClass SaemixParameter

# Variability level class - generic
setClass(Class = "SaemixParameter",
         representation=representation(
           name = "character", # parameter name
           log="character",		# A record of the warnings and messages during the creation of the object
           distribution = "character", # distribution (currently one of normal, lognormal, logit, probit)
           transform="function", # function to transform from the normal scale (h) to the scale of the model parameters, currently one of id, exp, inverse logit, inverse probit
           invtransform = "function", # inverse transformation (h-1), currently one of id, log, logit, probit
           dtransform="function", # derivative of transform, automatically set for the four default functions
           mu.estim = "character", # one of estimated, fixed or prior (currently only estimated and fixed are valid)
           mu.init = "numeric", # initial value (psi0)
           varlevel = "character", # variability levels (vector of grouping levels)
           sd.init = "numeric", # initial value of variability for each level of variability 
           sd.estim = "character", # a vector with for each parameter in varlevel, one of estimated, fixed or prior (currently only estimated and fixed are valid)
           corr = "list", # a list of vectors of strings, with for each level of variability, which parameters is the parameter correlated to
           corr.init = "list", # a list of vectors of numerical values between -1 and 1, giving the initial value of correlations in the model (same order as corr)
           corr.estim = "list", # a list of vectors of numerical values between -1 and 1, giving the initial value of correlations in the model (same order as corr)
           ncovariate = "numeric", # number of covariates in the model
           covariate = "character", # the vector of covariate names entering the model
           covariate.estim = "character", # vector of estimated, fixed or prior (currently only estimated and fixed are valid)
           covariate.init = 'numeric', # initial value for covariate effect (one per covariate even if later a covariate happens to be categorical)
           covariate.varlevel = "character" # level of variability to which each covariate is associated
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
  definition=function(.Object, name="theta", distribution="lognormal", transform="",invtransform="", dtransform="", mu.estim='estimated', mu.init=1, varlevel="iiv", covariate=c(), covariate.estim=c(),covariate.init=c(), covariate.varlevel=c(), verbose=FALSE){
    .Object@name <- name
    distribution<-tolower(distribution)
    .Object@log<-""
    if(distribution %in% c("log","log-normal")) {
      distribution<-"lognormal"
      .Object@log<-"Recognising the distribution as log-normal\n"
    }
    .Object@distribution <- distribution
    if(!(mu.estim %in% c("estimated", "fixed","prior"))) {
      if(verbose) message("The mu.estim argument must be one of estimated, fixed or prior\n")
      msg<-"The mu.estim argument must be one of estimated, fixed or prior, setting to estimated\n"
      .Object@log<-paste0(.Object@log,msg)
      mu.estim<-"estimated"
    }
    if(mu.estim=="prior") {
      msg<-"Prior function currently unavailable, setting parameter to estimated\n"
      .Object@log<-paste0(.Object@log,msg)
      mu.estim<-"estimated"
    }
    if(length(mu.init)>length(mu.estim)) {
      mu.estim<-c(mu.estim,rep("estimated",length(mu.init)-length(mu.estim)))
    }
    .Object@mu.estim <- mu.estim
    .Object@mu.init <- mu.init
#    if(!is.list(varlevel)) varlevel<-list(varlevel)           
    if(!(distribution %in% c("normal","lognormal","logit","probit"))) {
      if(transform=="" | invtransform=="" | dtransform=="") {
        if(verbose) message(paste("Please specify the transform, inverse transform function and derivative of transform function for distribution",distribution,"\n"))
      return("Creation of SaemixParameter failed")
      }
    }
    if(transform=="") {
      msg<-paste0("Automatically adding the transform and inverse transform functions for known distribution: ",distribution,"\n")
      if(verbose) cat(msg)
      if(invtransform!="") .Object@log<-paste0(.Object@log,msg)
      invtransform=switch(distribution,
                       normal=function(x) x,
                       lognormal=log,
                       logit=function(x) log(x/(1-x)),
                       probit=norminv)
      transform=switch(distribution,
                          normal=function(x) x,
                          lognormal=exp,
                          logit=function(x) 1/(1+exp(-x)),
                          probit=normcdf)
      dtransform=switch(distribution,
                       normal=function(x) 1,
                       lognormal=exp,
                       logit=function(x) 1/(2+exp(-x)+exp(x)),
                       probit=1/dnorm(qnorm(x)) )
    }
    .Object@transform <- transform
    .Object@invtransform <- invtransform
    .Object@dtransform <- dtransform
    if(length(varlevel)==0) varlevel<-"iiv"
    .Object@varlevel <- as.character(varlevel)
    if(length(covariate)>0) {
      beta.init<-rep(0,length(covariate))
      if(length(covariate.init)>0) beta.init[1:length(covariate.init)]<-covariate.init
      beta.estim<-rep("estimated",length(covariate))
      if(length(covariate.estim)>0) {
        covariate.estim[!(covariate.estim %in% c("estimated", "fixed","prior"))]<-"estimated"
        beta.estim[1:length(covariate.estim)]<-covariate.estim
      }
      if(length(covariate.varlevel)!=length(covariate) & length(covariate.varlevel)>1) {
        msg<-"   If given, covariate.varlevel must be either of length 1 or of the same size as the name.covariates vector, defaulting to the first grouping factor.\n"
        .Object@log<-paste0(.Object@log,msg)
        if(verbose) cat(msg)
        covariate.varlevel<-character()
      }
      if(length(covariate.varlevel)==1) covariate.varlevel<-rep(covariate.varlevel,length(covariate))
      if(length(covariate.varlevel)==0) covariate.varlevel<-rep(varlevel[1],length(covariate))
      
    } else {
      covariate<-character()
      covariate.varlevel<-character()
      beta.init<-numeric()
      beta.estim<-character()
    }
    .Object@covariate <- covariate
    .Object@covariate.init <- beta.init
    .Object@covariate.estim <- beta.estim
    .Object@ncovariate = length(covariate)
    .Object@covariate.varlevel <- covariate.varlevel
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
            "log"={return(x@log)},
            "distribution"={return(x@distribution)},
            "transform"={return(x@transform)},
            "dtransform"={return(x@dtransform)},
            "invtransform"={return(x@invtransform)},
            "mu.estim"={return(x@mu.estim)},
            "mu.init"={return(x@mu.init)},
            "varlevel"={return(x@varlevel)},
            "corr"={return(x@corr)},
            "sd.init"={return(x@sd.init)},
            "sd.estim"={return(x@sd.estim)},
            "corr.init"={return(x@corr.init)},
            "corr.estim"={return(x@corr.estim)},
            "covariate"={return(x@covariate)},
            "ncovariate"={return(x@ncovariate)},
            "covariate.estim"={return(x@covariate.estim)},
            "covariate.init"={return(x@covariate.init)},
            "covariate.varlevel"={return(x@covariate.varlevel)},
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
            "log"={x@log<-value},
            "distribution"={x@distribution<-value},
            "transform"={x@transform<-value},
            "invtransform"={x@invtransform<-value},
            "dtransform"={x@dtransform<-value},
            "mu.estim"={x@mu.estim<-value},
            "mu.init"={x@mu.init<-value},
            "varlevel"={x@varlevel<-value},
            "corr"={x@corr<-value},
            "sd.init"={x@sd.init<-value},
            "sd.estim"={x@sd.estim<-value},
            "corr.init"={x@corr.init<-value},
            "corr.estim"={x@corr.estim<-value},
            "covariate"={x@covariate<-value},
            "ncovariate"={x@ncovariate<-value},
            "covariate.estim"={x@covariate.estim<-value},
            "covariate.init"={x@covariate.init<-value},
            "covariate.varlevel"={x@covariate.varlevel<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)


################################################################################
# Print/show functions 

setMethod("print","SaemixParameter",
          function(x,nlines=10,...) {
            nampar <- deparse(substitute(x)) # writes x instead of the parameter name
            # x1<-match.call(expand.dots=FALSE) # passes x
            # print(x1$x)
            nampar <-x@name
            distsymb <- switch(x@distribution,
                               normal="N",
                               lognormal="LN",
                               logit="logit",
                               probit="probit")
            htrans <- switch(x@distribution,
                             normal="phi",
                             lognormal="exp(phi)",
                             logit="1/(1+exp(-phi))",
                             probit="normcdf(phi)")
            cat("Parameter:",nampar,"~", paste0(distsymb,"(",",",ifelse(length(x@varlevel)>0,"Omega",0),")"),"\n")
            form1 <- paste0(htrans,"(mu.",nampar,")")
            if(x@ncovariate>0) 
              form1<-paste0(form1, paste0(" + beta.",x@covariate," ",x@covariate, collapse=""))
            if(length(x@varlevel)>0) {
              for(ilev in x@varlevel) form1<-paste0(form1," + eta.",ilev)
            }
            cat("            phi=",form1,"\n")
            cat(x@log)
          }
)


setMethod("show","SaemixParameter",
          function(object) {
            # x1<-match.call(expand.dots=FALSE) # passes the object name... maddening
            # print(x1$object)
#            nampar <- deparse(substitute(object)) # works
            nampar<-object@name
            cat("Saemix parameter with",object@distribution," distribution ")
            if(length(object@varlevel)>0) cat("associated with variability on variable(s)",object@varlevel) else cat("without variability")
            if(object@ncovariate>0) cat(", associated with covariates:",object@covariate)
            cat("\n")
            }
)

setMethod("showall","SaemixParameter",
          function(object) {
            nampar<-object@name
            distsymb <- switch(object@distribution,
                               normal="normal",
                               lognormal="lognormal",
                               logit="logit",
                               probit="probit")
            htrans <- switch(object@distribution,
                             normal="phi",
                             lognormal="exp(phi)",
                             logit="1/(1+exp(-phi))",
                             probit="normcdf(phi)")
            form1 <- paste0(htrans,"\n              with phi=h(mu.",nampar,")")
            if(object@ncovariate>0) 
              form1<-paste0(form1, paste0(" + beta.",object@covariate," ",object@covariate, collapse=""))
            if(length(object@varlevel)>0) {
              for(ilev in object@varlevel) form1<-paste0(form1," + eta.",ilev)
            }
            cat("Parameter",object@name,":",distsymb,"distribution\n")
            if(length(object@varlevel)>0) cat("            variability levels:",paste(object@varlevel,sep=", "),"\n") else
              cat("            no variability\n")
            if(length(object@corr)>0) {
              cat("            correlated with: ")
              for(ilev in 1:length(object@corr)) 
                cat(paste(object@corr[[ilev]],collapse=", "),paste0("(",names(object@corr)[ilev],"); "))
              cat("\n")
            }
            if(object@ncovariate>0) cat("            covariates:",paste(object@covariate,sep=", "),"\n")
            cat("            formula: ",nampar,"=",form1,"\n")
            if(length(object@varlevel)>0) {
              for(ilev in object@varlevel) 
                cat(paste0("                   eta.",ilev,"~N(0, omega.",ilev,")\n"))
            }
            cat(object@log)
          }
)


