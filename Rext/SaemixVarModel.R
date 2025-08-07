#' @include aaa_generics.R
#' @include SaemixParameter.R
NULL

#' Class "SaemixVarModel"
#' 
#' An object of the SaemixVarModel class, representing the model for a variability level. 
#' The SaemixVarModelHat class inherits from the SaemixVarModel class and adds slots for the estimated parameters.
#' 
#' @name SaemixVarModel-class 
#' @docType class
#' @aliases SaemixVarModel SaemixVarModel-class 
#' @aliases print,SaemixVarModel showall,SaemixVarModel show,SaemixVarModel
#' 
#' @aliases SaemixVarModelHat SaemixVarModelHat-class 
#' @aliases print,SaemixVarModelHat showall,SaemixVarModelHat show,SaemixVarModelHat
#' 
#' @section Objects from the Class: 
#' An object of the SaemixVarModel class contains the following slots:
#' @slot name.level string giving the name of variability level (defaults to "iiv" as the subject-level grouping)
#' @slot variable string giving the corresponding name of the variable (in the dataset). Defaults to "id" as the subject-level grouping [maybe remove]
#' @slot nphi number of parameters in the model (size of omega matrix)
#' @slot omega variance-covariance matrix (initialised to the CI for the variance and covariance parameters)
#' @slot omega.model model (0/1 for parameters in the model or nor)
#' @slot omega.estim etc...
#' 
#' 
#' @section Methods:
#'   \describe{
#'     \item{[<-}{\code{signature(x = "SaemixVarModel")}: replace elements of object}
#'     \item{[}{\code{signature(x = "SaemixVarModel")}: access elements of object}
#'     \item{initialize}{\code{signature(.Object = "SaemixVarModel")}: internal function to initialise object, not to be used}
#'     \item{print}{\code{signature(x = "SaemixVarModel")}: prints details about the object (more extensive than show)}
#'     \item{showall}{\code{signature(object = "SaemixVarModel")}: shows all the elements in the object}
#'     \item{show}{\code{signature(object = "SaemixVarModel")}: prints details about the object}
#' 	 }
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' 
#' @seealso \code{\link{SaemixData}} \code{\link{SaemixModel}} \code{\link{saemixControl}} \code{\link{saemix}}
#' @examples
#' showClass("SaemixVarModel")
#' 
#' @keywords classes
#' @exportClass SaemixVarModel
#' @exportClass SaemixVarModelHat


# Variability level class - generic
setClass(Class = "SaemixVarModel",
         representation=representation(
           name.level = "character", # name of variability level
           variable = "character", # level of grouping = which variable (in the dataset) is the variability associated to  (maybe remove ? it should be created when associated to the dataset... or transfer it to SaemixVarModelHat ?)
           log = "character", # warning messages
           nphi = "numeric", # number of parameters (size of omega.model)
           param = "numeric", # vector of parameters in the model (variance and covariances)
           param.names="character", # names of the parameters to be estimated (variances and covariances)
           omega.model = "matrix", # variance-covariance matrix (square matrix) of 0/1, 1 indicating that the parameter is present in the model
           omega.estim = "matrix", # variance-covariance matrix with elements estimated, fixed or prior; 
           omega = "matrix", # variance-covariance matrix omega
           subomega ="matrix", # positive-definite submatrix of omega
           idvec.var= "numeric", # indices of variance terms in param
           idvec.cov= "numeric", # indices of covariance terms in param
           idvec.estim= "numeric", # indices of elements estimated in the matrix, as vector
           idcol.eta= "numeric", # which parameters have variability (subomega=omega[idcol.eta, idcol.eta])
           idcol.eta.fix= "numeric", # which parameters have variability and this variability is fixed
#           idcol.eta.prior= "numeric", # which parameters have variability and this variability is defined via a prior (future dev)
           idmat.var= "numeric", # indices of variance term in lower triangular matrix (with diag=TRUE)
           idmat.cov= "numeric", # indices of covariance term in lower triangular matrix (with diag=TRUE)
           idmat.estim= "numeric" # indices of elements estimated in the matrix, as vector
         ),
         validity=function(object){
           return(TRUE)
         }
)

setClass(Class = "SaemixVarModelHat",
         contains = "SaemixVarModel",
         representation=representation(
           omega.hat = "matrix", # estimated variance-covariance matrix
           omega.var = "matrix", # estimated variance of estimation for the variance-covariance matrix
           param.hat = "numeric", # estimated parameters (same order as param.names)
           param.se = "numeric", # estimated SE for the parameters (same order as param.names)
           conf.int = "data.frame" # Table giving the estimates, SE, CV and confidence intervals (assuming normality of the estimators)
         ),
         validity=function(object){
           return(TRUE)
         }
)


setMethod( 
  f="initialize",
  signature="SaemixVarModel",
  definition=function(.Object, omega.model, name.level="iiv", verbose=FALSE){
    if(missing(omega.model)) {
      if(verbose) cat("Please supply a valid omega.model matrix\n")
      return("Creation of VarModel object failed \n")
    }
    .Object@omega.model <- omega.model
    .Object@name.level <- name.level
    validObject(.Object)
    return(.Object)
  }
)

# setMethod( 
#   f="initialize",
#   signature="SaemixVarModelHat",
#   definition=function(.Object, omega.model, name.level="iiv", verbose=FALSE){
#     .Object <- callNextMethod(.Object,omega.model=omega.model, name.level=name.level, verbose=verbose)
# #    .Object@variable <- variable
#     validObject(.Object)
#     return(.Object)
#   }
# )

setMethod( 
  f="initialize",
  signature="SaemixVarModelHat",
  definition=function(.Object, saemix.varmodel, verbose=FALSE){
    .Object <- callNextMethod(.Object,omega.model=saemix.varmodel@omega.model, name.level=saemix.varmodel@name.level, verbose=verbose)
    for(i in slotNames(saemix.varmodel))
      slot(.Object,i)<-slot(saemix.varmodel,i)
    #    .Object@variable <- variable
    validObject(.Object)
    return(.Object)
  }
)



# Getteur
setMethod(
  f ="[",
  signature = "SaemixVarModel" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.level"={return(x@name.level)},
            "variable"={return(x@variable)},
            "log"={return(x@log)},
            "nphi"={return(x@nphi)},
            "param"={return(x@param)},
            "param.names"={return(x@param.names)},
            "omega"={return(x@omega)},
            "omega.estim"={return(x@omega.estim)},
            "omega.model"={return(x@omega.model)},
            "subomega"={return(x@subomega)},
            "idvec.estim"={return(x@idvec.estim)},
            "idvec.var"={return(x@idvec.var)},
            "idvec.cov"={return(x@idvec.cov)},
            "idcol.eta"={return(x@idcol.eta)},
            "idcol.eta.fix"={return(x@idcol.eta.fix)},
            "idmat.var"={return(x@idmat.var)},
            "idmat.cov"={return(x@idmat.cov)},
            "idmat.estim"={return(x@idmat.estim)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixVarModel" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.level"={x@name.level<-value},
            "variable"={x@variable<-value},
            "log"={x@log<-value},
            "nphi"={x@nphi<-value},
            "param"={x@param<-value},
            "param.names"={x@param.names<-value},
            "omega"={x@omega<-value},
            "omega.estim"={x@omega.estim<-value},
            "omega.model"={x@omega.model<-value},
            "subomega"={x@subomega<-value},
            "idvec.var"={x@idvec.var<-value},
            "idvec.cov"={x@idvec.cov<-value},
            "idvec.estim"={x@idvec.estim<-value},
            "idcol.eta"={x@idcol.eta<-value},
            "idcol.eta.fix"={x@idcol.eta.fix<-value},
            "idmat.var"={x@idmat.var<-value},
            "idmat.cov"={x@idmat.cov<-value},
            "idmat.estim"={x@idmat.estim<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)


# Getteur
setMethod(
  f ="[",
  signature = "SaemixVarModelHat" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.level"={return(x@name.level)},
            "variable"={return(x@variable)},
            "log"={return(x@log)},
            "nphi"={return(x@nphi)},
            "param"={return(x@param)},
            "param.names"={return(x@param.names)},
            "omega"={return(x@omega)},
            "omega.estim"={return(x@omega.estim)},
            "omega.model"={return(x@omega.model)},
            "subomega"={return(x@subomega)},
            "idvec.estim"={return(x@idvec.estim)},
            "idcol.eta"={return(x@idcol.eta)},
            "idcol.eta.fix"={return(x@idcol.eta.fix)},
            "idvec.var"={return(x@idvec.var)},
            "idvec.cov"={return(x@idvec.cov)},
            "idmat.var"={return(x@idmat.var)},
            "idmat.cov"={return(x@idmat.cov)},
            "idmat.estim"={return(x@idmat.estim)},
            "omega.hat"={return(x@omega.hat)},
            "param.se"={return(x@param.se)},
            "omega.var"={return(x@omega.var)},
            "conf.int"={return(x@conf.int)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixVarModelHat" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.level"={x@name.level<-value},
            "variable"={x@variable<-value},
            "log"={x@log<-value},
            "nphi"={x@nphi<-value},
            "param"={x@param<-value},
            "param.names"={x@param.names<-value},
            "omega"={x@omega<-value},
            "omega.estim"={x@omega.estim<-value},
            "omega.model"={x@omega.model<-value},
            "subomega"={x@subomega<-value},
            "idvec.estim"={x@idvec.estim<-value},
            "idcol.eta"={x@idcol.eta<-value},
            "idcol.eta.fix"={x@idcol.eta.fix<-value},
            "idmat.var"={x@idmat.var<-value},
            "idmat.cov"={x@idmat.cov<-value},
            "idmat.estim"={x@idmat.estim<-value},
            "omega.hat"={x@omega.hat<-value},
            "param.se"={x@param.se<-value},
            "omega.var"={x@omega.var<-value},
            "conf.int"={x@conf.int<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

################################################################################
# Print/show functions for SaemixVarModel  (statistical model)

setMethod("print","SaemixVarModel",
          function(x,nlines=10,...) {
            cat("Variability level",x@name.level,"associated with variable",x@variable,"\n")
# if variable moves to SaemixModel
#            cat("Variability level:",x@name.level,":\n")
            print(x@omega.model)
          }
)

setMethod("show","SaemixVarModel",
          function(object) {
            cat("Model for variability level",object@name.level,":\n")
            print(object@omega.model)
          }
)
setMethod("showall","SaemixVarModel",
          function(object) {
            cat("Model for variability level",object@name.level,":\n")
            print(object@omega.model)
            cat("      Estimated elements:\n")
            print(object@omega.estim)
            cat("      CI:\n")
            print(object@omega)
            cat("      parameter names:",object@param.names,"\n")
          }
)
################################################################################
# Print/show functions for SaemixVarModelHat (estimates) 


setMethod("print","SaemixVarModelHat",
          function(x,nlines=10,...) {
            cat("Variability level",x@name.level,"associated with variable",x@variable,"\n")
            # if variable moves to SaemixModel
            #            cat("Variability level:",x@name.level,":\n")
            print(x@conf.int)
          }
)

setMethod("show","SaemixVarModelHat",
          function(object) {
            cat("Estimated parametres for variability level",object@name.level,":\n")
            print(object@conf.int)
          }
)
setMethod("showall","SaemixVarModelHat",
          function(object) {
            cat("Model for variability level",object@name.level,":\n")
            print(object@omega.model)
            cat("      Estimated elements:\n")
            print(object@omega.estim)
            cat("      CI:\n")
            print(object@omega)
            cat("      parameter names:",object@param.names,"\n")
            cat("      estimates:\n")
            print(object@conf.int)
          }
)

