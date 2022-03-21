
#' Class "SaemixVarLevel"
#' 
#' An object of the SaemixVarLevel class, representing an observation model for an outcome. Outcomes can be discrete or survival-type data
#' (SaemixDiscreteOutcome) or continuous (SaemixContinuousOutcome), with both subclasses inheriting from the SaemixVarLevel class
#' 
#' @name SaemixVarLevel-class 
#' @docType class
#' @aliases SaemixVarLevel SaemixVarLevel-class 
#' @aliases print,SaemixVarLevel showall,SaemixVarLevel show,SaemixVarLevel
#' @aliases SaemixDiscreteOutcome SaemixDiscreteOutcome-class 
#' @aliases print,SaemixDiscreteOutcome showall,SaemixDiscreteOutcome show,SaemixDiscreteOutcome
#' @aliases SaemixContinuousOutcome SaemixContinuousOutcome-class 
#' @aliases print,SaemixContinuousOutcome showall,SaemixContinuousOutcome show,SaemixContinuousOutcome
#' 
#' @section Objects from the Class: 
#' An object of the SaemixVarLevel class contains the following slots:
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
#' showClass("SaemixVarLevel")
#' 
#' @keywords classes
#' @exportClass SaemixVarLevel


# Variability level class - generic
setClass(Class = "SaemixVarLevel",
         representation=representation(
           name.level = "character", # name of variability level
           variable = "character" # which variable (in the dataset) is the variability associated to
         ),
         validity=function(object){
           return(TRUE)
         }
)


setMethod( 
  f="initialize",
  signature="SaemixVarLevel",
  definition=function(.Object, name.level="iiv", variable="id"){
    .Object@name.level <- name.level
    .Object@variable <- variable
    validObject(.Object)
    return(.Object)
  }
)


# Getteur
setMethod(
  f ="[",
  signature = "SaemixVarLevel" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.level"={return(x@name.level)},
            "variable"={return(x@variable)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixVarLevel" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.level"={x@name.level<-value},
            "variable"={x@variable<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)
