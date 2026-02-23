#' @include aaa_generics.R
NULL

#' Class "SaemixIndivPar"
#' 
#' An object of the SaemixIndivPar class, containing the individual parameters for a given level of variability, on the eta scale
#' 
#' @name SaemixIndivPar-class 
#' @docType class
#' @aliases SaemixIndivPar SaemixIndivPar-class 
#' @aliases print,SaemixIndivPar showall,SaemixIndivPar show,SaemixIndivPar
#' 
#' @section Objects from the Class: 
#' An object of the SaemixIndivPar class contains the following slots:
#' @slot phi individual parameters at the global level
#' @slot psi individual parameters at the global level, on the model scale (psi=h(phi))
#' @slot mean.phi mean (population) individual parameters (no random effects)
#' @slot log A string recording the warnings and messages during the creation of the object
#' @slot phi.samp sampled individual parameters 
#' @slot psi.samp sampled individual parameters, on the model scale
#' @slot varlevel vector of variability levels
#' @slot varphi list of SaemixVarPhi objects at each variability level
#' @slot index indicator to match the lines in phi in the dataframe to the corresponding lines in the full dataset (ex: for IIV, this would be 1 for all lines corresponding to subject 1, 2 for all lines of subject 2, etc...; for IOV, 1.1 would be occasion 1 in subject 1, 1.2 occasion 2 in subject 1, 2.1 occasion 1 in subject 2, etc...)
#' @section Methods:
#'   \describe{
#'     \item{[<-}{\code{signature(x = "SaemixIndivPar")}: replace elements of object}
#'     \item{[}{\code{signature(x = "SaemixIndivPar")}: access elements of object}
#'     \item{initialize}{\code{signature(.Object = "SaemixIndivPar")}: internal function to initialise object, not to be used}
#'     \item{print}{\code{signature(x = "SaemixIndivPar")}: prints details about the object (more extensive than show)}
#'     \item{showall}{\code{signature(object = "SaemixIndivPar")}: shows all the elements in the object}
#'     \item{show}{\code{signature(object = "SaemixIndivPar")}: prints details about the object}
#' 	 }
#' 
#' @details
#' Intermediate object holding the levels 
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' 
#' @examples
#' showClass("SaemixIndivPar")
#' 
#' @keywords classes
#' @exportClass SaemixIndivPar

# Variability level class - generic
setClass(Class = "SaemixIndivPar",
         representation=representation(
           phi = "matrix", # individual parameters (eta scale)
           psi = "matrix", # individual parameters (model scale)
           mean.phi = "matrix", # mean individual parameter
           log="character",		# A record of the warnings and messages during the creation of the object
           phi.samp = "array", # sampled individual parameters
           psi.samp = "array", # sampled individual parameters
           varlevel = "character", # vector giving the names of the variability level
           varphi = "SaemixVarPhi", # list of SaemixVarPhi objects # Check if list or saemixVarPhi
           index = "numeric" # index to match to the dataset
         ),
         validity=function(object){
           return(TRUE)
         }
)

# Getteur
setMethod(
  f ="[",
  signature = "SaemixIndivPar" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "phi"={return(x@phi)},
            "psi"={return(x@psi)},
            "mean.phi"={return(x@mean.phi)},
            "log"={return(x@log)},
            "psi.samp"={return(x@psi.samp)},
            "phi.samp"={return(x@phi.samp)},
            "varlevel"={return(x@varlevel)},
            "varphi"={return(x@varphi)},
            "index"={return(x@index)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixIndivPar" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "phi"={x@phi<-value},
            "psi"={x@psi<-value},
            "mean.phi"={x@mean.phi<-value},
            "log"={x@log<-value},
            "psi.samp"={x@psi.samp<-value},
            "phi.samp"={x@phi.samp<-value},
            "varlevel"={x@varlevel<-value},
            "varphi"={x@varphi<-value},
            "index"={x@index<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)


setMethod(
  f="initialize",
  signature="SaemixIndivPar",
  definition=function(.Object, varphi=list()) {
    # check all objects from varphi are SaemixVarPhi type and are not empty
    nlev <- length(varphi)
    # 
    l1<-c()
    for(i in 1:length(varphi)) l1<-c(l1, varphi[[i]]@name)
    .Object@varlevel <- l1
    .Object@varphi <- varphi
    # Compute phi from the phi in varphi elements by replicating and summing
    
    # index=index of the lowest level
    .Object@index <- varphi[[nlev]]@index
    validObject(.Object)
    return(.Object)
  }
)

