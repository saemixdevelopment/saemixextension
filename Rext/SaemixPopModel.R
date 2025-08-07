#' @include aaa_generics.R
#' @include SaemixParameter.R
#' @include SaemixParameter-methods.R
NULL

# (phi model without eta)

#' Class "SaemixPopModel"
#' 
#' An object of the SaemixPopModel class, representing the model for a variability level. 
#' The SaemixPopModelHat class inherits from the SaemixPopModel class and adds slots for the estimated parameters.
#' 
#' @name SaemixPopModel-class 
#' @docType class
#' @aliases SaemixPopModel SaemixPopModel-class 
#' @aliases print,SaemixPopModel showall,SaemixPopModel show,SaemixPopModel
#' 
#' @aliases SaemixPopModelHat SaemixPopModelHat-class 
#' @aliases print,SaemixPopModelHat showall,SaemixPopModelHat show,SaemixPopModelHat
#' 
#' @section Objects from the Class: 
#' An object of the SaemixPopModel class contains the following slots:
#' @slot name.level string giving the name of variability level (defaults to "iiv" as the subject-level grouping)
#' @slot nphi number of model parameters (number of columns of phi)
#' @slot phi.model mu and beta model matrix (0/1 for parameters in the model or nor)
#' @slot phi initial values for the parameters, as a matrix (initialised to the CI from mu.init and covariate.init given in parameters)
#' @slot phi.estim model
#' @slot param etc...
#' 
#' 
#' @section Methods:
#'   \describe{
#'     \item{[<-}{\code{signature(x = "SaemixPopModel")}: replace elements of object}
#'     \item{[}{\code{signature(x = "SaemixPopModel")}: access elements of object}
#'     \item{initialize}{\code{signature(.Object = "SaemixPopModel")}: internal function to initialise object, not to be used}
#'     \item{print}{\code{signature(x = "SaemixPopModel")}: prints details about the object (more extensive than show)}
#'     \item{showall}{\code{signature(object = "SaemixPopModel")}: shows all the elements in the object}
#'     \item{show}{\code{signature(object = "SaemixPopModel")}: prints details about the object}
#' 	 }
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' 
#' @seealso \code{\link{saemixData}} \code{\link{SaemixModel}} \code{\link{saemixControl}} \code{\link{saemix}}
#' @examples
#' showClass("SaemixPopModel")
#' 
#' @keywords classes
#' @exportClass SaemixPopModel
#' @exportClass SaemixPopModelHat

# Population model - mu, beta, covariate model for one level of variability
setClass(Class = "SaemixPopModel",
         representation=representation(
           name.level = "character", # name of variability level
           log = "character", # warning messages
           nphi = "numeric", # number of model parameters
           param = "numeric", # vector of population parameters in the model (mu and beta)
           param.names="character", # names of the parameters to be estimated (mu and beta)
           phi.model = "matrix", # parameter + covariate model matrix of 0/1, 1 indicating that the parameter is present in the model; the first row represents the fixed effect and subsequent rows the covariate models
           phi.estim = "matrix", # same matrix indicating whether elements are estimated, fixed or prior;
           phi = "matrix", # CI for mu and beta, as a matrix
           idvec.mu= "numeric", # indices of mu terms in param
           idvec.beta= "numeric", # indices of beta terms in param
           idvec.estim= "numeric", # indices of elements in param to be estimated (the other parameters are fixed)
           idmat.mu= "numeric", # indices of mu terms in phi.model
           idmat.beta= "numeric", # indices of beta terms in phi.model
           idmat.estim= "numeric" # indices of elements estimated in phi.model, as a vector (by column)
         ),
         validity=function(object){
           # Check all sizes are commensurate TODO
           return(TRUE)
         }
)
# Getteur
setMethod(
  f ="[",
  signature = "SaemixPopModel" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.level"={return(x@name.level)},
            "log"={return(x@log)},
            "nphi"={return(x@nphi)},
            "param"={return(x@param)},
            "param.names"={return(x@param.names)},
            "phi.model"={return(x@phi.model)},
            "phi.estim"={return(x@phi.estim)},
            "phi"={return(x@phi)},
            "idvec.mu"={return(x@idvec.mu)},
            "idvec.beta"={return(x@idvec.beta)},
            "idvec.estim"={return(x@idvec.estim)},
            "idmat.mu"={return(x@idmat.mu)},
            "idmat.beta"={return(x@idmat.beta)},
            "idmat.estim"={return(x@idmat.estim)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixPopModel" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.level"={x@name.level<-value},
            "log"={x@log<-value},
            "nphi"={x@nphi<-value},
            "param"={x@param<-value},
            "param.names"={x@param.names<-value},
            "phi.model"={x@phi.model<-value},
            "phi.estim"={x@phi.estim<-value},
            "phi"={x@phi<-value},
            "idvec.mu"={x@idvec.mu<-value},
            "idvec.beta"={x@idvec.beta<-value},
            "idvec.estim"={x@idvec.estim<-value},
            "idmat.mu"={x@idmat.mu<-value},
            "idmat.beta"={x@idmat.beta<-value},
            "idmat.estim"={x@idmat.estim<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

setMethod( 
  f="initialize",
  signature="SaemixPopModel",
  definition=function(.Object, phi.model, name.level="iiv", verbose=FALSE){
    if(missing(phi.model)) {
      if(verbose) cat("Please supply the model structure\n")
      return("Creation of SaemixPopModel object failed \n")
    }
    .Object@phi.model <- phi.model
    .Object@nphi<-dim(phi.model)[2]
    .Object@name.level <- name.level
    validObject(.Object)
    return(.Object)
  }
)
################################################################################
# Child class with estimates added

setClass(Class = "SaemixPopModelHat",
         contains = "SaemixPopModel",
         representation=representation(
           phi.hat = "matrix", # estimated parameters, in matrix form
           phi.se = "matrix", # estimated SE for parameters, in matrix form
           param.hat = "numeric", # estimated parameters, vector in order of param
           param.se = "numeric", # estimated SE on parameters, vector in order of param
           param.fim = "matrix", # estimated variance-covariance matrix of estimation for fixed effects block (square matrix of size length(estimated parameters))
           conf.int = "data.frame" # Table giving the estimates, SE, CV and confidence intervals (assuming normality of the estimators)
         ),
         validity=function(object){
           # Check all sizes are commensurate TODO
           return(TRUE)
         }
)
# Getteur
setMethod(
  f ="[",
  signature = "SaemixPopModelHat" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.level"={return(x@name.level)},
            "log"={return(x@log)},
            "nphi"={return(x@nphi)},
            "param"={return(x@param)},
            "param.names"={return(x@param.names)},
            "phi.model"={return(x@phi.model)},
            "phi.estim"={return(x@phi.estim)},
            "phi"={return(x@phi)},
            "idvec.mu"={return(x@idvec.mu)},
            "idvec.beta"={return(x@idvec.beta)},
            "idvec.estim"={return(x@idvec.estim)},
            "idmat.mu"={return(x@idmat.mu)},
            "idmat.beta"={return(x@idmat.beta)},
            "idmat.estim"={return(x@idmat.estim)},
            "phi.hat"={return(x@phi.hat)},
            "phi.se"={return(x@phi.se)},
            "param.hat"={return(x@param.hat)},
            "param.se"={return(x@param.se)},
            "param.fim"={return(x@param.fim)},
            "conf.int"={return(x@conf.int)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixPopModelHat" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.level"={x@name.level<-value},
            "log"={x@log<-value},
            "nphi"={x@nphi<-value},
            "param"={x@param<-value},
            "param.names"={x@param.names<-value},
            "phi.model"={x@phi.model<-value},
            "phi.estim"={x@phi.estim<-value},
            "phi"={x@phi<-value},
            "idvec.mu"={x@idvec.mu<-value},
            "idvec.beta"={x@idvec.beta<-value},
            "idvec.estim"={x@idvec.estim<-value},
            "idmat.mu"={x@idmat.mu<-value},
            "idmat.beta"={x@idmat.beta<-value},
            "idmat.estim"={x@idmat.estim<-value},
            "phi.hat"={x@phi.hat<-value},
            "phi.se"={x@phi.se<-value},
            "param.hat"={x@param.hat<-value},
            "param.se"={x@param.se<-value},
            "param.fim"={x@param.fim<-value},
            "conf.int"={x@conf.int<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

setMethod( 
  f="initialize",
  signature="SaemixPopModelHat",
  definition=function(.Object, saemix.popmodel, verbose=FALSE){
    .Object <- callNextMethod(.Object,phi.model=saemix.popmodel@phi.model, name.level=saemix.popmodel@name.level, verbose=verbose)
    for(i in slotNames(saemix.popmodel))
      slot(.Object,i)<-slot(saemix.popmodel,i)
    if(length(.Object@param.names)>0) {
      tabhat <- data.frame(Parameter=.Object@param.names)
      if(length(.Object@param.hat)>0) {
        tabhat$Estimate <- .Object@param.hat
      } else  {
        tabhat$Initial.Value=.Object@param
        status<-rep("fixed", length(.Object@param))
        if(length(.Object@idvec.estim)>0) status[.Object@idvec.estim]<-"estimated"
#        status[.Object@idvec.estim=="prior"]<-"prior"
        tabhat$Status<-status
      }
      if(length(.Object@param.hat)>0) {
        tabhat$SE <- .Object@param.SE
      }
      .Object@conf.int <- tabhat
    }
    validObject(.Object)
    return(.Object)
  }
)

################################################################################
# Print/show functions for SaemixPopModel and SaemixPopModelHat (phi model without eta)

setMethod("print","SaemixPopModel",
          function(x,nlines=10,...) {
            cat("Model of fixed effects for level",x@name.level,"\n")
            # if variable moves to SaemixVarLevelHat
            #            cat("Variability level:",x@name.level,":\n")
            print(x@phi.model)
          }
)

setMethod("show","SaemixPopModel",
          function(object) {
            cat("Model of fixed effects for level",object@name.level,":\n")
            print(object@phi.model)
          }
)
setMethod("showall","SaemixPopModel",
          function(object) {
            cat("Model of fixed effects for level",object@name.level,":\n")
            print(data.frame(Parameter=object@param.names, CI=object@param, Status=object@phi.estim[object@phi.model==1]),row.names=FALSE)
            cat(object@log)
          }
)


setMethod("print","SaemixPopModelHat",
          function(x,nlines=10,...) {
            cat("Estimated fixed effects for level",x@name.level,"\n")
            # if variable moves to SaemixVarLevelHat
            #            cat("Variability level:",x@name.level,":\n")
            print(x@conf.int, row.names=FALSE)
          }
)

setMethod("show","SaemixPopModelHat",
          function(object) {
            cat("Estimated fixed effects for level",object@name.level,":\n")
            print(object@phi.model)
          }
)
setMethod("showall","SaemixPopModelHat",
          function(object) {
            cat("Estimated fixed effects for level",object@name.level,":\n")
            print(object@conf.int, row.names=FALSE)
            cat(object@log)
          }
)

################################################################################
# Create individual model (mu and beta)

# create fixed effect model for each level of variability
extractFixedEffectModel <- function(parameters, varlevel=c()) {
  xcheck<-checkParameters(parameters=parameters, varlevel=varlevel)
  if(inherits(xcheck,"character")) {
    return("Extraction failed: parameters should be a list of SaemixParameter objects\n")
  }
  parameters <- xcheck$parameters
  varlevel<-xcheck$varlevel.order
  logmsg<-xcheck$logmsg
  
  nphi<-length(parameters)
  nampar <- names(parameters)
  list.indmodel<-vector(length=length(varlevel),mode="list")
  names(list.indmodel)<-varlevel
  # list of covariates
  covariate<-c()
  for(ipar in 1:nphi) 
    if(length(parameters[[ipar]]@covariate)>0) covariate<-c(covariate, parameters[[ipar]]@covariate)
  if(length(covariate)>0) covariate<-unique(covariate) else covariate<-c()
  for(ilev in 1:length(varlevel)) {
    xlev <- varlevel[ilev] # name of the variability level 
    if(ilev==1) {
      covmodel.init <- covmodel<-matrix(rep(1,nphi),nrow=1,dimnames=list(c("pop"), nampar))
      covmodel.estim<-matrix(rep("estimated",nphi),nrow=1,dimnames=list(c("pop"), nampar))
      } else {
        covmodel.init<-covmodel<-matrix(rep(0,nphi),nrow=1,dimnames=list(c("pop"), nampar))
        covmodel.estim<-matrix(rep("fixed",nphi),nrow=1,dimnames=list(c("pop"), nampar))
      }
    covmodel<-rbind(covmodel,
                    matrix(data=0,nrow=length(covariate),ncol=nphi, dimnames = list(covariate, nampar)))
    covmodel.init<-rbind(covmodel.init,
                    matrix(data=0,nrow=length(covariate),ncol=nphi, dimnames = list(covariate, nampar)))
    covmodel.estim<-rbind(covmodel.estim,
                    matrix(data="estimated",nrow=length(covariate),ncol=nphi, dimnames = list(covariate, nampar)))
    for(ipar in 1:nphi) {
      if(length(parameters[[ipar]]@mu.init)>=ilev) {
        if(parameters[[ipar]]@mu.estim[ilev]=="estimated") {
          covmodel[1,ipar]<-1
          covmodel.init[1,ipar]<-parameters[[ipar]]@invtransform(parameters[[ipar]]@mu.init[ilev])
          covmodel.estim[1,ipar]<-"estimated"
        }
        if(parameters[[ipar]]@mu.estim[ilev]=="fixed" & parameters[[ipar]]@invtransform(parameters[[ipar]]@mu.init[ilev])!=0) {
          covmodel[1,ipar]<-1
          covmodel.init[1,ipar]<-parameters[[ipar]]@invtransform(parameters[[ipar]]@mu.init[ilev])
        }
      }
      if(length(parameters[[ipar]]@covariate)>0) {
        idx.cov<-which(parameters[[ipar]]@covariate.varlevel==varlevel[ilev])
        if(length(idx.cov)>0) {
          covariate2 <- parameters[[ipar]]@covariate[idx.cov]
          icov<-match(covariate2, rownames(covmodel))
          covmodel[icov,ipar]<-1
          covmodel.estim[icov,ipar]<-parameters[[ipar]]@covariate.estim[idx.cov]
          covmodel.init[icov,ipar]<-parameters[[ipar]]@covariate.init[idx.cov]
        }
      }
    }
    # Remove empty lines
    all0<-which(rowSums(covmodel)==0)
    if(length(all0)>0) {
      covmodel<-covmodel[-c(all0),,drop=FALSE]
      covmodel.estim<-covmodel.estim[-c(all0),,drop=FALSE]
      covmodel.init<-covmodel.init[-c(all0),,drop=FALSE]
    }
    vlev <-  new(Class="SaemixPopModel", phi.model=covmodel, name.level=varlevel[ilev])
    vlev@phi.estim <- covmodel.estim
    vlev@phi <- covmodel.init
    vlev <- addSaemixIndices(vlev)
    list.indmodel[[ilev]] <- vlev
    
  }
  return(list.indmodel)
}

## idvec.mu: index of mu parameters in param
## idvec.beta: index of beta parameters in param
## idvec.estim: index of estimated parameters in param (others=fixed)

## idcol.eta.fix: parameters with fixed variability
## idvec.var: variance parameters in param
## idvec.cov: covariance parameters 

setMethod("addSaemixIndices",
          signature="SaemixPopModel",
          function(object, verbose=FALSE) {
            matname<-matrix(data="",nrow=dim(object@phi.model)[1], ncol=dim(object@phi.model)[2])
            for(i in 1:dim(object@phi.model)[1]) {
              if(rownames(object@phi.model)[i]=="pop") matname[i,]<-paste0("mu.",colnames(object@phi.model)) else matname[i,]<-paste0("beta.",rownames(object@phi.model)[i],".",colnames(object@phi.model))
            }
            idmat.par<-which(object@phi.model==1)
            object@param.names<-matname[idmat.par]
            object@param<-object@phi[idmat.par]
            x1 <- grep("mu",object@param.names)
            if(length(x1)>0) object@idvec.mu<-x1
            x1 <- grep("beta",object@param.names)
            if(length(x1)>0) object@idvec.beta<-1
            x1<-which(object@phi.estim[idmat.par]=="estimated")
            if(length(x1)>0) object@idvec.estim<-x1
            
            # not sure we need those indices
            x1<-grep("mu",matname)
            if(length(x1)>0) object@idmat.mu<-x1
            x1<-grep("beta",matname)
            if(length(x1)>0) object@idmat.beta<-x1
            x1<-which(object@phi.estim=="estimated")
            if(length(x1)>0) object@idmat.estim<-x1
            return(object)
          }
)
