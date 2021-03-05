####################################################################################
####      SaemixModel class - definition        ####
####################################################################################
#' Class "SaemixModel"
#' 
#' An object of the SaemixModel class, representing a nonlinear mixed-effect
#' model structure, used by the SAEM algorithm.
#' 
#' @name SaemixModel-class
#' @docType class
#' @aliases SaemixModel-class SaemixModel [<-,SaemixModel-method
#' @aliases print,SaemixModel showall,SaemixModel show,SaemixModel summary,SaemixModel 
#' @section Objects from the Class: 
#' An object of the SaemixModel class can be created by using the function \code{\link{saemixModel}} and contain the following slots:
#'   \describe{
#'     \item{\code{model}:}{Object of class \code{"function"}: name of the function used to get predictions from the model (see the User Guide and the online examples for the format and what this function should return).}
#'     \item{\code{description}:}{Object of class \code{"character"}: an optional text description of the model}
#'     \item{\code{psi0}:}{Object of class \code{"matrix"}: a matrix with named columns containing the initial estimates for the parameters in the model (first line) and for the covariate effects (second and subsequent lines, optional). The number of columns should be equal to the number of parameters in the model.}
#'     \item{\code{transform.par}:}{Object of class \code{"numeric"}: vector giving the distribution for each model parameter (0: normal, 1: log-normal, 2: logit, 3: probit). Its length should be equal to the number of parameters in the model.}
#'     \item{\code{fixed.estim}:}{Object of class \code{"numeric"}: for each parameter, 0 if the parameter is fixed and 1 if it should be estimated. Defaults to a vector of 1 (all parameters are estimated). Its length should be equal to the number of parameters in the model.}
#'     \item{\code{error.model}:}{Object of class \code{"character"}: name of the error model. Valid choices are "constant" (default), "proportional" and "combined" (see equations in User Guide)}
#'     \item{\code{covariate.model}:}{Object of class \code{"matrix"}: a matrix of 0's and 1's, with a 1 indicating that a parameter-covariate relationship is included in the model (and an associated fixed effect will be estimated). The nmuber of columns should be equal to the number of parameters in the model and the number of rows to the number of covariates.}
#'     \item{\code{covariance.model}:}{Object of class \code{"matrix"}: a matrix f 0's and 1's giving the structure of the variance-covariance matrix. Defaults to the Identity matrix (diagonal IIV, no correlations between parameters)}
#'     \item{\code{omega.init}:}{Object of class \code{"matrix"}: a matrix giving the initial estimate for the variance-covariance matrix}
#'     \item{\code{error.init}:}{Object of class \code{"numeric"}: a vector giving the initial estimate for the parameters of the residual error}
#'   }
#'   Additional elements are added to the model object after a call to \code{saemix} and are used in the algorithm.
#' @section Methods:
#'   \describe{
#'     \item{[<-}{\code{signature(x = "SaemixModel")}: replace elements of object}
#'     \item{[}{\code{signature(x = "SaemixModel")}: access elements of object}
#'     \item{initialize}{\code{signature(.Object = "SaemixModel")}: internal function to initialise object, not to be used}
#'     \item{plot}{\code{signature(x = "SaemixModel")}: plot predictions from the model}
#'     \item{print}{\code{signature(x = "SaemixModel")}: prints details about the object (more extensive than show)}
#'     \item{showall}{\code{signature(object = "SaemixModel")}: shows all the elements in the object}
#'     \item{show}{\code{signature(object = "SaemixModel")}: prints details about the object}
#'   }
#' @references Comets  E, Lavenu A, Lavielle M. Parameter estimation in nonlinear mixed effect models using saemix, an R implementation of the SAEM algorithm. Journal of Statistical Software 80, 3 (2017), 1-41.
#' 
#' Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.
#' 
#' Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece (2011), Abstr 2173.
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' @author Audrey Lavenu
#' @author Marc Lavielle.
#' @seealso \code{\link{SaemixData}} \code{\link{SaemixObject}} \code{\link{saemixControl}} \code{\link{saemix}}
#' \code{\link{plot.saemix}}
#' @keywords classes
#' @exportClass SaemixModel
#' @examples
#' 
#' showClass("SaemixModel")
#' 

setClass(
  Class="SaemixModel",
  representation=representation(
    model="function",     # name of model function
    description="character",  # model description
    modeltype="character",     # type of model (structural, for continuous responses, or likelihood)
    psi0="matrix",    # CI for parameter estimates
    transform.par="numeric",  # distribution for model parameters
    fixed.estim="numeric",  # 1 for fixed parameters estimated
    error.model="character",  # residual error model
    covariate.model="matrix", # covariate model
    betaest.model="matrix", # 1st line=ones, next lines=covariate model
    covariance.model="matrix",  # covariance model
    omega.init="matrix",  # CI for Omega
    error.init="numeric", # CI for residual error
    nb.parameters="integer",  # nb of parameters in the model
    name.modpar="character",  # name of parameters in the model (columns of psi0)
    name.fixed="character", # name of fixed parameters
    name.random="character",  # name of random parameters
    name.sigma="character", # name of residual parameters (maybe not necessary)
    name.predictors="character",# name of predictors 
    name.X="character", # name of X 
    name.response="character",  # name of response
    name.cov="character", # name of covariates
    indx.fix="numeric",   # index of mean param estimated (was indx.betaI)
    indx.cov="numeric",   # index of cov param estimated (was indx.betaC)
    indx.omega="numeric", # index of random param estimated (was i1.omega2)
    indx.res="numeric",   # index of param of residual errors estimated (was indx.res)
    Mcovariates="data.frame"  # matrix of individual covariates in the model
  ),
  validity=function(object){
#    cat ("--- Checking SaemixModel object ---\n")
    if (dim(object@psi0)[1]==0) {
      message("[ SaemixModel : Error ] Please provide initial estimates for the fixed effect (a matrix with columns named after the parameters in the model).")
      return("Missing psi0")
    }
    isize<-0
    npar<-dim(object@psi0)[2]
    if(npar!=length(object@transform.par)) isize<-1
    if(npar!=length(object@fixed.estim)) isize<-1
    if (npar!=dim(object@covariate.model)[2]) isize<-1
    if (npar!=dim(object@covariance.model)[1]) isize<-1
    if (npar!=dim(object@omega.init)[1]) isize<-1
#    cat("npar=",npar,length(object@transform.par),length(object@fixed.estim), dim(object@covariate.model)[2],dim(object@covariance.model)[1],dim(object@omega.init)[1],"\n")
    if(isize==1) {
      message("[ SaemixModel : Error ] The number of parameters should be the same in the following elements: psi0 (initial conditions), transform.par, fixed.estim, covariate.model, and the matrices covariance.model and omega.init should be square matrices of size equal to the number of parameters. Please check the input.")
      return("Size mismatch")
    }
    if(npar<2) {
      message("[ SaemixModel : Error ] SAEM needs at least two parameters to work on.")
      return("Psi0 has size 1")
    }
    if(sum(object@fixed.estim*mydiag(object@covariance.model))==0) {
      message("[ SaemixModel : Error ] ")
      if(sum(mydiag(object@covariance.model))==0) message("At least one parameter with IIV must be included in the model.") else message("At least one parameter with IIV must be estimated and not fixed in the model.")
      return("Invalid IIV structure")
    }
    if(is.na(sum(match(object@error.model,c('constant','proportional','combined', 'exponential'))))) {
      message("[ SaemixModel : Error ] Invalid residual error model")
      return("Invalid residual error model")
    }
    if(is.na(match(object@modeltype,c("structural","likelihood")))) {
      message("[ SaemixModel : Error ] Invalid type of model")
      return("Invalid model type")
    }
    return(TRUE)
  }
)

#' @rdname initialize-methods
#' 
#' @param model name of the function used to compute the structural model. The
#' function should return a vector of predicted values given a matrix of
#' individual parameters, a vector of indices specifying which records belong
#' to a given individual, and a matrix of dependent variables (see example
#' below).
#' @param description a character string, giving a brief description of the
#' model or the analysis
#' @param modeltype a character string, giving the type of the model for the analysis (one of "structural" or "likelihood", defaults to structural)
#' @param psi0 a matrix with a number of columns equal to the number of
#' parameters in the model, and one (when no covariates are available) or two
#' (when covariates enter the model) giving the initial estimates for the fixed
#' effects. The column names of the matrix should be the names of the
#' parameters in the model, and will be used in the plots and the summaries.
#' When only the estimates of the mean parameters are given, psi0 may be a
#' named vector.
## #' @param name.response a character string or a column number specifying which column of the data contains the dependent variable
#' @param name.sigma a vector of character string giving the names of the residual error parameters (defaults to "a" and "b")
#' @param transform.par the distribution for each parameter (0=normal,
#' 1=log-normal, 2=probit, 3=logit). Defaults to a vector of 1s (all parameters
#' have a log-normal distribution)
#' @param fixed.estim whether parameters should be estimated (1) or fixed to
#' their initial estimate (0). Defaults to a vector of 1s
#' @param error.model type of residual error model (valid types are constant,
#' proportional, combined and exponential). Defaults to constant
#' @param covariate.model a matrix giving the covariate model. Defaults to no
#' covariate in the model
#' @param covariance.model a square matrix of size equal to the number of parameters in the model, 
#' giving the variance-covariance matrix of the model: 1s correspond to estimated variances (in the diagonal) 
#' or covariances (off-diagonal elements). Defaults to the identity matrix
#' @param omega.init a square matrix of size equal to the number of parameters
#' in the model, giving the initial estimate for the variance-covariance matrix
#' of the model. Defaults to the identity matrix
#' @param error.init a vector of size 2 giving the initial value of a and b in
#' the error model. Defaults to 1 for each estimated parameter in the error model
#' @param name.modpar names of the model parameters, if they are not given as
#' the column names (or names) of psi0
#' 
#' @exportMethod initialize

setMethod(
  f="initialize",
  signature="SaemixModel",
  definition=function(.Object,model,description,modeltype,psi0, name.response, name.sigma, transform.par,fixed.estim, error.model,covariate.model,covariance.model,omega.init,error.init, name.modpar, verbose=TRUE){
#    cat ("--- initialising SaemixModel Object --- \n")
    if(missing(model)) {
#      cat("Error initialising SaemixModel object:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
      return(.Object)
    }
    .Object@model<-model
    if(missing(description)) description<-""
    .Object@description<-description
    if(missing(modeltype)) modeltype<-"structural"
    .Object@modeltype<-modeltype
    if(missing(psi0) || length(psi0)==0) {
      if(verbose) message("Error initialising SaemixModel object:\n   Please provide initial estimates for the fixed effect (a matrix with columns named after the parameters in the model).\n")
      return(.Object)
    }
    npar<-dim(psi0)[2]
    if(missing(name.modpar) || length(name.modpar)==0) {
      y1<-try(name.modpar<-colnames(psi0))
      if(inherits(y1,"try-error")) {
        if(verbose) message("     Can't find parameter names.\n")
        name.modpar<-paste("theta",1:npar)
      }
    }
    if(is.null(colnames(psi0))) {
      y1<-try(colnames(psi0)<-name.modpar)
      if(inherits(y1,"try-error")) {
        if(verbose) message("Warning:\n   Problem with names of psi0\n")
        colnames(psi0)<-name.modpar<-paste("theta",1:npar)
      }
    }
    if(missing(name.response)) name.response<-""
    .Object@name.response<-name.response
    if(!missing(covariate.model)) {
      if(dim(psi0)[1]<2 & sum(covariate.model)>0){
        psi0<-rbind(psi0,rep(0,dim(psi0)[2]))
      }
    }
    if(is.null(rownames(psi0))) {
      rownames(psi0)<-rep("",dim(psi0)[1])
      rownames(psi0)[1]<-"Pop.CondInit"
      if(dim(psi0)[1]>1) rownames(psi0)[2:dim(psi0)[1]]<-"Cov.CondInit"
    }
    .Object@psi0<-psi0    
    .Object@name.modpar<-name.modpar
    if(missing(error.model) || length(error.model)==0) error.model<-"constant"
    if(sum(!error.model %in% c('constant','proportional','combined', 'exponential'))) {
      message("Invalid error model, switching to constant")
      error.model[!error.model %in% c('constant','proportional','combined', 'exponential')] <- "constant"
    }
    if(length(error.model)<length(name.response)) error.model<-rep(error.model,length.out=length(name.response))
    .Object@error.model<-error.model
# Checking sizes
    .Object@nb.parameters<-npar
    if(missing(transform.par) || length(transform.par)==0) transform.par<-rep(0,npar)
    .Object@transform.par<-transform.par
    if(missing(fixed.estim) || length(fixed.estim)==0) fixed.estim<-rep(1,npar)
    .Object@fixed.estim<-fixed.estim
    if(missing(covariate.model) || length(covariate.model)==0 || sum(covariate.model)==0) covariate.model<-matrix(nrow=0,ncol=npar)
    if(is.null(dim(covariate.model)) & length(covariate.model)>0) covariate.model<-matrix(covariate.model,byrow=T,ncol=npar) # Covariate model given as a vector
    if(is.null(colnames(covariate.model))) colnames(covariate.model)<-colnames(psi0)
    .Object@covariate.model<-covariate.model
    if(missing(covariance.model) || length(covariance.model)==0) {
      covariance.model<-diag(nrow=npar,ncol=npar)
    } else {
      if(dim(covariance.model)[1]!=dim(covariance.model)[2]) {
        if(verbose) message("Error initialising SaemixModel object:\n   The covariance model needs to be a square matrix, please check dimensions.\n")
      return(.Object)
      }
    }
    nomg<-dim(covariance.model)[1]
    if(nomg!=npar) {
      if(verbose) message("Error initialising SaemixModel object:\n   The covariance model needs to have the same size as the number of parameters.\n")
      return(.Object)
    }
    if(!valid_covariance.model(covariance.model, verbose)) return(.Object)
    if(is.null(colnames(covariance.model))) colnames(covariance.model)<-rownames(covariance.model)<-colnames(psi0)
    .Object@covariance.model<-covariance.model
    indx.omega<-which(diag(covariance.model)>0)
    .Object@indx.omega<-indx.omega
    if(!missing(omega.init) && length(omega.init)>0) {
      if(dim(omega.init)[1]!=dim(omega.init)[2]) {
        if(verbose) message("Warning:   the matrix giving the initial conditions for the covariance model (omega.init) needs to be a square matrix. Changing it to the diagonal matrix\n")
        omega.init<-NULL
      }
    }
    if(missing(omega.init) || length(omega.init)==0) {
      omega.init<-diag(fixed.estim)
      diag.omegi<-rep(1,npar)
      j1<-which(transform.par==0)
      if(length(j1)>0) {
        diag.omegi[j1]<-sapply(psi0[1,j1]**2,function(x) { x[x<1]<-1; return(x)})
#      for(i in j1) d[i]<-max(psi0[i]^2,1)
      }
      omega.init<-diag(diag.omegi,nrow=npar)
    }
    if(is.null(colnames(omega.init))) {
      if(dim(omega.init)[1]==length(colnames(psi0))) colnames(omega.init)<-rownames(omega.init)<-colnames(psi0) else message("The dimensions of omega.init don't agree with the number of parameters")
    }
    .Object@omega.init<-omega.init
    if(sum(.Object@fixed.estim*mydiag(.Object@covariance.model))==0) {
      message("Error initialising SaemixModel object:\n")
#     if(sum(mydiag(.Object@covariance.model))==0) cat("At least one parameter with IIV must be included in the model.\n") else cat("At least one parameter with IIV must be estimated and not fixed in the model.\n")
      return(.Object)
    }

## Residual Error model.
# error models are a + bf described by [a b]
# error models :
#   constant            y = f + a*e
#   proportional        y = f + b*f*e
#   combined            y = f + (a+b*f)*e
#   exponential         y = f*exp(a*e)    ( <=>  log(y) = log(f) + a*e )
    if(missing(error.init) || length(error.init)!=2*length(.Object@error.model)) {
      error.init<-c()
      for(i in 1:length(.Object@error.model)) {
        error.init<-c(error.init,switch(error.model[i],
        "constant"=c(1,0),
        "exponential"=c(1,0),
        "proportional"=c(0,1),
            "combined"=c(1,1)))
     }
    }
    xres<-c()
    if(missing(name.sigma)) mis.sig<-TRUE else mis.sig<-FALSE
    if(missing(name.sigma) || length(name.sigma)!=2) name.sigma<-c("a","b")
    if(!mis.sig) { # & .Object@name.response!="" # pb if response has more than 1 element
      for(i in 1:length(.Object@name.response)) {
        if(.Object@name.response[i]!="") xres<-c(xres,paste(name.sigma,.Object@name.response[i],sep=".")) else xres<-c(xres,paste(name.sigma,i,sep=".")) # xres<-c(xres,name.sigma) ?
      }
    } else xres<-rep(name.sigma,length(.Object@name.response))
    .Object@name.sigma<-xres
    .Object@error.init<-error.init
    indx.res<-c()
    for(i in 1:length(.Object@error.model)) {
      if(.Object@error.model[i]=='constant') {
        indx.res1<-1
    } else {
        if(.Object@error.model[i]=='proportional') {
          indx.res1<-2
      } else {
          if(.Object@error.model[i]=='combined') {
            indx.res1<-c(1,2) 
        } else {
            if(.Object@error.model[i]=='exponential') {
              indx.res1<-1
           }
        }
      }
    }
      indx.res<-c(indx.res,indx.res1+2*(i-1))
    }
    .Object@error.init[-indx.res]<-0
    .Object@indx.res<-indx.res
    .Object@betaest.model<-matrix(c(rep(1,.Object@nb.parameters), c(t(.Object@covariate.model))),ncol=.Object@nb.parameters,byrow=TRUE)
    colnames(.Object@betaest.model)<-colnames(.Object@covariate.model)
    if(!is.null(rownames(.Object@covariate.model))) {
      rownames(.Object@betaest.model)<-c("Fixed",rownames(.Object@covariate.model))
    } else {
      rownames(.Object@betaest.model)<-rep("",dim(.Object@betaest.model)[1])
      rownames(.Object@betaest.model)[1]<-"Fixed"
    }
# Object validation
    validObject(.Object)
    return (.Object)
  }
)

####################################################################################
####      saemixModel class - accesseur       ####
####################################################################################

#' Get/set methods for SaemixModel object
#' 
#' Access slots of a SaemixModel object using the object["slot"] format
#' 
#' @param x object
#' @param i element to be replaced
#' @param j element to replace with
#' @param drop whether to drop unused dimensions
#' @keywords methods
#' @exportMethod [
#' @exportMethod [<-


# Getteur
setMethod(
  f ="[",
  signature = "SaemixModel" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "model"={return(x@model)},
    "description"={return(x@description)},
    "modeltype"={return(x@modeltype)},
    "psi0"={return(x@psi0)},
    "transform.par"={return(x@transform.par)},
    "fixed.estim"={return(x@fixed.estim)},
    "error.model"={return(x@error.model)},
    "covariate.model"={return(x@covariate.model)},
    "betaest.model"={return(x@betaest.model)},
    "covariance.model"={return(x@covariance.model)},
    "omega.init"={return(x@omega.init)},
    "error.init"={return(x@error.init)},
    "nb.parameters"={return(x@nb.parameters)},
    "name.modpar"={return(x@name.modpar)},
    "name.fixed"={return(x@name.fixed)},
    "name.random"={return(x@name.random)},
    "name.sigma"={return(x@name.sigma)},
    "name.X"={return(x@name.X)},
    "name.response"={return(x@name.response)},
    "name.predictors"={return(x@name.predictors)},
    "name.cov"={return(x@name.cov)},
    "indx.fix"={return(x@indx.fix)},
    "indx.cov"={return(x@indx.cov)},
    "indx.omega"={return(x@indx.omega)},
    "indx.res"={return(x@indx.res)},
    "Mcovariates"={return(x@Mcovariates)},
    stop("No such attribute\n")
   )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixModel" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "model"={x@model<-value},
    "description"={return(x@description)},
    "modeltype"={return(x@modeltype)},
    "psi0"={x@psi0<-value},
    "transform.par"={x@transform.par<-value},
    "fixed.estim"={x@fixed.estim<-value},
    "error.model"={x@error.model<-value},
    "covariate.model"={x@covariate.model<-value},
    "betaest.model"={x@betaest.model<-value},
    "covariance.model"={x@covariance.model<-value},
    "omega.init"={x@omega.init<-value},
    "error.init"={x@error.init<-value},
    "nb.parameters"={x@nb.parameters<-value},
    "name.modpar"={x@name.modpar<-value},
    "name.fixed"={x@name.fixed<-value},
    "name.random"={x@name.random<-value},
    "name.sigma"={x@name.sigma<-value},
    "name.X"={x@name.X<-value},
    "name.response"={x@name.response<-value},
    "name.predictors"={x@name.predictors<-value},
    "name.cov"={x@name.cov<-value},
    "indx.fix"={x@indx.fix<-value},
    "indx.cov"={x@indx.cov<-value},
    "indx.omega"={x@indx.omega<-value},
    "indx.res"={x@indx.res<-value},
    "Mcovariates"={x@Mcovariates<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

####################################################################################
####      SaemixModel class - method to print/show data   ####
####################################################################################

#' @rdname print-methods
#' @exportMethod print

setMethod("print","SaemixModel",
  function(x,...) {
    cat("Nonlinear mixed-effects model\n")
    distrib<-c("normal","log-normal","probit","logit")
    cat("  Model function")
    if(length(x@description)>0 && nchar(x@description)>0) cat(": ",x@description)
    cat("  Model type")
    if(length(x@modeltype)>0 && nchar(x@modeltype)>0) cat(": ",x@modeltype)
    cat("\n")
    print(x@model)
    cat("  Nb of parameters:",x@nb.parameters,"\n")
    cat("      parameter names: ",x@name.modpar,"\n")
    cat("      distribution:\n")
    tab<-cbind(Parameter=x@name.modpar,Distribution=distrib[x@transform.par+1], Estimated=ifelse(x@fixed.estim==1,"Estimated","Fixed"))
    print(tab,quote=FALSE)
    cat("  Variance-covariance matrix:\n")
    tab<-x@covariance.model
#    try(colnames(tab)<-rownames(tab)<-x@name.modpar)
    print(tab,quote=FALSE)
    st1<-paste(x@name.sigma,x@error.init,sep="=")
    if (x@modeltype=="structural"){
      cat("  Error model:",x@error.model,", initial values:",st1[x@indx.res],"\n")
    }
   if(dim(x@covariate.model)[1]>0) {
      cat("  Covariate model:")
      if(sum(x@covariate.model)==0) cat(" none\n") else {
        cat("\n")
        print(x@covariate.model)
    }
  } else cat("    No covariate in the model.\n")
    cat("    Initial values\n")
    print(x@psi0)
  }
)

#' @rdname show-methods
#' @exportMethod show

setMethod("show","SaemixModel",
  function(object) {
    cat("Nonlinear mixed-effects model\n")
    cat("  Model function")
    if(length(object@description)>0 && nchar(object@description)>0) {
      cat(": ",object@description,"\n")}
    else {
      cat("\n")
      print(object@model)
    }
    fix1<-ifelse(object@fixed.estim==1,""," [fixed]")
    cat("    ",object@nb.parameters,"parameters:", paste(object@name.modpar,fix1,sep=""),"\n")
    cat("     error model:",object@error.model,"\n")
    if(dim(object@covariate.model)[1]>0) {
      cat("     covariate model:\n")
      print(object@covariate.model) 
    } else cat("No covariate\n")
  }
)

#' @rdname showall-methods
#' @exportMethod showall

setMethod("showall","SaemixModel",
  function(object) {
    cat("Nonlinear mixed-effects model\n")
    distrib<-c("normal","log-normal","probit","logit")
    cat("  Model function")
    if(length(object@description)>0 && nchar(object@description)>0) cat(": ",object@description)
    cat("\n")
    print(object@model)
    cat("  Nb of parameters:",object@nb.parameters,"\n")
    cat("      parameter names: ",object@name.modpar,"\n")
    if(length(object@name.fixed)>0) cat("      fixed parameters: ",object@name.fixed,"\n")
    if(length(object@name.random)>0) cat("      random parameters: ",object@name.random,"\n")
    if(length(object@name.sigma)>0) cat("      parameters of residual variability: ",object@name.sigma,"\n")
    if(length(object@name.predictors)>0) cat("      predictors: ",object@name.predictors,"\n")
    if(length(object@name.X)>0) cat("      X predictor: ",object@name.X,"\n")
    if(length(object@name.cov)>0) cat("      covariates: ",object@name.cov,"\n")
    cat("      distribution:\n")
    tab<-cbind(Parameter=object@name.modpar, Distribution=distrib[object@transform.par+1], Estimated=ifelse(object@fixed.estim==1, "Estimated","Fixed"))
    print(tab,quote=FALSE)
    cat("  Variance-covariance matrix:\n")
    tab<-object@covariance.model
    print(tab,quote=FALSE)
    cat("  Initial estimate for variance-covariance matrix:\n")
    print(object@omega.init)
    st1<-paste(object@name.sigma,object@error.init,sep="=")
    cat("  Error model:",object@error.model,", initial values:",st1[object@indx.res],"\n")
   if(dim(object@covariate.model)[1]>0) {
      cat("  Covariate model:")
      if(sum(object@covariate.model)==0) cat(" none\n") else {
        cat("\n")
        print(object@covariate.model)
    }
  } else cat("  No covariate in the model.\n")
    cat("    Initial values\n")
    print(object@psi0)
    if(length(object@indx.fix)>0) cat("      index for fixed parameters: ", object@indx.fix,"\n")
    if(length(object@indx.cov)>0) cat("      index for covariate parameters: ", object@indx.cov,"\n")
    if(length(object@indx.omega)>0) cat("      index for random parameters: ", object@indx.omega,"\n")
    if(length(object@indx.res)>0) cat("      index for parameters of residual variability: ", object@indx.res,"\n")
    if(length(object@Mcovariates)>0) print(object@Mcovariates)
  }
)


####################################################################################
####        Summary method for SaemixModel      ####
####################################################################################

#' @rdname summary-methods
#' @exportMethod summary

setMethod("summary","SaemixModel",
  function(object, print=TRUE, ...) {
    if(print) {
      cat("Nonlinear mixed-effects model\n")
    cat("  Model function")
    if(length(object@description)>0 && nchar(object@description)>0) {
      cat(": ",object@description,"\n")}
    else {
      cat("\n")
      print(object@model)
    }
    fix1<-ifelse(object@fixed.estim==1,""," [fixed]")
    cat("    ",object@nb.parameters,"parameters:", paste(object@name.modpar,fix1,sep=""),"\n")
    cat("     error model:",object@error.model,"\n")
    if(dim(object@covariate.model)[1]>0) {
      cat("     covariate model:\n")
      print(object@covariate.model) 
    } else cat("No covariate\n")
    }
    distrib<-c("normal","log-normal","probit","logit")
    tab.par<-data.frame(Parameter=object@name.modpar, Distribution=distrib[object@transform.par+1], Estimated=ifelse(as.numeric(object@betaest.model[1,])==1,"estimated","fixed"), Initial.value=object@psi0[1,])
    tab.res<-data.frame(parameters=object@name.sigma,Initial.value=object@error.init)   
    res<-list(model=list(modeltype=object@modeltype,model.function=object@model, error.model=object@error.model),parameters=list(fixed=tab.par, residual.error=tab.res),covariance.model=object@covariance.model, covariate.model=object@covariate.model)
    invisible(res)
 }
)

####################################################################################
####      SaemixModel class - method to plot      ####
####################################################################################

#' Plot model predictions using an SaemixModel object
#' 
#' This function will plot predictions obtained from an SaemixModel object over a given range of X. Additional predictors may be passed on to the function using the predictors argument.
#' 
#' @name plot-SaemixModel
#' 
#' @param x an SaemixData object or an SaemixSimData object
#' @param y unused, present for compatibility with base plot function
#' @param range range of X over which the model is to be plotted
#' @param psi parameters of the model 
#' @param predictors additional predictors needed to pass on to the model
#' @param ... additional arguments to be passed on to plot (titles, legends, ...)
#' 
#' @aliases plot,SaemixModel-methods 
#' @aliases plot,SaemixModel
#' @keywords methods
### #' @docType methods
#' @exportMethod plot
#' @rdname plot-SaemixModel

# Plot simulations from the model
# ECO TODO: adjust to multiple responses

setMethod("plot","SaemixModel",
  function(x,y,range=c(0,1),psi,predictors,...) {
    if(missing(psi)) psi<-x@psi0[1,]
    psi<-matrix(psi,nrow=1)
    npred<-length(x@name.predictors)
    if(npred==0 & missing(predictors)) npred<-1 else {
      if(npred==0 & !missing(predictors)) {
        npred<-1+length(predictors)
      } else {
        if(npred>1 & (missing(predictors) || length(predictors)<(npred-1))) {
        cat("Please provide the value of the predictors other than X\n")
        return()
      }
     }
    }
    if(length(x@name.response)>1) {
      cat("Currently the plot can only be obtained for single-response models.\n")
      return()
    }
    npts<-100
    psi<-matrix(rep(psi,npts+1),byrow=T,nrow=(npts+1))
    id<-matrix(rep(1,npts+1),ncol=1)
    xval<-range[1]+(range[2]-range[1])*c(0:100)/100
    if(npred==1) {
      xdep<-matrix(xval,ncol=1)
    } else {
      xdep<-cbind(xval,matrix(rep(predictors[1:(npred-1)],(npts+1)), byrow=T,nrow=(npts+1)))
      if(length(x@name.X)>0) {
        colnames(xdep)<-c(x@name.X,x@name.predictors[x@name.predictors!=x@name.X])
        xdep<-xdep[,match(x@name.predictors,colnames(xdep))]
      } else colnames(xdep)<-paste("Predictor",1:npred)
    }
    ypred<-try(x@model(psi,id,xdep))
    if(!is.numeric(ypred)) {
      cat("Problem when attempting to obtain predictions from the model.\n")
      cat("Usage: plot(x,range=c(0,1),psi,predictors) \n")
      cat("Possible solutions can be:\n")
      cat("   1. provide suitable values for X (option range=c(<lower bound>, <upper bound>))\n")
      cat("   2. provide values for additional predictors (option predictors=c(<value for predictor 1>, <value for predictor 2>, ...)).\n")
      cat("   3. check values for the model parameters (defaults to component psi0[1,] of the model).\n")
      cat("   4. the predictor used the X-axis is assumed to be in the first column; please check your model is written in a compatible way.\n")
    } else {
      if(length(x@name.X)==0 | length(x@name.predictors)==0) cat("Warning: X predictor supposed to be on the first axis\n")
      cat("Plot characteristics:\n")
      if(npred>1) {
        for(j in 1:dim(xdep)[2]) {
    if(length(x@name.X)==0) {
      if(j>1) cat("   predictor:",colnames(xdep)[j],"=",xdep[1,j],"\n")
    } else {
      if(colnames(xdep)[j]!=x@name.X) cat("    predictor:",colnames(xdep)[j],"=",xdep[1,j],"\n")
    }
      }}
      cat("   range for X-axis:",min(xval),"-",max(xval),"\n")
      cat("   parameters used in the simulation:", paste(x@name.modpar,"=",psi[1,],collapse=", "),"\n")
      plot(xval,ypred,type="l",xlab=ifelse(length(x@name.X)==0, "X",x@name.X),ylab=ifelse(length(x@name.response)==0, "Response",x@name.response))
    }
  }
)

####################################################################################
####      SaemixModel class - User-level function     ####
####################################################################################

#' Function to create a SaemixModel object
#' 
#' This function creates a SaemixModel object. The two mandatory arguments are
#' the name of a R function computing the model in the SAEMIX format (see
#' details and examples) and a matrix psi0 giving the initial estimates of the
#' fixed parameters in the model, with one row for the population mean
#' parameters and one row for the covariate effects (see documentation).
#' 
#' This function is the user-friendly constructor for the SaemixModel object
#' class.
#' 
#' @param model name of the function used to compute the structural model. The
#' function should return a vector of predicted values given a matrix of
#' individual parameters, a vector of indices specifying which records belong
#' to a given individual, and a matrix of dependent variables (see example
#' below).
#' @param psi0 a matrix with a number of columns equal to the number of
#' parameters in the model, and one (when no covariates are available) or two
#' (when covariates enter the model) giving the initial estimates for the fixed
#' effects. The column names of the matrix should be the names of the
#' parameters in the model, and will be used in the plots and the summaries.
#' When only the estimates of the mean parameters are given, psi0 may be a
#' named vector.
#' @param description a character string, giving a brief description of the
#' model or the analysis
#' @param modeltype a character string, giving model type (structural or likelihood)
#' @param name.response the name of the dependent variable
#' @param name.sigma a vector of character string giving the names of the residual error parameters
#' @param error.model type of residual error model (valid types are constant,
#' proportional, combined and exponential). Defaults to constant
#' @param transform.par the distribution for each parameter (0=normal,
#' 1=log-normal, 2=probit, 3=logit). Defaults to a vector of 1s (all parameters
#' have a log-normal distribution)
#' @param fixed.estim whether parameters should be estimated (1) or fixed to
#' their initial estimate (0). Defaults to a vector of 1s
#' @param covariate.model a matrix giving the covariate model. Defaults to no
#' covariate in the model
#' @param covariance.model a square matrix of size equal to the number of
#' parameters in the model, giving the variance-covariance matrix of the model:
#' 1s correspond to estimated variances (in the diagonal) or covariances
#' (off-diagonal elements). Defaults to the identity matrix
#' @param omega.init a square matrix of size equal to the number of parameters
#' in the model, giving the initial estimate for the variance-covariance matrix
#' of the model. Defaults to the identity matrix
#' @param error.init a vector of size 2 giving the initial value of a and b in
#' the error model. Defaults to 1 for each estimated parameter in the error
#' model
#' @param name.modpar names of the model parameters, if they are not given as
#' the column names (or names) of psi0
#' @param verbose a boolean, controlling whether information about the created should be printed out. Defaults to TRUE
#' @return A SaemixModel object (see \code{\link{saemixModel}}).
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{saemixControl}},\code{\link{saemix}}
#' @references Comets  E, Lavenu A, Lavielle M. Parameter estimation in nonlinear mixed effect models using saemix, an R implementation of the SAEM algorithm. Journal of Statistical Software 80, 3 (2017), 1-41.
#' 
#' Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.
#' 
#' Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm.
#' 20th meeting of the Population Approach Group in Europe, Athens, Greece
#' (2011), Abstr 2173.
#' @keywords models
#' @examples
#' 
#' model1cpt<-function(psi,id,xidep) { 
#'    dose<-xidep[,1]
#'    tim<-xidep[,2]  
#'    ka<-psi[id,1]
#'    V<-psi[id,2]
#'    CL<-psi[id,3]
#'    k<-CL/V
#'    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#'    return(ypred)
#' }
#' 
#' saemix.model<-saemixModel(model=model1cpt,
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' @export saemixModel

saemixModel<-function(model,psi0,description="",modeltype ="structural", name.response="", name.sigma=character(), error.model=character(), transform.par=numeric(),fixed.estim=numeric(),covariate.model=matrix(nrow=0,ncol=0), covariance.model=matrix(nrow=0,ncol=0),omega.init=matrix(nrow=0,ncol=0),error.init=numeric(), name.modpar=character(), verbose=TRUE) {
# Creating model from class
  if(missing(model)) {
    if(verbose) cat("Error in saemixModel:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
    return("Creation of SaemixModel failed")  
  }
  xcal<-try(typeof(model))
  if(inherits(xcal,"try-error")) {
    if(verbose) cat("Error in saemixModel:\n   the model function does not exist.\n")
    return("Creation of SaemixModel failed")  
  }
  if(typeof(model)=="character") {
    if(exists(model)) model<-get(model) else {
      if(verbose) cat("Error in saemixModel:\n   The argument model to saemixModel must be a valid function.\n")
      return("Creation of SaemixModel failed")
    }
  }
  if(!is.function(model)) {
    if(verbose) cat("Error in saemixModel:\n   The argument model to saemixModel must be a valid function.\n")
    return("Creation of SaemixModel failed")
  }
  if(length(formals(model))!=3) {
    if(verbose) cat("Error in saemixModel:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
    return("Creation of SaemixModel failed")
  }
  if(missing(psi0) || length(psi0)==0) {
    if(verbose) cat("Error in saemixModel:\n   please provide initial estimates psi0 for at least the fixed effects.\n")
    return("Creation of SaemixModel failed")  
  }
  if(is.null(dim(psi0))) {
    psi1<-matrix(psi0,nrow=1)
    if(!is.null(names(psi0))) colnames(psi1)<-names(psi0)
    psi0<-psi1
    if(verbose) cat("Warning: psi0 given as a vector, reshaping it to a matrix.\n")
  }
  if(is.null(colnames(psi0))) {
    if(verbose) cat("Warning: no names given for the parameters in the model, please consider including parameter names.\n")
  }
  xmod<-try(new(Class="SaemixModel",model=model,description=description , modeltype=modeltype,psi0=psi0, name.response=name.response, name.sigma=name.sigma, error.model=error.model, transform.par=transform.par,fixed.estim=fixed.estim, covariate.model=covariate.model,covariance.model=covariance.model, omega.init=omega.init,error.init=error.init,name.modpar=name.modpar))
  if(class(xmod)=="SaemixModel") x1<-try(validObject(xmod),silent=FALSE) else x1<-xmod
  if(!inherits(x1,"try-error")) {
    if(verbose) cat("\n\nThe following SaemixModel object was successfully created:\n\n")
    } else xmod<-"Creation of SaemixModel failed"
  if(verbose) print(xmod)
  return(xmod)
}

####################################################################################
####  Auxiliary function

#' Matrix diagonal
#' 
#' Extract or replace the diagonal of a matrix, or construct a diagonal matrix (replace diag function from R-base)
#' 
#' @param x a matrix, vector or 1D array, or missing.
#' @param nrow, ncol Optional dimensions for the result when x is not a matrix. 
#' @param value either a single value or a vector of length equal to that of the current diagonal. Should be of a mode which can be coerced to that of x.
#' @return If x is a matrix then diag(x) returns the diagonal of x. The resulting vector will have names if the matrix x has matching column and rownames.
#' @seealso \code{diag}
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @keywords models
#' @examples
#' 
#' mydiag(1)
#' mydiag(c(1,2))
#' 
#' @export mydiag

# Redefining diag function, too many problems with the R version
mydiag <- function (x = 1, nrow, ncol) {
  if (is.matrix(x)) {
    if (nargs() > 1L) 
      stop("'nrow' or 'ncol' cannot be specified when 'x' is a matrix")
    if ((m <- min(dim(x))) == 0L) 
      return(vector(typeof(x), 0L))
    y <- c(x)[1L + 0L:(m - 1L) * (dim(x)[1L] + 1L)]
    nms <- dimnames(x)
    if (is.list(nms) && !any(sapply(nms, is.null)) && identical((nm <- nms[[1L]][seq_len(m)]), 
                                                                nms[[2L]][seq_len(m)])) 
      names(y) <- nm
    return(y)
  }
  if (is.array(x) && length(dim(x)) != 1L) 
    stop("'x' is an array, but not 1D.")
  if (missing(x)) 
    n <- nrow
  else n <- length(x)
  if (!missing(nrow)) 
    n <- nrow
  if (missing(ncol)) 
    ncol <- n
  p <- ncol
  y <- array(0, c(n, p))
  if ((m <- min(n, p)) > 0L) 
    y[1L + 0L:(m - 1L) * (n + 1L)] <- x
  y
}

############ VALIDITY OF COVARIANCE MODEL 
valid_covariance.model <- function(covariance.model, verbose=TRUE){
  cm <- covariance.model
  
  #non-square matrix
  if(dim(cm)[1]!=dim(cm)[2]) {
    if(verbose) message("Error initialising SaemixModel object:\n   The covariance model needs to be a square matrix, please check dimensions.\n")
    return(FALSE)
  }
  
  
  #values other than 1 or 0
  s <- sum(cm[cm!=1 & cm!=0])
  if (s>0){
    if(verbose) message("Error initialising SaemixModel object:\n  Invalid covariance model, only 0 or 1 values accepted, please change covariance model.\n")
    return(FALSE)
  }
  
  #asymmetrical
  if (!all(t(cm)==cm)){
    if(verbose) message("Error initialising SaemixModel object:\n  The matrix of covariance model is not symmetrical, please change covariance model.\n")
    return(FALSE)
  }
  
  #values other than 0 when diagonal number is 0
  for (i in 1:nrow(cm)){
    for(j in 1:ncol(cm)){
      if (cm[i,j]!=0){
        if(cm[i,i]==0 |cm[j,j]==0){
          if(verbose) message("Error initialising SaemixModel object:\n  The covariance model is invalid, please change covariance model.\n")
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)
}




