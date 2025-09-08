## ToDo:
### Create a legacy constructor function with the same input as previously (same arguments)
### but output=object with new SaemixModel class

####################################################################################
####			SaemixModel class - User-level function			####
####################################################################################

#' Function to create an SaemixModel object
#' 
#' This function creates an SaemixModel object. The two mandatory arguments are
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
#' @param simulate.function for non-Gaussian data models, defined as modeltype='likelihood', 
#' the name of the function used to simulate from the structural model. The
#' function should have the same header as the model function, and should return 
#' a vector of simulated values given a matrix of individual parameters, 
#' a vector of indices specifying which records belong to a given individual, 
#' and a matrix of dependent variables (see example in the documentation, section
#' discrete data examples)
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
#' of the model.
#' @param error.init a vector of size 2 giving the initial value of a and b in
#' the error model. Defaults to 1 for each estimated parameter in the error
#' model
#' @param name.modpar names of the model parameters, if they are not given as
#' the column names (or names) of psi0
#' @param verbose a boolean, controlling whether information about the created should be printed out. Defaults to TRUE
#' @return An SaemixModel object (see \code{\link{saemixModel}}).
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{saemixControl}},\code{\link{saemix}}
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
#' @keywords models
#' @examples
#' 
#' model1cpt<-function(psi,id,xidep) { 
#' 	  dose<-xidep[,1]
#' 	  tim<-xidep[,2]  
#' 	  ka<-psi[id,1]
#' 	  V<-psi[id,2]
#' 	  CL<-psi[id,3]
#' 	  k<-CL/V
#' 	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#' 	  return(ypred)
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

saemixModel.legacy<-function(model,psi0,description="",modeltype ="structural", name.response="", name.sigma=character(), error.model=character(), transform.par=numeric(),fixed.estim=numeric(),covariate.model=matrix(nrow=0,ncol=0), covariance.model=matrix(nrow=0,ncol=0),omega.init=matrix(nrow=0,ncol=0),error.init=numeric(), name.modpar=character(), simulate.function=NULL, verbose=TRUE) {
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
  if(is(model,"character")) {
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
  has.sim<-FALSE
  if(!is.null(simulate.function)) {
    xcal<-try(typeof(simulate.function))
    if(inherits(xcal,"try-error")) {
      if(verbose) message("The simulate.function does not exist, ignoring.\n")
    } else {
      if(is(simulate.function,"character")) {
        if(exists(simulate.function)) {
          simulate.function<-get(simulate.function)
          has.sim <- TRUE
        }
      }
      if(!is.function(simulate.function) || length(formals(simulate.function))!=3) {
        if(verbose) cat("The simulate.function should have the same format as the model function, ignoring.\n")
        has.sim <- FALSE
      } else has.sim <- TRUE
    }
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
  if(is(xmod,"SaemixModel")) x1<-try(validObject(xmod),silent=FALSE) else x1<-xmod
  if(!inherits(x1,"try-error")) {
    if(verbose) cat("\n\nThe following SaemixModel object was successfully created:\n\n")
  } else xmod<-"Creation of SaemixModel failed"
  if(has.sim) xmod@simulate.function <- simulate.function
  if(verbose) print(xmod)
  return(xmod)
}

