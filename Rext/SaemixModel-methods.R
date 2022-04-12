######################################################################################
####			S4 methods: show, print, showall	
######################################################################################

######################################################################################
####			S4 methods: plot	
######################################################################################


######################################################################################
####			SaemixModel class - User-level S3 function to create an object from the class
######################################################################################

#' Function to create a SaemixModel object
#' 
#' This function is the user-friendly constructor for the SaemixModel object class.
#' 
#' Two elements need to be passed to saemixModel():
#' - model: the name of a R function computing the model in SAEMIX format, with a set header 
#' (see details and examples) 
#' - either npar (the number of parameters in the model) or parameter (a named vector of starting
#' values)
#'   * if npar is given and parameter is omitted, the parameters will be named theta1, theta2, ... 
#'   and their starting value will be set to 1
#'   * if parameter is given as a list of values without names, the parameters will be named 
#'   theta1, theta2, ... and the values will be used as starting values
#'   * if parameter is given as a named vector, the names will be used to name the parameters
#'   and the values will be used as starting values
#' If the model contains more than one response, these nned to be specified using the new outcome argument.
#' 
#' All other arguments are optional (see defaults in the list below) but can be set to create a full 
#' model (see examples)
#' 
#' Note that since version 3.0 the call to saemixModel() has changed considerably. A legacy
#' function, saemixModel.legacy(), is available to avoid breaking previous code, but we 
#' encourage users to take advantage of the new definition of models which is much more 
#' structured, consistent and flexible.
#' 
#' @param model name of the function used to compute the structural model. The
#' function should return a vector of predicted values given a matrix of
#' individual parameters, a vector of indices specifying which records belong
#' to a given individual, and a matrix of dependent variables (see example
#' below).
#' @param npar the number of parameters in the model (if only npar is given, the parameters
#' will be named theta1, theta2,..., thetanpar and initialised to 1)
#' @param parameter a list or vector of named values, which will be used in the plots and the summaries.
#' Will evolve in future versions to a list of SaemixParameter objects containing the information
#' about parameter distribution and covariate models
#' @param description a character string, giving a brief description of the
#' model or the analysis (for user reference only)
#' @param outcome either a vector giving the name of the responses, or a list of (named) outcomes
#' defined through the functions continuousOutcome() and discreteOutcome() 
#' (see \code{\link{continuousOutcome}} and \code{\link{discreteOutcome}}).
#' If given as a vector of names, each response will then be assumed to be
#' continuous responses with a constant residual error model.
#' If given as a list of outcomes, additional elements can be specified for each response, 
#' such as residual error model and starting values for continuous responses, 
#' or the type of response for non-continuous responses.
#' @param transform.par the distribution for each parameter (0=normal,
#' 1=log-normal, 2=probit, 3=logit). Defaults to a vector of 1s (all parameters
#' have a log-normal distribution)
#' @param mu.fix set to 1 for parameters which are fixed (to their starting value) 
#' during the estimation (defaults to 0)
#' @param covariate.model a matrix giving the covariate model. Defaults to no
#' covariate in the model.
#' @param covariate.model.fix a matrix giving the covariate model. Defaults to empty
#' @param beta.start starting values for the covariate effect parameters (corresponding to 1's
#' in the matrix covariate.model)
#' @param covariance.model a square matrix of size equal to the number of
#' parameters in the model, giving the variance-covariance matrix of the model:
#' 1s correspond to estimated variances (in the diagonal) or covariances
#' (off-diagonal elements). Defaults to the identity matrix.
#' @param covariance.model.fix a matrix of the same size as covariance.model with 1's for 
#' parameters fixed during the estimation (defaults to a 0x0 matrix)
#' @param covariance.model.start a square matrix of size equal to the number of parameters
#' in the model, giving the initial estimate for the variance-covariance matrix
#' of the model (was omega.init)
#' @param parameter names of the model parameters, if they are not given as
#' the column names (or names) of psi0
#' @param verbose a boolean, controlling whether information about the created should be printed out. 
#' Defaults to FALSE (changed from version 3.0)
#' 
#' @return A SaemixModel object (see \code{\link{saemixModel}}).
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{saemixControl}},\code{\link{saemix}}
#' 
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
#'   outcome=continuousOutcome(model="proportional", start=0.5),
#'   description="One-compartment model with first-order absorption", 
#'   parameter=c(ka=1,V=20,CL=0.5),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE), beta.start=,
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE))
#' 
#' @export saemixModel

# covariate.model=matrix(nrow=0,ncol=0), covariate.model.fix=matrix(nrow=0,ncol=0), beta.start=numeric(),
# covariance.model=matrix(nrow=0,ncol=0), covariance.model.fix=matrix(nrow=0,ncol=0), 
# covariance.model.start=matrix(nrow=0,ncol=0),

# Eco 12/04/2022 - integrate these tests into the creator function, to keep the initialize of the classes as simple as can be
if(!missing(mu.start) & !(missing(parameter))) {
  if(!is.null(dim(mu.start))) {
    mu.start<-mu.start[1,] 
    if(missing(beta.start)) beta.start<-c(t(mu.start[-c(1),]))
  } else {
    mu.start<-mu.start
    if(missing(beta.start)) beta.start<-numeric()
  }
  if(length(mu.start)<length(parameter)) {
    if(verbose) message("mu.start argument should have the same length as the parameter vector, the starting value of additional parameters will be set to 1")
    mu.start<-c(mu.start, rep(1,length(parameter)-length(mu.start)))
  }
  mu.start<-mu.start[1:length(parameter)] # if mu.start longer than the number of parameters
  npar<-length(mu.start)
}
if(missing(parameter)) { # if mu.start is a named vector, use the names for the parameters (legacy)
  if(is.null(dim(mu.start))) { # mu.start is a vector
    npar<-length(mu.start)
    if(is.null(names(mu.start))) parameter<-paste0("theta",1:npar) else parameter<-names(mu.start)
    mu.start<-mu.start
    if(missing(beta.start)) beta.start<-numeric()
  } else {
    npar<-dim(mu.start)[2] # mu.start is an array (legacy)
    if(is.null(colnames(mu.start))) parameter<-paste0("theta",1:npar) else parameter<-colnames(mu.start)
    if(missing(beta.start)) beta.start<-c(t(mu.start[-c(1),]))
    mu.start<-mu.start[1,]
  }
} else {
  if(is.numeric(parameter)) {
    if(missing(mu.start)) {
      mu.start<-unname(parameter)
    }
    npar<-length(parameter)
    parameter<-names(parameter)
    if(is.null(parameter)) parameter<-paste0("theta",1:npar)
  } # here add a test if we decide to give parameter a more structured class
}
if(missing(mu.start)) {# if the name of the parameters are given but not their initial value, initialise to 1
  npar<-length(parameter)
  mu.start<-mu.start<-rep(1,npar)
  beta.start<-numeric()
}
if(missing(transform.par)) transform.par<-rep(0, npar)
if(length(transform.par)!=npar) transform.par<-rep(transform.par, length.out=npar)
.Object@transform.par<-transform.par
if(!missing(mu.fix)) {
  if(length(mu.fix)==npar) .Object@mu.fix<-mu.fix else 
    if(verbose) message("Size mismatch between mu.fix and the number of parameters, ignoring")
}


saemixModel<-function(model, parameter, outcome, var.model, description="",
                      npar=0, transform.par=numeric(), mu.fix=numeric(), 
                      covariate.model=NULL, covariate.model.fix=NULL, beta.start=numeric(),
                      covariance.model=NULL, covariance.model.fix=NULL, covariance.model.start=NULL,
                      verbose=FALSE) {
  # Testing whether model is given
  if(missing(model)) {
    if(verbose) message("Error in saemixModel:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
    return("Creation of SaemixModel failed")  
  }
  # Parameters
  misnpar<-(npar==0)
  miscipar<-missing(parameter)
  if(misnpar & miscipar) { # need one of npar (nb of parameters), mu.start or parameter (defining parameters)
    if(verbose) message("Error in saemixModel:\n Please give either the number of parameters in the model or the parameters with their starting values.\n")
    return("Creation of SaemixModel failed")  
  } else {
    if(miscipar) { # npar given
      mu.start<-rep(1,npar)
    } else { # parameter given
      npar<-length(parameter)
      if(is(parameter,"numeric")) {
        mu.start<-unname(parameter)
      } else {
        if(is(parameter,"character")) { #only names given
          mu.start<-rep(1,npar)
          names(parameter)<-parameter
        } else {
          if(verbose) message("Error in saemixModel:\n The argument parameter should be a named vector with the starting values of the parameters.\n")
          return("Creation of SaemixModel failed")  
        }
      }
    }
  }
  if(is.null(names(parameter))) parnames<-paste0("theta",1:npar) else parnames<-names(parameter)
  
  # Outcome - check format and create SaemixOutcome model objects
  if(!missing(outcome)) print(outcome)
  if(!missing(outcome)) {
    if(is(outcome, "character")) { # only the names of the responses => assume they are all continuous
      smx.out<-vector(length(outcome), mode="list")
      for(i in 1:length(outcome)) smx.out[[i]]<-createSaemixOutcome(continuousOutcome())
      names(smx.out)<-outcome
    } else {
      if(is(outcome, "list")) {
        smx.out<-vector(length(outcome), mode="list")
        for(i in 1:length(outcome)) {
          if(is(outcome[[i]], "SaemixOutcome")) smx.out[[i]]<-outcome[[i]] else { # lists created by the creator functions
            x1<-try(createSaemixOutcome(outcome[[i]]))
            if(is(x1,"try-error")) {
              if(verbose) message("Error in saemixModel:\n Outcome must either be a vector of (continuous) response names, a list of outcomes defined through continuousOutcome() and discreteOutcome(), or a list of outcomes defined through createSaemixOutcome().\n")
              return("Creation of SaemixModel failed") 
            } else smx.out[[i]]<-x1
          } 
        }
      } else {
        if(verbose) message("Error in saemixModel:\n Outcome must either be a vector of (continuous) response names, a list of outcomes defined through continuousOutcome() and discreteOutcome(), or a list of outcomes defined through createSaemixOutcome().\n")
        return("Creation of SaemixModel failed") 
      }
    }
  } else {
    outcome<-list(createSaemixOutcome(continuousOutcome()))
    names(outcome)<-outcome[[1]]@name
  }

  # Variability model
  if(missing(var.model)) {
    var.model<-saemixVarModel(size=npar, omega=covariance.model.start, omega.model=covariance.model, omega.model.fix=covariance.model.fix)
  } else {
    if(is(var.model,"list")) {
      if(is(var.model[[1]],"SaemixVarLevel")) { # Will need to add a check on all elements when more levels of variability are used in the model
        if(verbose & length(var.model)>1) message("Currently only one level of variability is used in the algorithm, setting variability structure to the first element of var.model")
        var.model<-var.model[[1]]
      } else {
        x1<-try(saemixVarLevel(var.model[[1]]))
        if(is(x1,"try-error")) {
          if(verbose) cat("Error in saemixModel:\n var.model must either be a list of SaemixVarLevel objects defined using saemixVarModel, or a list containing the elements to create it.\n")
          return("Creation of SaemixModel failed") 
        } else var.model<-x1
      }
    }
  }
  x<-try(new(Class="SaemixModel", model=model, parameter=parnames,
         outcome=outcome, mu.start=mu.start, transform.par=transform.par, mu.fix=mu.fix,
         covariate.model=covariate.model, covariate.model.fix=covariate.model.fix, beta.start=beta.start,
         var.model=var.model, verbose=verbose))
  return(x)
}

#' Function to create a SaemixModel object, legacy version
#' 
#' This function creates a SaemixModel object using the same arguments as in the previous version of the package.
#' We advise users to move to the new way of specifying models
#' 
#' The two mandatory arguments to the legacy function are
#' the name of a R function computing the model in the SAEMIX format (see
#' details and examples) and a matrix psi0 giving the initial estimates of the
#' fixed parameters in the model, with one row for the population mean
#' parameters and one row for the covariate effects (see documentation).
#' 

# Legacy arguments to be deprecated in future versions
## modeltype
## psi0: preferred way will be through parameter
## name.sigma, error.model

saemixModel.legacy<-function(model, psi0, description="", modeltype ="structural", name.response="", 
                             name.sigma=character(), error.model=character(), name.modpar=character(), 
                             transform.par=numeric(), fixed.estim=numeric(), covariate.model=matrix(nrow=0,ncol=0),
                             covariance.model=matrix(nrow=0,ncol=0), omega.init=matrix(nrow=0,ncol=0),error.init=numeric(), 
                             verbose=FALSE) {
  # Testing input
  if(missing(model)) {
    if(verbose) cat("Error in saemixModel:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
    return("Creation of SaemixModel failed")  
  }
  misnpar<-missing(npar)
  mispsi0<-missing(psi0)
  miscipar<-missing(parameter)
  if(misnpar & miscipar & mispsi0) { # need one of npar (nb of parameters), psi0 or parameter (defining parameters)
    if(verbose) cat("Error in saemixModel:\n Please give either the number of parameters in the model or the parameters with their initial values.\n")
    return("Creation of SaemixModel failed")  
  }    
  if(!(misnpar) & miscipar & mispsi0) {
    psi0<-rep(1,npar)
    parameter<-paste0("theta",1:npar)
  }
    
    if(missing(parameter) & missing(npar)) {
      if(missing(psi0) || length(psi0)==0) {
        if(verbose) cat("Error in saemixModel:\n Please give the.\n")
        return("Creation of SaemixModel failed")  
        
      }
      
    }
  
  # Creating model from class
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
  xmod<-try(new(Class="SaemixModel", model=model, description=description, outcome=outcome, psi0=psi0, 
                name.modpar=name.modpar, transform.par=transform.par, fixed.estim=fixed.estim, covariate.model=covariate.model,
                iiv.model=iiv.model, iiv.init=iiv.init))
  if(class(xmod)=="SaemixModel") x1<-try(validObject(xmod),silent=FALSE) else x1<-xmod
  
}
