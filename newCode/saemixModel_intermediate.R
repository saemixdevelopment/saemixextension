
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
#' model (see examples). 
#' 
#' Note that since version 3.0 the call to saemixModel() has changed considerably. A legacy
#' function, \code{\link{saemixModel.legacy()}}, is available to avoid breaking previous code, but we 
#' encourage users to take advantage of the new definition of models which is much more 
#' structured, consistent and flexible.
#' 
#' The structure of the variability model can be specified in one of two ways:
#' - using the new object SaemixVarLevel, which defines a variability structure for the 
#' random effects, and passing it using the var.model argument
#' - using the arguments covariance.model (model for the variance-covariance matrix), covariance.model.start 
#' (the starting value of the matrix, previously omega.init) and covariance.model.fix (new feature, a matrix
#' specifying which variance-covariance parameters are to be fixed in the estimation).
#' If the var.model argument is missing or is invalid, the program will look for the other
#' arguments, but if the var.model is given and valid then the other arguments will be ignored.
#' 
#' @aliases sanitiseSaemixVarmodel
#' 
#' @param model name of the function used to compute the structural model. The
#' function should return a vector of predicted values given a matrix of
#' individual parameters, a vector of indices specifying which records belong
#' to a given individual, and a matrix of dependent variables (see example
#' below).
#' @param parameter a list or vector of named values, which will be used in the plots and the summaries.
#' Will evolve in future versions to a list of SaemixParameter objects containing the information
#' about parameter distribution and covariate models
#' @param npar the number of parameters in the model (if only npar is given, the parameters
#' will be named theta1, theta2,..., thetanpar and initialised to 1)
#' @param outcome either a vector giving the name of the responses, or a list of (named) outcomes
#' defined through the functions continuousOutcome() and discreteOutcome() 
#' (see \code{\link{continuousOutcome}} and \code{\link{discreteOutcome}}).
#' If given as a vector of names, each response will then be assumed to be
#' continuous responses with a constant residual error model.
#' If given as a list of outcomes, additional elements can be specified for each response, 
#' such as residual error model and starting values for continuous responses, 
#' or the type of response for non-continuous responses.
#' @param var.model list of SaemixVarLevel objects defining the structure of the variance-covariance matrix
#' for the different variability levels (currently only the first element of the list will be used, and 
#' will be considered to be attached to the grouping variable)
#' @param description a character string, giving a brief description of the
#' model or the analysis (for user reference only)
#' @param transform.par the distribution for each parameter (0=normal,
#' 1=log-normal, 2=probit, 3=logit). Defaults to a vector of 1s (all parameters
#' have a log-normal distribution)
#' @param mu.fix set to 1 for parameters which are fixed (to their starting value) 
#' during the estimation (defaults to rep(0,npar), all parameters estimated)
#' @param covariate.model a matrix giving the covariate model. Defaults to no
#' covariate in the model.
#' @param covariate.model.fix a matrix with the same dimensions as covariate.model indicating
#' which elements of the covariate model should be fixed in the estimation. Defaults to empty (all 
#' parameters estimated)
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
#' @param verbose a boolean, controlling whether information about the created should be printed out. 
#' Defaults to FALSE (changed from version 3.0)
#' 
#' @details 
#' If not empty, the outcome element can be given in different formats
#' - a number nout: will create nout continuous responses with a constant residual error model
#'   example: outcome=3 will create 3 such responses
#' - a vector of valid types: will create the corresponding responses with names y1, y2,...
#' example: outcome=c("continuous","binary","event") will create 3 responses (first: continuous; second: binary; third: event)
#' - a named list of types: will create the corresponding responses with their name
#' example: outcome=c(conc="continuous",response="binary",death="event") will create the same responses with their names
#' - a list of named elements resulting of calls to continuousOutcome() and discreteOutcome() allowing to pass additional arguments (see \code{\link{continuousOutcome}} and \code{\link{discreteOutcome}})
#' example: outcome=list(conc=continuousOutcome(model="combined2",start=c(2,0.5)),effect=continuousOutcome(),response=discreteOutcome(type="categorical", levels=c(0:4))) creates 3 responses: conc: a continuous response, with a combined2 error model and starting values of 2 and 0.5; effect: a continuous response with a constant error model and default starting values; response: a categorical response with 5 categories from 0 to 4
#' - a SaemixOutcome or a list of SaemixOutcome created by createSaemixOutcome [useful if outcomes have been created as objects before]
#' example:
#'   out1<-createSaemixOutcome(continuousOutcome())
#'   out2<-createSaemixOutcome(discreteOutcome())
#'   outcome=list(conc=out1, response=out2)
#' The last format is used to store the list of outcomes in the model object returned by saemixModel().
#' Note: No check is made with respect to the model itself as to whether the outcomes are actually compatible with the model.
#' 
#' @return A SaemixModel object (see \code{\link{saemixModel}}).
#' 
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


#' The structure of the variability model can be specified in one of two ways:
#' - using the new object SaemixVarLevel, which defines a variability structure for the 
#' random effects, and passing it using the var.model argument
#' - using the arguments covariance.model (model for the variance-covariance matrix), covariance.model.start 
#' (the starting value of the matrix, previously omega.init) and covariance.model.fix (new feature, a matrix
#' specifying which variance-covariance parameters are to be fixed in the estimation).
#' If the var.model argument is missing or is invalid, the program will look for the other
#' arguments, but if the var.model is given and valid then the other arguments will be ignored.


saemixModel.intermediate<-function(model, parameter, outcome, var.model, description="",
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
        if(is(parameter,"character")) { # only names given
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
  #  if(!missing(outcome)) print(outcome)
  if(missing(outcome)) {
    outcome<-list(createSaemixOutcome(continuousOutcome()))
    names(outcome)<-outcome[[1]]@name
  } else {
    outcome<-sanitiseSaemixOutcome(outcome=outcome)
    if(is(outcome,"character")) return(outcome)
  }
  # Variability model
  if(missing(var.model)) {
    var.model<-saemixVarModel(size=npar, omega=covariance.model.start, omega.model=covariance.model, omega.model.fix=covariance.model.fix)
  } else {
    var.model<-sanitiseSaemixVarmodel(var.model=var.model)
    if(is(var.model,"character")) return(var.model)
  }
  x<-try(new(Class="SaemixModel", model=model, parameter=parnames,
             outcome=outcome, mu.start=mu.start, transform.par=transform.par, mu.fix=mu.fix,
             covariate.model=covariate.model, covariate.model.fix=covariate.model.fix, beta.start=beta.start,
             var.model=var.model, verbose=verbose))
  return(x)
}


sanitiseSaemixVarmodel <- function(var.model) {
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
  return(var.model)
}

# Sanitise outcomes, replaced
# Function to validate and format the list of outcomes of the model (internal)

sanitiseSaemixOutcome <- function(outcome) {
  if(is(outcome,"numeric")) { # only a number of responses => assume they are all continuous
    smx.out<-vector(outcome[1], mode="list")
    for(i in 1:outcome[1]) smx.out[[i]]<-createSaemixOutcome(continuousOutcome(), name=paste0("y",i))
    names(smx.out)<-paste0("y",1:outcome[1])
    return(smx.out)
  }
  if(is(outcome, "character")) { # only the names of the responses => assume they are all continuous
    smx.out<-vector(length(outcome), mode="list")
    if(is.null(names(outcome)) & sum(is.na(match(outcome, c("continuous", "binary", "categorical","event"))))==0) 
      names(outcome)<-outcome
    if(is.null(names(outcome))) {
      for(i in 1:length(outcome)) 
        smx.out[[i]]<-createSaemixOutcome(continuousOutcome(), outcome[i])
      names(smx.out)<-outcome
    } else  {
      for(i in 1:length(outcome)) {
        if(!(outcome[i] %in% c("continuous", "binary", "count", "categorical","event"))) {
          if(verbose) message(paste("Unknown outcome type",outcome[i],", changing to continuous"))
          outcome[i]<-"continuous"
        }
        if(outcome[i]=="continuous") smx.out[[i]]<-createSaemixOutcome(continuousOutcome(), name=names(outcome)[i])
        if(outcome[i]=="binary") smx.out[[i]]<-createSaemixOutcome(discreteOutcome(type="binary"), names(outcome)[i])
        if(outcome[i]=="categorical") smx.out[[i]]<-createSaemixOutcome(discreteOutcome( type="categorical"), names(outcome)[i])
        if(outcome[i]=="count") smx.out[[i]]<-createSaemixOutcome(discreteOutcome( type="count"), names(outcome)[i])
        if(outcome[i]=="event") smx.out[[i]]<-createSaemixOutcome(discreteOutcome(type="event"), names(outcome)[i])
      }
      names(smx.out)<-names(outcome)
    }
    outcome<-smx.out
  } else {
    if(is(outcome, "list")) {
      smx.out<-vector(length(outcome), mode="list")
      names(smx.out)<-names(outcome)
      for(i in 1:length(outcome)) {
        if(is(outcome[[i]], "SaemixOutcome")) {
          smx.out[[i]]<-outcome[[i]]
          smx.out[[i]]@name<-names(outcome)[i]
        } else { # lists created by the creator functions
          x1<-try(createSaemixOutcome(outcome[[i]], name=names(outcome)[i]))
          if(is(x1,"try-error")) {
            if(verbose) message("Error in saemixModel:\n Outcome must either be a vector of (continuous) response names, a list of outcomes defined through continuousOutcome() and discreteOutcome(), or a list of outcomes defined through createSaemixOutcome().\n")
            return("Creation of SaemixModel failed") 
          } else smx.out[[i]]<-x1
        } 
      }
      outcome<-smx.out
    } else {
      if(verbose) message("Error in saemixModel:\n Outcome must either be a vector of (continuous) response names, a list of outcomes defined through continuousOutcome() and discreteOutcome(), or a list of outcomes defined through createSaemixOutcome().\n")
      return("Creation of SaemixModel failed") 
    }
  }
  return(outcome)
}



# Function to validate and format the list of variability levels of the model (internal)

# Eco 12/04/2022 - integrate these tests into the creator function, to keep the initialize of the classes as simple as can be
# if(!missing(mu.start) & !(missing(parameter))) {
#   if(!is.null(dim(mu.start))) {
#     mu.start<-mu.start[1,] 
#     if(missing(beta.start)) beta.start<-c(t(mu.start[-c(1),]))
#   } else {
#     mu.start<-mu.start
#     if(missing(beta.start)) beta.start<-numeric()
#   }
#   if(length(mu.start)<length(parameter)) {
#     if(verbose) message("mu.start argument should have the same length as the parameter vector, the starting value of additional parameters will be set to 1")
#     mu.start<-c(mu.start, rep(1,length(parameter)-length(mu.start)))
#   }
#   mu.start<-mu.start[1:length(parameter)] # if mu.start longer than the number of parameters
#   npar<-length(mu.start)
# }
# if(missing(parameter)) { # if mu.start is a named vector, use the names for the parameters (legacy)
#   if(is.null(dim(mu.start))) { # mu.start is a vector
#     npar<-length(mu.start)
#     if(is.null(names(mu.start))) parameter<-paste0("theta",1:npar) else parameter<-names(mu.start)
#     mu.start<-mu.start
#     if(missing(beta.start)) beta.start<-numeric()
#   } else {
#     npar<-dim(mu.start)[2] # mu.start is an array (legacy)
#     if(is.null(colnames(mu.start))) parameter<-paste0("theta",1:npar) else parameter<-colnames(mu.start)
#     if(missing(beta.start)) beta.start<-c(t(mu.start[-c(1),]))
#     mu.start<-mu.start[1,]
#   }
# } else {
#   if(is.numeric(parameter)) {
#     if(missing(mu.start)) {
#       mu.start<-unname(parameter)
#     }
#     npar<-length(parameter)
#     parameter<-names(parameter)
#     if(is.null(parameter)) parameter<-paste0("theta",1:npar)
#   } # here add a test if we decide to give parameter a more structured class
# }
# if(missing(mu.start)) {# if the name of the parameters are given but not their initial value, initialise to 1
#   npar<-length(parameter)
#   mu.start<-mu.start<-rep(1,npar)
#   beta.start<-numeric()
# }
# if(missing(transform.par)) transform.par<-rep(0, npar)
# if(length(transform.par)!=npar) transform.par<-rep(transform.par, length.out=npar)
# .Object@transform.par<-transform.par
# if(!missing(mu.fix)) {
#   if(length(mu.fix)==npar) .Object@mu.fix<-mu.fix else 
#     if(verbose) message("Size mismatch between mu.fix and the number of parameters, ignoring")
# }
# 
