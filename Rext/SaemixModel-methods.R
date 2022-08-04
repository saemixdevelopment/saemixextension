########################################################################
# Extract variance-covariance matrix - should be a method applying to SaemixObject and SaemixModel

getOmega<-function(object, level=1) {
  if(is(object, "SaemixObject"))
    object<-object@model
  if(is.character(level)) level<-match(level, names(object@var.model))
  if(is.na(level)) {
    message("Level",level,"not found, returning first variability level")
    level<-1
  }
  if(!is(object, "SaemixModel") || length(object@var.model)<level)
    return()
  return(object@var.model[[level]]@omega)
}

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
#' Three elements need to be passed to saemixModel():
#' - model: the name of a R function computing the model in SAEMIX format, with a set header 
#' (see details and examples)
#' 
#' - outcome: the list of outcomes in the model
#'   * if outcome is missing, it is assumed that only one response is present in the model, with the default
#'   outcome being a continuous response with a constant error model with a starting value of 1
#'   If the model contains more than one response, these need to be specified using the new outcome argument.
#'   * a list of SaemixOutcome objects of different types
#'   (SaemixContinuousOutcome, SaemixDiscreteOutcome or SaemixEventOutcome): see examples and \code{\link{saemixOutcome}}
#'   for details
#'   * outcome can be also specified as name, name+type, type, or as 
#'   * the different specifications can be mixed (see examples)
#'   
#' - parameter: the model parameters
#'   * a list of SaemixParameter objects created using the saemixPar() or one of the specific distribution functions (see \code{\link{saemixPar})
#'   * if parameter is missing and npar is given, the parameters will be named theta1, theta2, ... 
#'   and their starting value will be set to 1 (if parameter is given, npar will be ignored and set to length(parameter))
#'   * if only a vector or list of parameter names is given, the parameters will be assumed to have a lognormal distribution
#'   with IIV and no covariance
#'   * if a vector or list of numeric values is given, the parameters will be assumed to have a lognormal distribution
#'   with IIV and no covariance and their starting value will be set to the numeric value given
#'   * if a vector or list of types is given, the parameters will take the distribution type specified and will 
#'   be named  theta1, theta2, ... 
#'   * if a named list or a vector is used in the two cases above, the names will be used as parameter names
#'   * the different specifications can be mixed (see examples)
#' 
#' Note that since version 3.0 the call to saemixModel() has changed considerably. A legacy
#' function, saemixModel.legacy(), is available to avoid breaking previous code, but we 
#' encourage users to take advantage of the new definition of models which is much more 
#' structured, consistent and flexible. (see  \code{\link{saemixModel.legacy}} for the old way of 
#' specifying models)
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
#' # Single response model
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
#' saemix.model<-saemixModel(model=model1cpt, outcome="concentration",
#'   description="One-compartment model with first-order absorption", 
#'   parameter=c(ka=lognormalPar(mu.start=1),V=lognormalPar(mu.start=20),CL=lognormalPar(mu.start=0.5)))
#' print(saemix.model)
#'   
#' # PK/PD model with 2 responses
#' saemix.model<-saemixModel(model=model1cptdirect, 
#'   outcome=c(concentration=continuousOutcome(model="combined1", unit="mg/L"),effect="continuous"),
#'   description="One-compartment model with first-order absorption, direct PD model", 
#'   parameter=c(ka=lognormalPar(),V=lognormalPar(mu.start=20),CL=lognormalPar(mu.start=0.5)), 
#'     ic50=normalPar(mu.start=2, rho.param="CL"))
#' print(saemix.model)
#' 
#' @export saemixModel


saemixModel<-function(model, parameter, outcome, description="", npar=0, verbose=FALSE) {
  # Testing whether model is given
  if(missing(model)) {
    if(verbose) message("Error in saemixModel:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
    return("Creation of SaemixModel failed")  
  }
  # Parameters
  misnpar<-(npar==0)
  mispardef<-missing(parameter)
  if(misnpar & mispardef) { # need one of npar (nb of parameters), mu.start or parameter (defining parameters)
    if(verbose) message("Error in saemixModel:\n Please give either the number of parameters in the model or the parameters with their starting values.\n")
    return("Creation of SaemixModel failed")  
  } 
  # Checking the type of outcome, converting to a list of SaemixOutcome
  if(missing(outcome)) outcome<-list(saemixOutcome(type="continuous")) else outcome <- convertArg2Outcome(outcome, verbose=verbose)
  # Parameters
  if(mispardef) { # no parameters are defined explicitely, so we assume they are all defaut continuous type
    parameter<-vector(mode="list", length=npar)
    for(i in 1:npar)
      parameter[[i]]<-saemixPar(name=paste0("theta",i))
  } else {
    # Checking the type of parameter, converting to a list of SaemixParameter
    parameter <- convertArg2Parameter(parameter, verbose=verbose)
    if(is(parameter,"character")) {
      if(verbose) message("Error in saemixModel:\n Parameters should be given either as a vector/list of SaemixParameter objects created by saemixPar(), or as a vector of names, types or of items names='type'.\n")
      return("Creation of SaemixModel failed")  
    }
  }
  if(!misnpar && npar>length(parameter)) { # add parameters if npar given and >length(parameter)
    for(ipar in length(parameter):npar) {
      nampar<-paste0("theta",ipar)
      parameter[[ipar]] <- saemixPar(name=nampar)
      names(parameter)[ipar]<-nampar
    }
  }
  distribution<-c("normal","lognormal","logit","probit")
  mu.start<-mu.fix<-transform.par<-c()
  for(ipar in 1:length(parameter)) {
    mu.start<-c(mu.start,parameter[[ipar]]@mu)
    mu.fix<-c(mu.fix,parameter[[ipar]]@mu.fix)
    transform.par<-c(transform.par, match(parameter[[ipar]]@distribution,distribution)-1)
  }
  # Creating covariance model
  var.level<-c()
  for(ipar in 1:length(parameter)) {
    var.level<-c(var.level,parameter[[ipar]]@omega.level)
  }
  var.level <- unique(var.level)
  # When more than 1 level of variability - create a checkNested function
  ## to check that variability levels are nested and order the levels
  ## to add omega levels to each variability level defined in the model
  # if(length(var.level)>1) {
  #   parameter<-checkNested(parameter, var.level)
  # }
  nampar<-names(parameter)
  npar<-length(nampar)
  nvarlevel<-length(var.level)
  if(nvarlevel>0) {
    var.model<-vector(mode='list', length=nvarlevel)
    for(i in 1:length(var.level)) {
      var.model[[i]] <- getVarianceModel(parameter, level=var.level[i])
      print(var.model[[i]])
    }
  }
  # Sanitise covariate definitions for each parameter - remove duplicated covariate definitions for each parameter
  parameter<-removeDuplicateCovDef(parameter)
  # Creating covariate model - same covariate model at each variability level
  lcov<-getCovariateModel(parameter)
  
  # Creating a SaemixModel object
  xobj<-new(Class="SaemixModel", model=model, outcome=outcome, parameter=nampar, 
            mu.start = mu.start, transform.par=transform.par, mu.fix=mu.fix,
            covariate.model=lcov$covariate.model, covariate.model.fix=lcov$covariate.model.fix, beta.start=lcov$beta.start,
            var.model=var.model, verbose=FALSE)
  return(xobj)
}

################################################################################################
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
                             name.sigma=character(), name.modpar=character(), 
                             transform.par=numeric(), fixed.estim=numeric(), covariate.model=matrix(nrow=0,ncol=0),
                             covariance.model=matrix(nrow=0,ncol=0), omega.init=matrix(nrow=0,ncol=0),
                             error.model=character(), error.init=numeric(), 
                             verbose=FALSE) {
  # Testing input
  if(missing(model)) {
    if(verbose) message("Error in saemixModel:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
    return("Creation of SaemixModel failed")  
  }
  if(missing(psi0) || length(psi0)==0) {
    if(verbose) cat("Error in saemixModel:\n   please provide initial estimates psi0 for at least the fixed effects.\n")
    return("Creation of SaemixModel failed")  
  }
  if(is(psi0,"matrix")) {
    if(length(name.modpar)==0) name.modpar<-colnames(psi0)
    parameter<-psi0[1,]
    names(parameter)<-name.modpar
  }  else {
    if(length(name.modpar)==0) name.modpar<-names(psi0)
    parameter<-psi0
    names(parameter)<-name.modpar
  }
  if(is.null(names(parameter))) {
    if(length(name.modpar)==0) name.modpar<-paste0("theta",1:length(parameter))
    names(parameter)<-name.modpar
  }
  npar<-length(parameter)
  lpar<-vector(mode="list",length = npar)
  if(length(fixed.estim)>0) mu.fix<-(1-fixed.estim) else mu.fix<-c(0)
  length(mu.fix)<-npar
  mu.fix[is.na(mu.fix)]<-0
  if(length(covariance.model)==0 || nrow(covariance.model)!=npar || ncol(covariance.model)!=npar) 
    covariance.model<-diag(npar)
  if(length(omega.init)==0 || nrow(omega.init)!=npar || ncol(omega.init)!=npar) 
    omega.init<-diag(npar)
  distribution<-c("normal","lognormal","logit","probit")
  if(length(transform.par)==0) transform.par<-rep(2,npar)
  length(transform.par)<-npar
  transform.par[is.na(transform.par)]<-2
  for(ipar in 1:npar) {
    lpar[[ipar]]<-saemixPar(name=name.modpar[ipar], mu.start=parameter[ipar], estimated=1-mu.fix[ipar], 
         distribution = distribution[match(transform.par[ipar],1:4)], omega.level = ifelse(omega.init[ipar,ipar]==1,"id",c()))
  }
  names(lpar)<-name.modpar
  # Response
  if(modeltype=='structural') {
    error.start<-error.init
    if(length(error.model)==0 || error.model=="constant") {
      error.model<-"constant"
      if(length(error.init)==0) error.start<-NULL else error.start<-error.init[1]
    }
    if(length(error.model)>0 & error.model=="proportional") {
      if(length(error.init)==0) error.start<-NULL else error.start<-error.init[2]
    }
    if(length(error.model)>0 & error.model=="combined") {
      error.model<-"combined2"
      if(length(error.init)==0) error.start<-NULL else error.start<-error.init
    }
    if(length(error.model)>0 || error.model=="exponential") {
      if(length(error.init)==0) error.start<-NULL else error.start<-error.init[1]
    }
    out1<-list(continuousOutcome(model=error.model, start=error.start))
  } else {
    out1<-list(discreteOutcome())
  }
  if(length(name.response)>0) names(out1)<-name.response

  xmod<-try(saemixModel(model, parameter=lpar, outcome=out1, description=description, verbose=verbose))
  if(class(xmod)!="SaemixModel") return()
  # Covariate model
  if(length(covariate.model)==0 | dim(psi0)[1]==1) beta.start<-numeric() else {
    psi.beta<-psi0[-c(1),,drop=FALSE]
    psi.beta<-do.call(rbind,rep(list(psi.beta), dim(covariate.model)[1]))
    psi.beta<-psi.beta[1:dim(covariate.model)[1],,drop=FALSE]
    beta.start<-t(psi.beta)[t(covariate.model==1)]
  }
  if(length(covariate.model)==0 || ncol(covariate.model)!=npar) covariate.model<-NULL
  if(length(covariate.model)>0) {
    xmod@covariate.model<-covariate.model
    xmod@covariate.model.fix<-covariate.model*0
    xmod@beta.start <- beta.start
  }
  if(length(rownames(covariate.model))>0) xmod@name.cov<-rownames(covariate.model)
  
#  if(class(xmod)=="SaemixModel") x1<-try(validObject(xmod),silent=FALSE) else x1<-xmod
  return(xmod)
}
