
setClass(
  Class="SaemixParameterModel",
  representation=representation(
    # model parameters
    nb.parameters="integer",	# nb of parameters in the model
    mu.start="numeric", # initial value for population parameters
    mu.fixed="numeric",	# 1 when parameter is fixed to its initial value, 0 if estimated [=1-fixed.estim=1-mu.estimated]
    transform.par="numeric",	# distribution for model parameters
    # covariate model
    covariate.model="matrix",	# covariate model
    covariate.model.fix="matrix",	# fixed parameters in the covariate model (defaults to all 0's)
    beta.start="numeric", # initial value for covariate parameters
    
    betaest.model="matrix",	# 1st line=ones, next lines=covariate model [MAYBE REMOVE]
    index.mu="numeric",		# index of mean param estimated (was indx.betaI then indx.fix)
    index.beta="numeric",		# index of covariate param estimated (was indx.betaC then indx.cov)
    index.res="numeric"		# index of param of residual errors estimated (was indx.res)
  ),
  validity=function(object){
    #    cat ("--- Checking SaemixModel object ---\n")
    if(object@nb.outcome<1 | length(object@outcome)==0) {
      message("[ SaemixModel : Error ] Please specify at least one outcome")
      return("No outcome given")
    }
    if(!is.function(object@model) || !identical(names(formals(object@model)),c("psi","id","xidep"))) {
      message("[ SaemixModel : Error ] Invalid type of model")
      return("Invalid model type")
    }
    return(TRUE)
  }
)




# New version of the SaemixModel class

setClass(
  Class="SaemixModel",
  representation=representation(
    # Outcome
    nb.outcome="integer", # number of outcome in the model (set to length(outcome))
    outcome="list", # list of outcomes in the model (class SaemixOutcome, either discrete SaemixDiscreteOutcome or continuous SaemixContinuousOutcome)
    # Model structure
    model="function", 		# name of model function
    description="character",	# model description
    modeltype="character",     # type of model (structural, for continuous responses, likelihood, for discrete responses, combined, when both types are present in the model) => defaults to structural if nb.outcome not given or 1 [MAYBE REMOVE]
    psi0="matrix",		# CI for parameter estimates
    transform.par="numeric",	# distribution for model parameters
    mu.estimated="numeric",	# 1 when parameter is estimated, 0 if fixed to its initial value [REMOVE]
    mu.fixed="numeric",	# 1 when parameter is fixed to its initial value, 0 if estimated [=1-fixed.estim=1-mu.estimated]
    covariate.model="matrix",	# covariate model
    covariate.model.fix="matrix",	# fixed parameters in the covariate model (defaults to all 0's)
    betaest.model="matrix",	# 1st line=ones, next lines=covariate model [MAYBE REMOVE]
    nb.parameters="integer",	# nb of parameters in the model
    index.mu="numeric",		# index of mean param estimated (was indx.betaI then indx.fix)
    index.beta="numeric",		# index of covariate param estimated (was indx.betaC then indx.cov)
    index.res="numeric",		# index of param of residual errors estimated (was indx.res)
    # Statistical model 
    iiv.model="matrix",	# covariance model
    iiv.init="matrix",	# CI for Omega
    index.iiv="numeric",	# index of random param estimated (was i1.omega2 then indx.omega)
    # Names
    name.modpar="character",	# name of parameters in the model (columns of psi0)
    name.mu="character",	# name of fixed effects
    name.random="character",	# name of random parameters
    name.sigma="character",	# name of residual parameters (maybe not necessary)
    name.predictors="character",# name of predictors 
    name.X="character",	# name of variable to plot on X-axis for graphs
    name.response="character",	# name of response
    name.cov="character",	# name of covariates
    Mcovariates="data.frame"	# matrix of individual covariates in the model
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
    if (npar!=dim(object@iiv.model)[1]) isize<-1
    if (npar!=dim(object@iiv.init)[1]) isize<-1
    #    cat("npar=",npar,length(object@transform.par),length(object@fixed.estim), dim(object@covariate.model)[2],dim(object@iiv.model)[1],dim(object@iiv.init)[1],"\n")
    if(isize==1) {
      message("[ SaemixModel : Error ] The number of parameters should be the same in the following elements: psi0 (initial conditions), transform.par, fixed.estim, covariate.model, and the matrices iiv.model and iiv.init should be square matrices of size equal to the number of parameters. Please check the input.")
      return("Size mismatch")
    }
    if(sum(object@fixed.estim*mydiag(object@iiv.model))==0) {
      message("[ SaemixModel : Error ] ")
      if(sum(mydiag(object@iiv.model))==0) message("At least one parameter with IIV must be included in the model.") else message("At least one parameter with IIV must be estimated and not fixed in the model.")
      return("Invalid IIV structure")
    }
    if(is.na(match(object@modeltype,c("structural","likelihood","combined")))) {
      message("[ SaemixModel : Error ] Invalid type of model")
      return("Invalid model type")
    }
    return(TRUE)
  }
)

setMethod(
  f="initialize",
  signature="SaemixModel",
  definition=function(.Object, outcome, model, description, parameter.init, parameter.distribution, 
                      covariate.model, mu.fixed, covariate.model.fixed, variability.model, verbose=TRUE){
    
    
    
    
    # Object validation
    validObject(.Object)
    return (.Object)
  }
)



# New call to create a model

xmod<-try(new(Class="SaemixModel", model=model, description=description, outcome=outcome, psi0=psi0, 
              name.modpar=name.modpar, transform.par=transform.par, fixed.estim=fixed.estim, covariate.model=covariate.model,
              covariance.model=covariance.model, omega.init=omega.init))



