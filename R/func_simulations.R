#######################	Model simulations ########################

#' Perform simulations under the model for an saemixObject object
#' 
#' This function is used to simulate data under the empirical design,
#' using the model and estimated parameters from a fit.
#' 
#' The simulated data can then be used to produce Visual Predictive Check graphs, as well as
#' to compute the normalised prediction distribution errors (npde).
#' 
#' @param object an saemixObject object returned by the \code{\link{saemix}} function
#' @param nsim Number of simulations to perform. Defaults to the nb.sim
#' element in options
#' @param seed if non-null, seed used to initiate the random number generator (defaults to NULL)
#' @param predictions Whether the simulated parameters should be used to compute predictions. 
#' Defaults to TRUE for continuous data, and to FALSE for non-Gaussian data models.
#' If FALSE, only individual parameters are simulated.
#' @param outcome the type of outcome (used to specify TTE or RTTE models)
#' @param res.var Whether residual variability should be added to the
#' predictions. Defaults to TRUE
#' @param uncertainty Uses uncertainty (currently not implemented). Defaults to FALSE
#' @param \dots additional arguments, unused (included for compatibility with the generic)
#' 
#' @details This function replaces the previous function (simul.saemix), which will be deprecated in future versions
#' but can still be called as previously for compatibility purposes.
#' 
#' @aliases simul.saemix simulateContinuousSaemix simulateIndividualParameters 
#' 
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' 
#' @seealso \code{\link{SaemixObject}},\code{\link{saemix}},
#' \code{\link{saemix.plot.data}}, \code{\link{saemix.plot.convergence}},
#' \code{\link{saemix.plot.llis}}, \code{\link{saemix.plot.randeff}},
#' \code{\link{saemix.plot.obsvspred}}, \code{\link{saemix.plot.fits}},
#' \code{\link{saemix.plot.parcov}}, \code{\link{saemix.plot.distpsi}},
#' \code{\link{saemix.plot.scatterresiduals}}, \code{\link{saemix.plot.vpc}}
#' 
#' @references Brendel, K, Comets, E, Laffont, C, Laveille, C, Mentre, F.
#' Metrics for external model evaluation with an application to the population
#' pharmacokinetics of gliclazide, Pharmaceutical Research 23 (2006),
#' 2036-2049.
#' 
#' Holford, N. The Visual Predictive Check: superiority to standard diagnostic
#' (Rorschach) plots (Abstract 738), in: 14th Meeting of the Population
#' Approach Group in Europe, Pamplona, Spain, 2005.
#' @keywords model
#' 
#' @importFrom stats simulate
#' @export 

simulate.SaemixObject<-function(object, nsim, seed, predictions, outcome="continuous", res.var=TRUE,uncertainty=FALSE,...) {
  # Simulate individual parameters from the population distribution
  # predictions: if TRUE, use the parameters to predict observations (only for continous data)
  # res.var: if TRUE, add residual error to the predictions to obtain simulated data (only for continuous data)
  # uncertainty: if TRUE, add uncertainty when simulating (not implemented yet)
  if(!is(object,"SaemixObject")) {
    message("The argument object should be a fitted saemixObject object\n")
    return(NULL)
  }
  if(missing(predictions))
    predictions <- ifelse(object@model@modeltype=="likelihood",FALSE, TRUE)
  if(missing(seed)) seed<-NULL
  if(!is.null(seed)) set.seed(seed)
  if(missing(nsim)) nsim<-object["options"]$nb.sim
  if(object@model@modeltype=="likelihood") {
    if(tolower(outcome) %in% c("tte","rtte"))
      object<-simulateTTESaemix(object=object, nsim=nsim, seed=seed, predictions=predictions, uncertainty=uncertainty, outcome="tte") else
        object<-simulateDiscreteSaemix(object=object, nsim=nsim, seed=seed, predictions=predictions, uncertainty=uncertainty)
  } else  
    object<-simulateContinuousSaemix(object=object, nsim=nsim, seed=seed, predictions=predictions, uncertainty=uncertainty, res.var=res.var, ...)
    return(object)
}

simulateIndividualParameters<-function(object, nsim, seed, uncertainty=FALSE) {
  saemix.model<-object["model"]
  saemix.res<-object["results"]
  if(missing(nsim)) nsim<-object["options"]$nb.sim
  if(missing(seed)) seed<-NULL
  if(!is.null(seed)) set.seed(seed)
  
  ind.eta<-saemix.model["indx.omega"]
  nb.etas<-length(ind.eta)
  NM <- object["data"]["N"]*nsim  
  domega<-cutoff(mydiag(saemix.res["omega"][ind.eta, ind.eta]),.Machine$double.eps)
  omega.eta<-saemix.res["omega"][ind.eta,ind.eta]
  omega.eta<-omega.eta-mydiag(mydiag(saemix.res["omega"][ind.eta,ind.eta]))+mydiag(domega)
  chol.omega<-chol(omega.eta)
  
  phiM<-mean.phiM<-do.call(rbind,rep(list(saemix.res["mean.phi"]),nsim))
  etaM<-matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
  phiM[,ind.eta]<-mean.phiM[,ind.eta]+etaM
  psiM<-transphi(phiM,saemix.model["transform.par"])
  return(psiM)
}

simulateContinuousSaemix<-function(object, nsim, seed, predictions, res.var=TRUE,uncertainty=FALSE,...) {
  if(!is(object,"SaemixObject")) {
    message("The argument object should be a fitted saemixObject object\n")
    return(NULL)
  }
  
  # Model for Gaussian outcomes
  saemix.model<-object["model"]
  saemix.data<-object["data"]
  saemix.res<-object["results"]
  xind<-saemix.data["data"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
  if(missing(nsim)) nsim<-object["options"]$nb.sim
  if(missing(seed)) seed<-NULL
  if(!is.null(seed)) set.seed(seed)
  
  psiM <- simulateIndividualParameters(object=object, nsim=nsim, seed=seed, uncertainty=uncertainty)
  
  N<-saemix.data["N"]
  index<-rep(1:N,times=saemix.data["nind.obs"])
  if(predictions) {
    IdM<-kronecker(c(0:(nsim-1)),rep(N,saemix.data["ntot.obs"]))+rep(index,nsim)
    XM<-do.call(rbind,rep(list(xind), nsim))
    pres<-saemix.res["respar"]
    sim.pred<-sim.data<-NULL
    fpred<-saemix.model["model"](psiM, IdM, XM)
    sim.pred<-fpred
    datasim<-data.frame(idsim=rep(saemix.data@data[,saemix.data@name.group],nsim),irep=rep(1:nsim, each=saemix.data["ntot.obs"]),ypred=sim.pred)
    if(res.var) {
      ind.exp<-which(saemix.model["error.model"]=="exponential")
      for(ityp in ind.exp) fpred[XM$ytype==ityp]<-log(cutoff(fpred[XM$ytype==ityp]))
      gpred<-error(fpred,pres,XM$ytype)
      eps<-rnorm(length(fpred))
      sim.data<-fpred+gpred*eps
      datasim$ysim<-sim.data
    }
  } else {
    sim.pred<-sim.data<-IdM<-c()
    datasim<-data.frame(idsim=rep(saemix.data@data[,saemix.data@name.group],nsim),irep=rep(1:nsim, each=saemix.data["ntot.obs"]))
  }
  sim.psi<-data.frame(id=rep(unique(saemix.data["data"][, saemix.data["name.group"]]),nsim),psiM)
  colnames(sim.psi)<-c(saemix.data["name.group"],saemix.model["name.modpar"])
  ysim<-new(Class="SaemixSimData",saemix.data,datasim)
  ysim["sim.psi"]<-sim.psi
  ysim["nsim"]<-nsim
  object["sim.data"]<-ysim
  
  return(object)
}

#' @export 

simul.saemix<-function(object, nsim, seed, predictions=TRUE,res.var=TRUE,uncertainty=FALSE,...) {
  object<-simulate.SaemixObject(object=object, nsim=nsim, seed=seed, predictions=predictions,res.var=res.var,uncertainty=uncertainty,...)
  return(object)
}
  
#######################	Simulations for models defined through their probability ########################

#' Perform simulations under the model for an saemixObject object defined by its log-likelihood
#' 
#' This function is used to simulate from discrete response models. 
#' 
#' To call this function, the user needs to define a simulate.function matching the model function
#' in the object. The function will then be used to simulate data under the empirical design,
#' using the model and estimated parameters from a fit.
#' 
#' @param object an saemixObject object returned by the \code{\link{saemix}} function. 
#' The model must contain a slot simulate.function, containing a function with the same structure as the model functions 
#' (arguments (psi, id, xidep)) which returns simulated responses when given a set of individual parameters (psi) 
#' and the corresponding predictors (xidep)
#' @param nsim Number of simulations to perform. Defaults to the nb.simpred
#' element in options
#' @param seed if non-null, seed used to initiate the random number generator (defaults to NULL)
#' @param predictions Whether the simulated parameters should be used to
#' compute predictions. Defaults to TRUE. If FALSE, only individual parameters are simulated.
#' @param uncertainty Uses uncertainty (currently not implemented). Defaults to FALSE
#' @param verbose if TRUE, prints messages (defaults to FALSE)
#' 
#' @details This function calls simulate.SaemixObject with the prediction=FALSE option to 
#' simulate individual parameters, then the simulate.function to obtain corresponding predictions.
#' 
#' @aliases simulateTTESaemix
#' 
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>
#' 
#' @seealso \code{\link{SaemixObject}},\code{\link{saemix}},
#' \code{\link{simulate.SaemixObject}}
#' 
#' @keywords model
#' 
#' @export 
#' 
  
simulateDiscreteSaemix <- function(object, nsim, seed, predictions=TRUE, uncertainty=FALSE, verbose=FALSE) {
  # Simulate individual parameters from the population distribution
  ## object: an SaemixObject object resulting from a call to saemix()
  ## simulate.function: a function matching the model object@model@model and simulating outcomes given predictors and individual parameters
  ## nsim: number of simulations
  ## uncertainty: if TRUE, add uncertainty when simulating (not implemented yet)
  if(!is(object,"SaemixObject")) {
    if(verbose) message("The argument object should be a fitted saemixObject object\n")
    return(NULL)
  }
  if(is.null(body(object@model@simulate.function)) || is(try(validObject(object@model)), "try-error")) {
    if(verbose) message("No simulate.function element in the model component, please add a simulate.function \nwith the same structure as the model functions (arguments (psi, id, xidep)), \nto return simulated responses when given a set of individual parameters (psi) 
and the corresponding predictors (xidep) \n")
    return()
    } else
      simulate.function <- object@model@simulate.function
  if(missing(nsim)) nsim<-object["options"]$nb.sim
  if(missing(seed)) seed<-NULL
  simpar <- simulateIndividualParameters(object=object, nsim=nsim, seed=seed, uncertainty=uncertainty)
  
#  object<-simulate.SaemixObject(object, nsim=nsim, seed=seed, predictions=FALSE, uncertainty=uncertainty)

  # Simulate observations using these parameters and the simulate.function to simulate from the same model
  xidep<-object@data@data[,object@data@name.predictors]
  id1<-object@data@data[,"index"]
#  id1<-object@data@data[,object@data@name.group]
  nsuj<-object@data@N
  datasim<-data.frame(idsim=rep(object@data@data[,object@data@name.group],nsim),irep=rep(1:nsim, each=object@data["ntot.obs"]))
  datasim$ysim<-NA
  for(irep in 1:nsim) {
    psi1<-simpar[(1+(irep-1)*nsuj):(irep*nsuj),,drop=FALSE]
    ysim<-simulate.function(psi1, id1, xidep)
#    if(sum(is.na(ysim))>0) cat(irep,"\n")
    datasim$ysim[datasim$irep==irep]<-ysim
  }
  sim.psi<-data.frame(id=rep(unique(object@data["data"][, object@data["name.group"]]),nsim),simpar)
  colnames(sim.psi)<-c(object@data["name.group"],object@model["name.modpar"])
  ysim<-new(Class="SaemixSimData",object@data,datasim)
  ysim["sim.psi"]<-sim.psi
  ysim["nsim"]<-nsim
  object["sim.data"]<-ysim
  return(object)
}

# Will require that time is a set variable...
# Currently doesn't work with RTTE
simulateTTESaemix <- function(object, nsim, seed, outcome="tte", predictions=TRUE, uncertainty=FALSE, verbose=FALSE) {
  # Simulate individual parameters from the population distribution
  ## object: an SaemixObject object resulting from a call to saemix()
  ## simulate.function: a function matching the model object@model@model and simulating outcomes given predictors and individual parameters
  ## nsim: number of simulations
  ## outcome: for this function, TTE or RTTE
  ## uncertainty: if TRUE, add uncertainty when simulating (not implemented yet)
  if(!is(object,"SaemixObject")) {
    if(verbose) message("The argument object should be a fitted saemixObject object\n")
    return(NULL)
  }
  if(is.null(body(object@model@simulate.function)) || is(try(validObject(object@model)), "try-error")) {
    if(verbose) message("No simulate.function element in the model component, please add a simulate.function \nwith the same structure as the model functions (arguments (psi, id, xidep)), \nto return simulated responses when given a set of individual parameters (psi) 
and the corresponding predictors (xidep) \n")
    return()
  } else
    simulate.function <- object@model@simulate.function
  if(missing(nsim)) nsim<-object["options"]$nb.sim
  if(missing(seed)) seed<-NULL
  simpar <- simulateIndividualParameters(object=object, nsim=nsim, seed=seed, uncertainty=uncertainty)
  
  #  object<-simulate.SaemixObject(object, nsim=nsim, seed=seed, predictions=FALSE, uncertainty=uncertainty)
  
  # Simulate observations using these parameters and the simulate.function to simulate from the same model
  xidep<-object@data@data[,object@data@name.predictors]
  id1<-object@data@data[,"index"]
  #  id1<-object@data@data[,object@data@name.group]
  nsuj<-object@data@N
  datasim<-data.frame(idsim=rep(object@data@data[,object@data@name.group],nsim),irep=rep(1:nsim, each=object@data["ntot.obs"]))
  datasim$ysim<-NA
  for(irep in 1:nsim) {
    psi1<-simpar[(1+(irep-1)*nsuj):(irep*nsuj),,drop=FALSE]
    ysim<-simulate.function(psi1, id1, xidep)
    #    if(sum(is.na(ysim))>0) cat(irep,"\n")
    datasim$ysim[datasim$irep==irep]<-ysim
  }
  sim.psi<-data.frame(id=rep(unique(object@data["data"][, object@data["name.group"]]),nsim),simpar)
  colnames(sim.psi)<-c(object@data["name.group"],object@model["name.modpar"])
  ysim<-new(Class="SaemixSimData",object@data,datasim)
  ysim["sim.psi"]<-sim.psi
  ysim["nsim"]<-nsim
  object["sim.data"]<-ysim
  return(object)
}
