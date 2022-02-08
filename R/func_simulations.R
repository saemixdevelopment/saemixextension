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
#' @param nsim Number of simulations to perform. Defaults to the nb.simpred
#' element in options
#' @param seed if non-null, seed used to initiate the random number generator (defaults to NULL)
#' @param predictions Whether the simulated parameters should be used to
#' compute predictions. Defaults to TRUE. If FALSE, only individual parameters are simulated.
#' @param res.var Whether residual variability should be added to the
#' predictions. Defaults to TRUE
#' @param uncertainty Uses uncertainty (currently not implemented). Defaults to FALSE
#' @param \dots additional arguments, unused (included for compatibility with the generic)
#' 
#' @details This function replaces the previous function (simul.saemix), which will be deprecated in future versions
#' but can still be called as previously for compatibility purposes.
#' 
#' @aliases simul.saemix
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

simulate.SaemixObject<-function(object, nsim, seed, predictions=TRUE,res.var=TRUE,uncertainty=FALSE,...) {
  # Simulate individual parameters from the population distribution
  # predictions: if TRUE, use the parameters to predict observations
  # res.var: if TRUE, add residual error to the predictions to obtain simulated data
  # uncertainty: if TRUE, add uncertainty when simulating (not implemented yet)
  saemix.model<-object["model"]
  saemix.data<-object["data"]
  saemix.res<-object["results"]
  xind<-saemix.data["data"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
  if(missing(nsim)) nsim<-object["options"]$nb.sim
  if(missing(seed)) seed<-NULL
  if(!is.null(seed)) set.seed(seed)
  
  N<-saemix.data["N"]
  ind.eta<-saemix.model["indx.omega"]
  nb.etas<-length(ind.eta)
  NM <- N*nsim  
  domega<-cutoff(mydiag(saemix.res["omega"][ind.eta, ind.eta]),.Machine$double.eps)
  omega.eta<-saemix.res["omega"][ind.eta,ind.eta]
  omega.eta<-omega.eta-mydiag(mydiag(saemix.res["omega"][ind.eta,ind.eta]))+mydiag(domega)
  chol.omega<-chol(omega.eta)
  
  phiM<-mean.phiM<-do.call(rbind,rep(list(saemix.res["mean.phi"]),nsim))
  etaM<-matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
  phiM[,ind.eta]<-mean.phiM[,ind.eta]+etaM
  psiM<-transphi(phiM,saemix.model["transform.par"])
  
  index<-rep(1:N,times=saemix.data["nind.obs"])
  if(predictions) {
    IdM<-kronecker(c(0:(nsim-1)),rep(N,saemix.data["ntot.obs"]))+rep(index,nsim)
    XM<-do.call(rbind,rep(list(xind), nsim))
    pres<-saemix.res["respar"]
    sim.pred<-sim.data<-NULL
    fpred<-saemix.model["model"](psiM, IdM, XM)
    sim.pred<-fpred
    if(res.var) {
      ind.exp<-which(saemix.model["error.model"]=="exponential")
      for(ityp in ind.exp) fpred[XM$ytype==ityp]<-log(cutoff(fpred[XM$ytype==ityp]))
      gpred<-error(fpred,pres,XM$ytype)
      eps<-rnorm(length(fpred))
      sim.data<-fpred+gpred*eps
    }
    datasim<-data.frame(idsim=rep(index,nsim),irep=rep(1:nsim, each=saemix.data["ntot.obs"]),ypred=sim.pred,ysim=sim.data)
  } else {
    sim.pred<-sim.data<-IdM<-c()
    datasim<-data.frame(idsim=rep(index,nsim),irep=rep(1:nsim, each=saemix.data["ntot.obs"]))
  }
  sim.psi<-data.frame(id=rep(unique(saemix.data["data"][, saemix.data["name.group"]]),nsim),psiM)
  colnames(sim.psi)<-c(saemix.data["name.group"],saemix.model["name.modpar"])
  ysim<-new(Class="SaemixSimData",saemix.data,datasim)
  ysim["sim.psi"]<-sim.psi
  object["sim.data"]<-ysim
  
  return(object)
}

#' @export 

simul.saemix<-function(object, nsim, seed, predictions=TRUE,res.var=TRUE,uncertainty=FALSE,...) {
  simulate.SaemixObject(object=object, nsim=nsim, seed=seed, predictions=predictions,res.var=res.var,uncertainty=uncertainty,...)
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
#' @param object an saemixObject object returned by the \code{\link{saemix}} function
#' @param nsim Number of simulations to perform. Defaults to the nb.simpred
#' element in options
#' @param simulate.function a function with the same structure as the model functions (arguments 
#' (psi, id, xidep)) which returns simulated responses when given a set of individual parameters (psi) 
#' and the corresponding predictors (xidep)
#' @param seed if non-null, seed used to initiate the random number generator (defaults to NULL)
#' @param uncertainty Uses uncertainty (currently not implemented). Defaults to FALSE
#' 
#' @details This function calls simulate.SaemixObject with the prediction=FALSE option to 
#' simulate individual parameters, then the simulate.function to obtain corresponding predictions.
#' 
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>
#' 
#' @seealso \code{\link{SaemixObject}},\code{\link{saemix}},
#' \code{\link{simulate.SaemixObject}}
#' 
#' @keywords model
#' 
#' @export 

simulateDiscreteSaemix <- function(object, simulate.function, nsim, seed, uncertainty=FALSE) {
  # Simulate individual parameters from the population distribution
  ## object: an SaemixObject object resulting from a call to saemix()
  ## simulate.function: a function matching the model object@model@model and simulating outcomes given predictors and individual parameters
  ## nsim: number of simulations
  ## uncertainty: if TRUE, add uncertainty when simulating (not implemented yet)
  if(missing(nsim)) nsim<-object["options"]$nb.sim
  if(missing(seed)) seed<-NULL
  object<-simulate.SaemixObject(object, nsim=nsim, seed=seed, predictions=FALSE, uncertainty=uncertainty)
  simpar<-object@sim.data@sim.psi
  
  # Simulate observations using these parameters and the simulate.function to simulate from the same model
  xidep<-object@data@data[,object@data@name.predictors]
  id1<-object@data@data[,"index"]
  nsuj<-object@data@N
  datasim<-object@sim.data@datasim
  datasim$ysim<-NA
  for(irep in 1:nsim) {
    psi1<-simpar[(1+(irep-1)*nsuj):(irep*nsuj),-c(1)]
    ysim<-simulate.function(psi1, id1, xidep)
#    if(sum(is.na(ysim))>0) cat(irep,"\n")
    datasim$ysim[datasim$irep==irep]<-ysim
  }
  object@sim.data@datasim<-datasim
  return(object)
}
