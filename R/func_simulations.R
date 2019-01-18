#######################	Model simulations ########################



#' Perform simulations under the model
#' 
#' This function is used to simulate from the model. It can be called with the
#' estimated parameters (the default), the initial parameters, or with a set of
#' parameters. The original design can be used in the simulations, or a
#' different dataset may be used with the same structure (covariates) as the
#' original design. This function is not yet implemented.
#' 
#' This function is used to produce Visual Predictive Check graphs, as well as
#' to compute the normalised prediction distribution errors (npde).
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' @param nsim Number of simulations to perform. Defaults to the nb.simpred
#' element in options
#' @param predictions Whether the simulated parameters should be used to
#' compute predictions. Defaults to TRUE
#' @param res.var Whether residual variability should be added to the
#' predictions. Defaults to TRUE
#' @param uncertainty Uses uncertainty (currently not implemented). Defaults to
#' FALSE
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixObject}},\code{\link{saemix}},
#' \code{\link{saemix.plot.data}}, \code{\link{saemix.plot.convergence}},
#' \code{\link{saemix.plot.llis}}, \code{\link{saemix.plot.randeff}},
#' \code{\link{saemix.plot.obsvspred}}, \code{\link{saemix.plot.fits}},
#' \code{\link{saemix.plot.parcov}}, \code{\link{saemix.plot.distpsi}},
#' \code{\link{saemix.plot.scatterresiduals}}, \code{\link{saemix.plot.vpc}}
#' @references Brendel, K, Comets, E, Laffont, C, Laveille, C, Mentre, F.
#' Metrics for external model evaluation with an application to the population
#' pharmacokinetics of gliclazide, Pharmaceutical Research 23 (2006),
#' 2036-2049.
#' 
#' Holford, N. The Visual Predictive Check: superiority to standard diagnostic
#' (Rorschach) plots (Abstract 738), in: 14th Meeting of the Population
#' Approach Group in Europe, Pamplona, Spain, 2005.
#' @keywords model
#' @export simul.saemix
simul.saemix<-function(saemixObject,nsim=saemixObject["options"]$nb.sim, predictions=TRUE,res.var=TRUE,uncertainty=FALSE) {
  # Simulate individual parameters from the population distribution
  # predictions: if TRUE, use the parameters to predict observations
  # res.var: if TRUE, add residual error to the predictions to obtain simulated data
  # uncertainty: if TRUE, add uncertainty when simulating (not implemented yet)
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  xind<-saemix.data["data"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
  
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
  
  if(predictions) {
    index<-rep(1:N,times=saemix.data["nind.obs"])
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
  } else {
    sim.pred<-sim.data<-IdM<-c()
  }
  sim.psi<-data.frame(id=rep(unique(saemix.data["data"][, saemix.data["name.group"]]),nsim),psiM)
  colnames(sim.psi)<-c(saemix.data["name.group"],saemix.model["name.modpar"])
  datasim<-data.frame(idsim=rep(index,nsim),irep=rep(1:nsim, each=saemix.data["ntot.obs"]),ypred=sim.pred,ysim=sim.data)
  ysim<-new(Class="SaemixSimData",saemix.data,datasim)
  ysim["sim.psi"]<-sim.psi
  saemixObject["sim.data"]<-ysim
  
  return(saemixObject)
}
