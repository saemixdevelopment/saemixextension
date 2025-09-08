#' Create an npdeObject from an saemixObject 
#'
#' This function uses the npde library to compute normalised prediction distribution errors (npde)
#' and normalised prediction discrepancies  (npd). The simulations can also be used to plot VPCs. 
#' As of version 3.0, the plot functions for diagnostics involving VPC, npde or npd are deprecated and
#' the user will be redirected to the current proposed method (creating an NpdeObject and using the 
#' plot functions from the library npde).
#' 
#' @param saemixObject a fitted object resulting from a call to saemix()
#' @param nsim the number of simulations used to compute npde (1000 by default, we suggest increasing
#' it for large datasets)
#' 
#' @return An object of class \code{\link{NpdeObject}}
#'
#' @details Since version 3.0, the saemix package depends on the npde package, which computes the npd/npde and
#' produces graphs. See the documentation for \code{\link{npde}} for details on the computation methods
#' See the PDF documentation and the bookdown \url{https://iame-researchcenter.github.io/npde_bookdown/}
#' for details on the different plots available.
#' 
#' @author Emmanuelle Comets <emmanuelle.comets@@bichat.inserm.fr>
#' @seealso \code{\link{npde.graphs}}, \code{\link{gof.test}}
#' @seealso \code{\link{NpdeObject}} \code{\link{npde.plot.select}}  \code{\link{autonpde}}
#' @seealso \code{\link{npde.plot.scatterplot}} \code{\link{npde.plot.dist}}
#' 
#' @references E Comets, K Brendel, F Mentre (2008). Computing normalised prediction distribution errors
#' to evaluate nonlinear mixed-effect models: the npde add-on package for R. 
#' Computer Methods and Programs in Biomedicine, 90:154-66.
#' 
#' K Brendel, E Comets, C Laffont, C Laveille, F Mentre (2006). Metrics for external model evaluation 
#' with an application to the population pharmacokinetics of gliclazide. Pharmaceutical Research,
#' 23:2036--49
#' 
#' @references PDF documentation for npde 3.0: \url{https://github.com/ecomets/npde30/blob/main/userguide_npde_3.1.pdf}
#' 
#' @import npde
#' @export
#' 
#' @examples
#' data(theo.saemix)
#' 
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
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
#' saemix.options<-list(algorithm=c(1,0,0),seed=632545,save=FALSE,save.graphs=FALSE, 
#' displayProgress=FALSE)
#' # Not run
#' # Works interactively but not in the contained environment of CRAN (it looks for a datafile 
#' # instesad of finding the dataset in the environment)
#' # saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' # npde.obj<-npdeSaemix(saemix.fit)
#' # plot(npde.obj)
#' # plot(npde.obj, plot.type="vpc")
#' # plot(npde.obj, plot.type="covariates")
#' # plot(npde.obj, plot.type="cov.x.scatter")
#' # plot(npde.obj, plot.type="cov.ecdf")
#' 
#' @importFrom npde autonpde

npdeSaemix<-function(saemixObject, nsim=1000) {
  if(saemixObject@results@status!="fitted") {
    if(saemixObject@options$warnings) message("Please fit the model first\n")
    return()
  }
  if(saemixObject@model@modeltype!="structural") {
    if(saemixObject@options$warnings) message("The npde library currently applies only to continuous response models\n")
    return()
  }
  namobs<-saemixObject@data@data
  if(saemixObject@sim.data@nsim==0) {
    if(saemixObject@options$warnings) message("Simulating from the model to compute npde, nsim=",nsim,"\n")
    saemixObject<-simulate(saemixObject, nsim=nsim)
  }
  namsim<-data.frame(idsim=saemixObject@sim.data@datasim[,"idsim"],
                     xsim=rep(namobs[,saemixObject@data@name.X],saemixObject@sim.data@nsim),
                     ysim=saemixObject@sim.data@datasim[,"ysim"])
  npdeObject<-autonpde(namobs=namobs, namsim=namsim, iid=saemixObject@data@name.group, ix=saemixObject@data@name.X, 
                       iy=saemixObject@data@name.response, icens="cens", imdv="mdv", 
                       icov=ifelse(length(saemixObject@data@name.covariates)==0,0,saemixObject@data@name.covariates),
                       units=list(x=saemixObject@data@units$x, y=saemixObject@data@units$y))
  
  return(npdeObject)
}
