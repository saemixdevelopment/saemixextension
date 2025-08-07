
#' @include aaa_generics.R
NULL

####################################################################################
####			SaemixRes class - definition				####
####################################################################################

###############################
# Definition with initialise

#' Class "SaemixRes"
#' 
#' An object of the SaemixRes class, representing the results of a fit through the SAEM algorithm.
#' 
#' @name SaemixRes-class 
#' @docType class
#' @aliases SaemixRes SaemixRes-class
#' @aliases plot,SaemixRes print,SaemixRes showall,SaemixRes show,SaemixRes 
#' 
#' @section Objects from the Class: 
#' An object of the SaemixData class can be created by using the function \code{\link{saemixData}} and contain the following slots:
#' @slot modeltype string giving the type of model used for analysis
#' @slot status string indicating whether a model has been run successfully; set to "empty" at initialisation, used to pass on error messages or fit status
#' @slot name.fixed a vector containing the names of the fixed parameters in the model
#' @slot name.random a vector containing the names of the random parameters in the model
#' @slot name.sigma a vector containing the names of the parameters of the residual error model
#' @slot npar.est the number of parameters estimated (fixed, random and residual)
#' @slot nbeta.random the number of estimated fixed effects for the random parameters in the model
#' @slot nbeta.fixed the number of estimated fixed effects for the non random parameters in the model
#' @slot fixed.effects a vector giving the estimated h(mu) and betas
#' @slot fixed.psi  a vector giving the estimated h(mu)
#' @slot betas a vector giving the estimated mu
#' @slot betaC a vector with the estimates of the fixed effects for covariates
#' @slot omega the estimated variance-covariance matrix
#' @slot respar the estimated parameters of the residual error model
#' @slot fim the Fisher information matrix
#' @slot se.fixed  a vector giving the estimated standard errors of estimation for the fixed effect parameters
#' @slot se.omega a vector giving the estimated standard errors of estimation for Omega 
#' @slot se.cov a matrix giving the estimated SE for each term of the covariance matrix (diagonal elements represent the SE on the variances of the random effects and off-diagonal elements represent the SE on the covariance terms)
#' @slot se.respar a vector giving the estimated standard errors of estimation for the parameters of the residual variability
#' @slot conf.int a dataframe containing the estimated parameters, their estimation error (SE), coefficient of variation (CV), and the associated confidence intervals; the variabilities for the random effects are presented first as estimated (variances) then converted to standard deviations (SD), and the correlations are computed. For SD and correlations, the SE are estimated via the delta-method
#' @slot parpop a matrix tracking the estimates of the population parameters at each iteration
#' @slot allpar a matrix tracking the estimates of all the parameters (including covariate effects) at each iteration
#' @slot indx.fix the index of the fixed parameters (used in the estimation algorithm)
#' @slot indx.cov the index of the covariance parameters (used in the estimation algorithm)
#' @slot indx.omega the index of the random effect parameters (used in the estimation algorithm)
#' @slot indx.res the index of the residual error model parameters (used in the estimation algorithm)
#' @slot MCOV a matrix of covariates (used in the estimation algorithm)
#' @slot cond.mean.phi a matrix giving the conditional mean estimates of phi (estimated as the mean of the conditional distribution)
#' @slot cond.mean.psi a matrix giving the conditional mean estimates of psi (h(cond.mean.phi))
#' @slot cond.var.phi a matrix giving the variance on the conditional mean estimates of phi (estimated as the variance of the conditional distribution)
#' @slot cond.mean.eta  a matrix giving the conditional mean estimates of the random effect eta
#' @slot cond.shrinkage a vector giving the shrinkage on the conditional mean estimates of eta
#' @slot mean.phi a matrix giving the population estimate (Ci*mu) including covariate effects, for each subject
#' @slot map.psi a dataframe giving the MAP estimates of individual parameters
#' @slot map.phi a dataframe giving the MAP estimates of individual phi
#' @slot map.eta a matrix giving the individual estimates of the random effects corresponding to the MAP estimates
#' @slot map.shrinkage  a vector giving the shrinkage on the MAP estimates of eta
#' @slot phi individual parameters, estimated at the end of the estimation process as the average over the chains of the individual parameters sampled during the successive E-steps
#' @slot psi.samp a three-dimensional array with samples of psi from the conditional distribution
#' @slot phi.samp a three-dimensional array with samples of phi from the conditional distribution
#' @slot phi.samp.var a three-dimensional array with the variance of phi
#' @slot ll.lin log-likelihood computed by lineariation
#' @slot aic.lin Akaike Information Criterion computed by linearisation
#' @slot bic.lin Bayesian Information Criterion computed by linearisation
#' @slot bic.covariate.lin Specific Bayesian Information Criterion for covariate selection computed by linearisation
#' @slot ll.is log-likelihood computed by Importance Sampling
#' @slot aic.is Akaike Information Criterion computed by Importance Sampling
#' @slot bic.is Bayesian Information Criterion computed by Importance Sampling
#' @slot bic.covariate.is Specific Bayesian Information Criterion for covariate selection computed by Importance Sampling
#' @slot LL a vector giving the conditional log-likelihood at each iteration of the algorithm
#' @slot ll.gq log-likelihood computed by Gaussian Quadrature
#' @slot aic.gq Akaike Information Criterion computed by Gaussian Quadrature
#' @slot bic.gq Bayesian Information Criterion computed by Gaussian Quadrature
#' @slot bic.covariate.gq Specific Bayesian Information Criterion for covariate selection computed by Gaussian Quadrature
#' @slot predictions a data frame containing all the predictions and residuals in a table format
#' @slot ppred a vector giving the population predictions obtained with the population estimates
#' @slot ypred a vector giving the mean population predictions
#' @slot ipred a vector giving the individual predictions obtained with the MAP estimates
#' @slot icpred a vector giving the individual predictions obtained with the conditional estimates
#' @slot ires a vector giving the individual residuals obtained with the MAP estimates
#' @slot iwres a vector giving the individual weighted residuals obtained with the MAP estimates
#' @slot icwres a vector giving the individual weighted residuals obtained with the conditional estimates
#' @slot wres a vector giving the population weighted residuals
#' @slot npde  a vector giving the normalised prediction distribution errors
#' @slot pd a vector giving the prediction discrepancies
#' 
#' @section Methods:
#'   \describe{
#'     \item{[<-}{\code{signature(x = "SaemixRes")}: replace elements of object}
#'     \item{[}{\code{signature(x = "SaemixRes")}: access elements of object}
#'     \item{initialize}{\code{signature(.Object = "SaemixRes")}: internal function to initialise object, not to be used}
#'     \item{print}{\code{signature(x = "SaemixRes")}: prints details about the object (more extensive than show)}
#'     \item{read}{\code{signature(object = "SaemixRes")}: internal function, not to be used }
#'     \item{showall}{\code{signature(object = "SaemixRes")}: shows all the elements in the object}
#'     \item{show}{\code{signature(object = "SaemixRes")}: prints details about the object}
#'     \item{summary}{\code{signature(object = "SaemixRes")}: summary of the results. Returns a list with a number of elements extracted from the results ().}
#' 	 }
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' @author Audrey Lavenu
#' @author Marc Lavielle.
#' @seealso \code{\link{saemixData}} \code{\link{SaemixModel}} \code{\link{saemixControl}} \code{\link{saemix}}
#' @examples
#' methods(class="SaemixRes")
#' 
#' showClass("SaemixRes")
#' 
#' @keywords classes
#' @exportClass SaemixRes

setClass(
  Class="SaemixRes",
  representation=representation(
    modeltype="character", # string giving the type of the model used for analysis
    status="character", # string indicating the status of the fit; set to "empty" at initialisation, used to pass on error messages or fit status
    name.fixed="character",	# names of fixed parameters in the model
    name.random="character",	# names of random effects
    name.sigma="character",	# names of parameters of residual error model
    npar.est="numeric",		# nb of parameters estimated (fixed, random & resid)
    nbeta.random="numeric", #nb of estimated fixed effects for the random parameters in the model
    nbeta.fixed="numeric", #nb of estimated fixed effects for the non random parameters in the model
    fixed.effects="numeric",	# vector with h(mu) and betas in estimation order
    fixed.psi="numeric",	# h(mu)
    betas="matrix",		# estimated mu
    betaC="numeric",		# estimated fixed effects for covariates
    omega="matrix",		# estimated omega
    respar="numeric",		# estimated residual variability
    fim="matrix",		# Fisher information matrix
    se.fixed="numeric",		# estimated SE for fixed effects
    se.omega="numeric",	# estimated SE for Omega (diagonal terms of Omega, presented as  vector)
    se.cov="matrix",	# matrix giving the estimated SE for each term of the covariance matrix (diagonal elements represent the SE on the variances of the random effects and off-diagonal elements represent the SE on the covariance terms)
    se.respar="numeric",	# estimated SE for residual variability
    conf.int="data.frame",   # Estimates, SE, CV and confidence intervals (assuming normality of the estimators)
    parpop="matrix",	# population parameters at each iteration
    allpar="matrix",	# all parameters (including covariate effects) at each iteration
    indx.fix="numeric",	# index of mean param estimated (was indx.betaI)
    indx.cov="numeric",	# index of cov param estimated (was indx.betaC)
    indx.omega="numeric",	# index of random param estimated (was i1.omega2)
    indx.res="numeric",	# index of param of residual errors estimated (was indx.res)
    MCOV="matrix",		# ECO TODO define
# Individual parameters
    cond.mean.phi="matrix",	# Cond mean estimates of Phi (mean of conditional distribution)
    cond.mean.psi="matrix",	# Cond mean estimates of Psi (h(cond.mean.phi))
    cond.var.phi="matrix",  # Variability on cond.mean.phi (var of conditional distribution)
    cond.mean.eta="matrix",	# Cond mean estimates of eta (CHECK FORMULA)
    cond.shrinkage="numeric",	# Shrinkage for cond mean estimates of Eta
    mean.phi="matrix",  # population estimate (Ci*mu) including covariate effects
    map.psi="data.frame",	# MAP estimates of individual parameters 
    map.phi="data.frame",	# MAP estimates of phi
    map.eta="matrix",		# ETAs corresponding to the MAP estimates (computed as map.phi-COV%*%MCOV in compute.eta.map)
    map.shrinkage="numeric", 	# shrinkage on MAP estimates
    phi="matrix",       # ECO TODO define
    psi.samp="array",		# nb.chains samples in the individual conditional distributions (psi)
    phi.samp="array",		# nb.chains samples in the individual conditional distributions (phi)
    phi.samp.var="array",	# variance of samples
# Statistical criteria
    ll.lin="numeric",           # for each method (linearisation, IS, GQ)
    aic.lin="numeric",          # ll=log-likelihood
    bic.lin="numeric",          # aic= Akaike Information Criterion
    bic.covariate.lin="numeric",# bic= Bayesian Information Criterion
    ll.is="numeric",        
    aic.is="numeric",
    bic.is="numeric",
    bic.covariate.is="numeric",
    LL="numeric",       # LL for each iteration in the IS algorithm
    ll.gq="numeric",
    aic.gq="numeric",
    bic.gq="numeric",
    bic.covariate.gq="numeric",
# Model predictions and residuals
		predictions="data.frame", # data frame containing all the predictions and residuals below
    ypred="numeric",        # vector of mean population predictions
    ppred="numeric",        # vector of population predictions with the population parameters
    ipred="numeric",		# vector of individual predictions with MAP
    icpred="numeric",		# vector of individual predictions with conditional estimates
    ires="numeric",		  # vector of individual residuals with MAP (ipred-x)
    iwres="numeric",		# vector of individual weighted residuals with MAP
    icwres="numeric",		# vector of individual weighted residuals with conditional estimates
    wres="numeric",		# vector of WRES (population weighted residuals)
    npde="numeric",		# vector of npde
    pd="numeric"		# vector of prediction discrepancies
  ),
  validity=function(object){
#    cat ("--- Checking SaemixRes object ---\n")
    return(TRUE)
  }
)

###############################
# initialize

#' @rdname initialize-methods
#' @param status string indicating whether a model has been run successfully; set to "empty" at initialisation, used to pass on error messages or fit status
#' @param modeltype a character string  giving the model used for analysis
#' @param name.fixed a character string  giving the name of the fixed parameters
#' @param name.random a character string  giving the name of the random parameters
#' @param fixed.effects vector with the estimates of h(mu) and betas in estimation order
#' @param fixed.psi vector with the estimates of h(mu)
#' @param betaC vector with the estimates of betas (estimated fixed effects for covariates)
#' @param betas vector with the estimates of mu
#' @param omega estimated variance-covariance matrix
#' @param respar vector with the estimates of the parameters of the residual error
#' @param cond.mean.phi matrix of size (number of subjects) x (nb of parameters) containing the conditional mean estimates of the, defined as the mean of the conditional distribution
#' @param cond.var.phi matrix of the variances on cond.mean.phi, defined as the variance of the conditional distribution
#' @param mean.phi matrix of size (number of subjects) x (nb of parameters) giving for each subject the estimates of the population parameters including covariate effects
#' @param phi matrix of size (number of subjects) x (nb of parameters) giving for each subject
#' @param phi.samp samples from the individual conditional distributions of the phi
#' @param parpop population parameters at each iteration
#' @param allpar all parameters (including covariate effects) at each iteration
#' @param MCOV design matrix C
#' 
#' @exportMethod initialize

setMethod(
  f="initialize",
  signature="SaemixRes",
  definition= function(.Object,status="empty",modeltype,name.fixed,name.random,name.sigma,fixed.effects, fixed.psi,betaC,betas,omega,respar,cond.mean.phi,cond.var.phi,mean.phi,phi, phi.samp,parpop, allpar,MCOV){
#    cat ("--- initialising SaemixRes Object --- \n")
    .Object@status<-status
    if(missing(modeltype)) modeltype<-character(0)
    .Object@modeltype<-modeltype
    if(missing(name.fixed)) name.fixed<-character(0)
    .Object@name.fixed<-name.fixed
    if(missing(name.random)) name.random<-character(0)
    .Object@name.random<-name.random
    if(missing(name.sigma)) name.sigma<-character(0)
    .Object@name.sigma<-name.sigma
    if(missing(fixed.effects)) fixed.effects<-numeric(0)
    .Object@fixed.effects<-fixed.effects
    if(missing(fixed.psi)) fixed.psi<-numeric(0)
    .Object@fixed.psi<-fixed.psi
    if(missing(betas)) betas<-matrix(nrow=0,ncol=0)
    .Object@betas<-betas
    if(missing(betaC)) betaC<-numeric(0)
    .Object@betaC<-betaC
    if(missing(omega)) omega<-matrix(nrow=0,ncol=0)
#    if(missing(omega)) omega<-matrix(data=NA,nrow=length(),ncol=length())
    .Object@omega<-omega
    if(missing(respar)) respar<-numeric(0)
    .Object@respar<-respar
    if(missing(cond.mean.phi)) cond.mean.phi<-matrix(nrow=0,ncol=0)
    .Object@cond.mean.phi<-cond.mean.phi
    if(missing(cond.var.phi)) cond.var.phi<-matrix(nrow=0,ncol=0)
    .Object@cond.var.phi<-cond.var.phi
    if(missing(mean.phi)) mean.phi<-matrix(nrow=0,ncol=0)
    .Object@mean.phi<-mean.phi
    if(missing(phi)) phi<-matrix(nrow=0,ncol=0)
    .Object@phi<-phi
    if(missing(phi.samp)) phi.samp<-matrix(nrow=0,ncol=0)
    .Object@phi.samp<-phi.samp
    if(missing(parpop)) parpop<-matrix(nrow=0,ncol=0)
    .Object@parpop<-parpop
    if(missing(allpar)) allpar<-matrix(nrow=0,ncol=0)
    .Object@allpar<-allpar
    if(missing(MCOV)) MCOV<-matrix(nrow=0,ncol=0)
    .Object@MCOV<-MCOV
#    if(missing()) <-
#    .Object@<-
# Object validation
#    validObject(.Object)
    return (.Object )
  }
)

####################################################################################
####			SaemixRes class - accesseur				####
####################################################################################

##' Get/set methods for SaemixRes object
##' 
##' Access slots of a SaemixRes object using the object\["slot"\] format
##' 
#' @param x object
#' @param i element to be replaced
#' @param j element to replace with
#' @param drop whether to drop unused dimensions
#' @keywords methods
#' @exportMethod [
#' @exportMethod [<-
#' @exportPattern "^[[:alpha:]]+"


# Getteur
setMethod(
  f ="[",
  signature = "SaemixRes" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "modeltype"={return(x@modeltype)},
    "status"={return(x@status)},
    "name.fixed"={return(x@name.fixed)},
    "name.sigma"={return(x@name.sigma)},
    "name.random"={return(x@name.random)},
    "npar.est"={return(x@npar.est)},
    "nbeta.random"={return(x@nbeta.random)},
    "nbeta.fixed"={return(x@nbeta.fixed)},
    "fixed.effects"={return(x@fixed.effects)},
    "fixed.psi"={return(x@fixed.psi)},
    "betas"={return(x@betas)},
    "betaC"={return(x@betaC)},
    "omega"={return(x@omega)},
    "respar"={return(x@respar)},
    "fim"={return(x@fim)},
    "se.fixed"={return(x@se.fixed)},
    "se.omega"={return(x@se.omega)},
    "se.cov"={return(x@se.cov)},
    "se.respar"={return(x@se.respar)},
    "conf.int"={return(x@conf.int)},
    "parpop"={return(x@parpop)},
    "allpar"={return(x@allpar)},
    "indx.fix"={return(x@indx.fix)},
    "indx.cov"={return(x@indx.cov)},
    "indx.omega"={return(x@indx.omega)},
    "indx.res"={return(x@indx.res)},
    "MCOV"={return(x@MCOV)},
    "cond.mean.phi"={return(x@cond.mean.phi)},
    "cond.mean.psi"={return(x@cond.mean.psi)},
    "cond.var.phi"={return(x@cond.var.phi)},
    "cond.mean.eta"={return(x@cond.mean.eta)},
    "cond.shrinkage"={return(x@cond.shrinkage)},
    "mean.phi"={return(x@mean.phi)},
    "map.psi"={return(x@map.psi)},
    "map.phi"={return(x@map.phi)},
    "map.eta"={return(x@map.eta)},
    "map.shrinkage"={return(x@map.shrinkage)},
    "phi"={return(x@phi)},
    "phi.samp"={return(x@phi.samp)},
    "psi.samp"={return(x@psi.samp)},
    "phi.samp.var"={return(x@phi.samp.var)},
    "ll.lin"={return(x@ll.lin)},
    "aic.lin"={return(x@aic.lin)},
    "bic.lin"={return(x@bic.lin)},
    "bic.covariate.lin"={return(x@bic.covariate.lin)},
    "ll.is"={return(x@ll.is)},
    "aic.is"={return(x@aic.is)},
    "bic.is"={return(x@bic.is)},
    "bic.covariate.is"={return(x@bic.covariate.is)},
    "LL"={return(x@LL)},
    "ll.gq"={return(x@ll.gq)},
    "aic.gq"={return(x@aic.gq)},
    "bic.gq"={return(x@bic.gq)},
    "bic.covariate.gq"={return(x@bic.covariate.gq)},
        "predictions"={return(x@predictions)},
    "ypred"={return(x@ypred)},
    "ppred"={return(x@ppred)},
    "ipred"={return(x@ipred)},
    "icpred"={return(x@icpred)},
    "ires"={return(x@ires)},
    "iwres"={return(x@iwres)},
    "icwres"={return(x@icwres)},
    "wres"={return(x@wres)},
    "npde"={return(x@npde)},
    "pd"={return(x@pd)},
    stop("No such attribute\n")
   )
  }
)
#paste("    ",slotNames(saemix.res),"={return(x@",slotNames(saemix.res),")},",sep="")

#' replace names of SaemixRes
#'
#' @name [
#' @aliases [<-,SaemixRes-method
#' @docType methods
#' @rdname extract-methods

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixRes" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "modeltype"={x@modeltype<-value},
    "status"={x@status<-value},
    "name.fixed"={x@name.fixed<-value},
    "name.random"={x@name.random<-value},
    "name.sigma"={x@name.sigma<-value},
    "npar.est"={x@npar.est<-value},
    "nbeta.random"={x@nbeta.random<-value},
    "nbeta.fixed"={x@nbeta.fixed<-value},
    "fixed.effects"={x@fixed.effects<-value},
    "fixed.psi"={x@fixed.psi<-value},
    "betas"={x@betas<-value},
    "betaC"={x@betaC<-value},
    "omega"={x@omega<-value},
    "respar"={x@respar<-value},
    "fim"={x@fim<-value},
    "se.fixed"={x@se.fixed<-value},
    "se.omega"={x@se.omega<-value},
    "se.cov"={x@se.cov<-value},
    "se.respar"={x@se.respar<-value},
    "conf.int"={x@conf.int<-value},
    "parpop"={x@parpop<-value},
    "allpar"={x@allpar<-value},
    "indx.fix"={x@indx.fix<-value},
    "indx.cov"={x@indx.cov<-value},
    "indx.omega"={x@indx.omega<-value},
    "indx.res"={x@indx.res<-value},
    "MCOV"={x@MCOV<-value},
    "cond.mean.phi"={x@cond.mean.phi<-value},
    "cond.mean.psi"={x@cond.mean.psi<-value},
    "cond.var.phi"={x@cond.var.phi<-value},
    "cond.mean.eta"={x@cond.mean.eta<-value},
    "cond.shrinkage"={x@cond.shrinkage<-value},
    "mean.phi"={x@mean.phi<-value},
    "map.phi"={x@map.phi<-value},
    "map.psi"={x@map.psi<-value},
    "map.eta"={x@map.eta<-value},
    "map.shrinkage"={x@map.shrinkage<-value},
    "phi"={x@phi<-value},
    "phi.samp"={x@phi.samp<-value},
    "psi.samp"={x@psi.samp<-value},
    "phi.samp.var"={x@phi.samp.var<-value},
    "ll.lin"={x@ll.lin<-value},
    "aic.lin"={x@aic.lin<-value},
    "bic.lin"={x@bic.lin<-value},
    "bic.covariate.lin"={x@bic.covariate.lin<-value},
    "ll.is"={x@ll.is<-value},
    "aic.is"={x@aic.is<-value},
    "bic.is"={x@bic.is<-value},
    "bic.covariate.is"={x@bic.covariate.is<-value},
    "LL"={x@LL<-value},
    "ll.gq"={x@ll.gq<-value},
    "aic.gq"={x@aic.gq<-value},
    "bic.gq"={x@bic.gq<-value},
    "bic.covariate.gq"={x@bic.covariate.gq<-value},
    "predictions"={x@predictions<-value},
    "ypred"={x@ypred<-value},
    "ppred"={x@ppred<-value},
    "ipred"={x@ipred<-value},
    "icpred"={x@icpred<-value},
    "ires"={x@ires<-value},
    "iwres"={x@iwres<-value},
    "icwres"={x@icwres<-value},
    "wres"={x@wres<-value},
    "npde"={x@npde<-value},
    "pd"={x@pd<-value},
    stop("No such attribute\n")
   )
#   validObject(x)
   return(x)
  }
)

####################################################################################
####			SaemixRes class - method to print/show data		####
####################################################################################

#' @rdname print-methods
#' 
#' @param digits number of digits to use for pretty printing
#' @param map when map is TRUE the individual parameter estimates are shown (defaults to FALSE)
#' 
#' @exportMethod print

setMethod("print","SaemixRes",
  function(x,digits=2,map=FALSE,...) {
#    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
#    cat("Dataset",x@name.data,"\n")
    if(x@status %in% c("empty")) {
      cat("No fit performed yet.\n")
      return()
    }
    if(x@status %in% c("initial")) {
      cat("No fit performed yet, parameters set by user.\n")
      if(length(grep("structural",x@modeltype))>0) {
        tab<-cbind(c(x@name.fixed,x@name.sigma[x@indx.res]), c(x@fixed.effects,x@respar[x@indx.res]))
      }else{
        tab<-cbind(c(x@name.fixed), c(x@fixed.effects))
      }
      tab<-rbind(tab,
                 cbind(x@name.random,diag(x@omega)[x@indx.omega]))
      colnames(tab)<-c("Parameter","Value")
      print(tab)
      return()
    }
    cat("----------------------------------------------------\n")
    cat("-----------------  Fixed effects  ------------------\n")
    cat("----------------------------------------------------\n")
    if(length(x@se.fixed)==0) {
      if(length(grep("structural",x@modeltype))>0) {
      tab<-cbind(c(x@name.fixed,x@name.sigma[x@indx.res]), c(x@fixed.effects,x@respar[x@indx.res]))
        }else{
            tab<-cbind(c(x@name.fixed), c(x@fixed.effects))
        }
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      if(length(grep("structural",x@modeltype))>0) {
        tab<-cbind(c(x@name.fixed,x@name.sigma[x@indx.res]), c(x@fixed.effects,x@respar[x@indx.res]),c(x@se.fixed,x@se.respar[x@indx.res]))
        }else{
            tab<-cbind(c(x@name.fixed), c(x@fixed.effects),c(x@se.fixed))
        }
      tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
      if(length(x@indx.cov)>0) {
      wstat<-as.double(tab[,2])/as.double(tab[,3])
      pval<-rep("-",length(wstat))
      pval[x@indx.cov]<-2*(1-normcdf(abs(wstat[x@indx.cov])))
      tab<-cbind(tab,"p-value"=pval)
      }
      is.not.est<-which(as.double(tab[,3])<=.Machine$double.xmin)
      ncol<-dim(tab)[2]-2
      tab[is.not.est,3:dim(tab)[2]]<-rep("-",ncol)
    }
    if(digits>0) {
      for(i in 2:dim(tab)[2]) {
       xcol<-as.double(as.character(tab[,i]))
       idx<-which(!is.na(xcol))
       tab[idx,i]<-format(xcol[idx],digits=digits)
      }
    }
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("-----------  Variance of random effects  -----------\n")
    cat("----------------------------------------------------\n")
#  cat("   ECO TODO: check if Omega or Omega2 (SD or variances) and can we choose ?\n") => returns omega2, and we can't choose
    if(length(x@se.omega)==0) {
      tab<-cbind(x@name.random,diag(x@omega)[x@indx.omega])
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-cbind(x@name.random,diag(x@omega)[x@indx.omega],x@se.omega[x@indx.omega])
      tab<-cbind(tab,100*as.double(tab[,3])/as.double(tab[,2]))
      nampar<-unlist(strsplit(x@name.random,"omega2."))
      nampar<-nampar[nampar!=""]
      if(length(nampar)>1) { # covariances only if more than one random effect in the model
      for(i in 1:(length(nampar)-1)) {
        for(j in (i+1):length(nampar)) {
          if(x@omega[i,j]!=0) {
            covar<-c(paste("cov",nampar[i],nampar[j],sep="."),x@omega[i,j],x@se.cov[i,j],x@se.cov[i,j]/x@omega[i,j]*100)
            tab<-rbind(tab,covar)
          }
        }
      }
      }
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
    }
    if(digits>0) {
      for(i in 2:dim(tab)[2]) 
         tab[,i]<-format(as.double(as.character(tab[,i])),digits=digits)
    }
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("------  Correlation matrix of random effects  ------\n")
    cat("----------------------------------------------------\n")
    tab<-cov2cor(x@omega[x@indx.omega,x@indx.omega,drop=FALSE])
    if(digits>0) {
      for(i in 1:dim(tab)[2]) tab[,i]<-format(as.double(as.character(tab[,i])),digits=digits)
    }
    try(colnames(tab)<-rownames(tab)<-x@name.random)
    print(tab,quote=FALSE)
    if(length(x@ll.lin)>0 | length(x@ll.is)>0 | length(x@ll.gq)>0) {
    cat("----------------------------------------------------\n")
    cat("---------------  Statistical criteria  -------------\n")
    cat("----------------------------------------------------\n")
    if(length(x@ll.lin)>0) {
    cat("Likelihood computed by linearisation\n")
    cat("      -2LL=",(-2*x@ll.lin),"\n")
    cat("      AIC =",x@aic.lin,"\n")
    cat("      BIC =",x@bic.lin,"\n")
#  cat("   ECO TODO: verifier si ca renvoie LL ou -2LL\n"): ok renvoie -2LL
    }
    if(length(x@ll.is)>0) {
    cat("\nLikelihood computed by importance sampling\n")
    cat("      -2LL=",(-2*x@ll.is),"\n")
    cat("      AIC =",x@aic.is,"\n")
    cat("      BIC =",x@bic.is,"\n")
    }  
    if(length(x@ll.gq)>0) {
    cat("\nLikelihood computed by Gaussian quadrature\n")
    cat("      -2LL=",(-2*x@ll.gq),"\n")
    cat("      AIC =",x@aic.gq,"\n")
    cat("      BIC =",x@bic.gq,"\n")
    }
    cat("----------------------------------------------------\n")
    }
    if(length(x@map.psi)>0 & map) {
    cat("----------------------------------------------------\n")
    cat("---------------  Individual parameters  ------------\n")
    cat("----------------------------------------------------\n")
      if(dim(x@map.psi)[1]<30)
      print(x@map.psi) else {
      cat("Individual estimates for the first 30 subjects:\n")
      print(x@map.psi[1:30,])
      }
    }
  }
)

#' @rdname show-methods
#' 
#' @exportMethod show

setMethod("show","SaemixRes",
  function(object) {
#    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
    if(object@status %in% c("empty")) {
      cat("No fit performed yet.\n")
      return()
    }
    if(object@status %in% c("initial")) {
      cat("No fit performed yet, parameters set by user.\n")
      if(length(grep("structural",object@modeltype))>0) {
        tab<-cbind(c(object@name.fixed,object@name.sigma[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]))
      }else{
        tab<-cbind(c(object@name.fixed), c(object@fixed.effects))
      }
      tab<-rbind(tab,
                 cbind(object@name.random,diag(object@omega)[object@indx.omega]))
      colnames(tab)<-c("Parameter","Value")
      print(tab)
      return()
    }
    cat("Fixed effects\n")
    if(length(object@se.fixed)==0) {
      if(length(grep("structural",object@modeltype))>0) {
        tab<-cbind(c(object@name.fixed,object@name.sigma[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]))
        }else{
            tab<-cbind(c(object@name.fixed), c(object@fixed.effects))
        }
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      if(length(grep("structural",object@modeltype))>0) {
        tab<-cbind(c(object@name.fixed,object@name.sigma[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]),c(object@se.fixed,object@se.respar[object@indx.res]))
      tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
        }else{
            tab<-cbind(c(object@name.fixed), c(object@fixed.effects),c(object@se.fixed))
            tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
        }
      colnames(tab)<-c("Parameter","Estimate","  SE"," CV(%)")
      if(length(object@indx.cov)>0) {
      wstat<-as.double(tab[,2])/as.double(tab[,3])
      pval<-rep("-",length(wstat))
      pval[object@indx.cov]<-2*(1-normcdf(abs(wstat[object@indx.cov])))
      tab<-cbind(tab,"p-value"=pval)
      }
      is.not.est<-which(as.double(tab[,3])<=.Machine$double.xmin)
      ncol<-dim(tab)[2]-2
      tab[is.not.est,3:dim(tab)[2]]<-rep("-",ncol)
    }
      for(i in 2:dim(tab)[2]) {
       xcol<-as.double(as.character(tab[,i]))
       idx<-which(!is.na(xcol))
       tab[idx,i]<-format(xcol[idx],digits=3)
      }
    rownames(tab)<-rep("",dim(tab)[1])
    print(tab,quote=FALSE)

    cat("\nVariance of random effects\n")
#  cat("   ECO TODO: check if Omega or Omega2 (SD orxres1@ezrzr variances) and can we choose ?\n")
    if(length(object@se.omega)==0) {
      tab<-cbind(object@name.random,diag(object@omega)[object@indx.omega])
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-cbind(object@name.random,diag(object@omega)[object@indx.omega],object@se.omega[object@indx.omega])
      tab<-cbind(tab,100*as.double(tab[,3])/as.double(tab[,2]))
      nampar<-unlist(strsplit(object@name.random,"omega2."))
      nampar<-nampar[nampar!=""]
      if(length(nampar)>1) { # covariances only if more than one random effect in the model
      for(i in 1:(length(nampar)-1)) {
        for(j in (i+1):length(nampar)) {
          if(object@omega[i,j]!=0) {
            covar<-c(paste("cov",nampar[i],nampar[j],sep="."),object@omega[i,j],object@se.cov[i,j],object@se.cov[i,j]/object@omega[i,j]*100)
            tab<-rbind(tab,covar)
          }
        }
      }
      }
      colnames(tab)<-c("Parameter","Estimate","  SE"," CV(%)")
    }
      for(i in 2:dim(tab)[2]) 
         tab[,i]<-format(as.double(as.character(tab[,i])),digits=3)
    rownames(tab)<-rep("",dim(tab)[1])
    print(tab,quote=FALSE)
    mat1<-object@omega
    if(sum(abs(mat1-mydiag(mydiag(mat1))))>0) {
    cat("\nCorrelation matrix of random effects\n")
    tab<-cov2cor(object@omega[object@indx.omega,object@indx.omega,drop=FALSE])
    for(i in 1:dim(tab)[2]) 
      tab[,i]<-format(as.double(as.character(tab[,i])),digits=3)
    try(colnames(tab)<-rownames(tab)<-object@name.random)
    print(tab,quote=FALSE)
    }
    if(length(object@ll.lin)>0 | length(object@ll.is)>0 | length(object@ll.gq)>0) {
      cat("\nStatistical criteria\n")
    }
    if(length(object@ll.lin)>0) {
    cat("Likelihood computed by linearisation\n")
    cat("      -2LL=",(-2*object@ll.lin),"\n")
    cat("       AIC=",object@aic.lin,"\n")
    cat("       BIC=",object@bic.lin,"\n")
#  cat("   ECO TODO: verifier si ca renvoie LL ou -2LL\n"): ok renvoie -2LL
    }
    if(length(object@ll.is)>0) {
    cat("Likelihood computed by importance sampling\n")
    cat("      -2LL=",(-2*object@ll.is),"\n")
    cat("       AIC=",object@aic.is,"\n")
    cat("       BIC=",object@bic.is,"\n")
    }  
    if(length(object@ll.gq)>0) {
    cat("Likelihood computed by Gaussian quadrature\n")
    cat("      -2LL=",(-2*object@ll.gq),"\n")
    cat("       AIC=",object@aic.gq,"\n")
    cat("       BIC=",object@bic.gq,"\n")
    }
#    cat("----------------------------------------------------\n")
  }
)


#' @rdname showall-methods
#' @exportMethod showall

# Could be print, with only head of data
setMethod("showall","SaemixRes",
  function(object) {
    if(object@status %in% c("empty")) {
      cat("No fit performed yet.\n")
      return()
    }
    if(object@status %in% c("initial")) {
      cat("No fit performed yet, parameters set by user.\n")
      if(length(grep("structural",object@modeltype))>0) {
        tab<-cbind(c(object@name.fixed,object@name.sigma[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]))
      }else{
        tab<-cbind(c(object@name.fixed), c(object@fixed.effects))
      }
      tab<-rbind(tab,
                 cbind(object@name.random,diag(object@omega)[object@indx.omega]))
      colnames(tab)<-c("Parameter","Value")
      print(tab)
      return()
    }
    cat("\n----------------------------------------------------\n")
    cat("-----------------  Fixed effects  ------------------\n")
    cat("----------------------------------------------------\n")
    if(length(object@se.fixed)==0) {
      if(length(grep("structural",object@modeltype))>0) {
        tab<-cbind(c(object@name.fixed,object@name.sigma[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]))
        }else{
            tab<-cbind(c(object@name.fixed), c(object@fixed.effects))
        }
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      if(length(grep("structural",object@modeltype))>0) {
        tab<-cbind(c(object@name.fixed,object@name.sigma[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]),c(object@se.fixed,object@se.respar[object@indx.res]))
      tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
            nampar<-unlist(strsplit(object@name.random,"omega2."))
            nampar<-nampar[nampar!=""]
            if(length(nampar)>1) { # covariances only if more than one random effect in the model
              for(i in 1:(length(nampar)-1)) {
                for(j in (i+1):length(nampar)) {
                  if(object@omega[i,j]!=0) {
                    covar<-c(paste("cov",nampar[i],nampar[j],sep="."),object@omega[i,j],object@se.cov[i,j],object@se.cov[i,j]/object@omega[i,j]*100)
                    tab<-rbind(tab,covar)
                  }
                }
              }
            }
        } else{
            tab<-cbind(c(object@name.fixed), c(object@fixed.effects),c(object@se.fixed))
            tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
        }
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
      if(length(object@indx.cov)>0) {
      wstat<-as.double(tab[,2])/as.double(tab[,3])
      pval<-rep("-",length(wstat))
      pval[object@indx.cov]<-2*(1-normcdf(abs(wstat[object@indx.cov])))
      tab<-cbind(tab,"p-value"=pval)
      }
      is.not.est<-which(as.double(tab[,3])<=.Machine$double.xmin)
      ncol<-dim(tab)[2]-2
      tab[is.not.est,3:dim(tab)[2]]<-rep("-",ncol)
    }
    for(i in 2:dim(tab)[2]) {
       xcol<-as.double(as.character(tab[,i]))
       idx<-which(!is.na(xcol))
       tab[idx,i]<-format(xcol[idx],digits=3)
      }
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("-----------  Variance of random effects  -----------\n")
    cat("----------------------------------------------------\n")
#  cat("   ECO TODO: check if Omega or Omega2 (SD or variances) and can we choose ?\n")
    if(length(object@se.omega)==0) {
      tab<-cbind(object@name.random,diag(object@omega)[object@indx.omega])
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-cbind(c(object@name.fixed,object@name.sigma[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]),c(object@se.fixed,object@se.respar[object@indx.res]))
      tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
      nampar<-unlist(strsplit(object@name.random,"omega2."))
      nampar<-nampar[nampar!=""]
      for(i in 1:(length(nampar)-1)) {
        for(j in (i+1):length(nampar)) {
          if(object@omega[i,j]!=0) {
            covar<-c(paste("cov",nampar[i],nampar[j],sep="."),object@omega[i,j],object@se.cov[i,j],object@se.cov[i,j]/object@omega[i,j]*100)
            tab<-rbind(tab,covar)
          }
        }
      }
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
    }
      for(i in 2:dim(tab)[2]) 
         tab[,i]<-format(as.double(as.character(tab[,i])),digits=3)
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("---------------  Statistical criteria  -------------\n")
    cat("----------------------------------------------------\n")
    if(length(object@ll.lin)>0) {
    cat("Likelihood computed by linearisation\n")
    cat("      -2LL=",(-2*object@ll.lin),"\n")
    cat("      AIC =",object@aic.lin,"\n")
    cat("      BIC =",object@bic.lin,"\n")
#  cat("   ECO TODO: verifier si ca renvoie LL ou -2LL\n"): ok renvoie -2LL
    }
    if(length(object@ll.is)>0) {
    cat("\nLikelihood computed by importance sampling\n")
    cat("      -2LL=",(-2*object@ll.is),"\n")
    cat("      AIC =",object@aic.is,"\n")
    cat("      BIC =",object@bic.is,"\n")
    }  
    if(length(object@ll.gq)>0) {
    cat("\nLikelihood computed by Gaussian quadrature\n")
    cat("      -2LL=",(-2*object@ll.gq),"\n")
    cat("      AIC =",object@aic.gq,"\n")
    cat("      BIC =",object@bic.gq,"\n")
    }
    cat("----------------------------------------------------\n")
  }
)

####################################################################################
####			SaemixRes class - extracting residuals	and fitted values		####
####################################################################################
# resid.lm
#           function (object, type = c("working", "response", "deviance", 
#                                      "pearson", "partial"), ...) 
#           {
#             type <- match.arg(type)
#             r <- object$residuals
#             res <- switch(type, working = , response = r, deviance = , 
#                           pearson = if (is.null(object$weights)) r else r * sqrt(object$weights), 
#                           partial = r)
#             res <- naresid(object$na.action, res)
#             if (type == "partial") 
#               res <- res + predict(object, type = "terms")
#             res
#           }

#' Extract Model Residuals
#' 
#' residuals is a generic function which extracts model residuals from objects returned by modelling functions. The  abbreviated form resid is an alias for residuals
#' 
#' @name resid.saemix
#' @aliases residuals.SaemixRes resid.SaemixRes
#' @aliases residuals.saemix residuals resid resid.SaemixObject residuals.SaemixObject
#' 
#' @param object an SaemixRes or an SaemixObject object
#' @param type string determining which residuals are extracted. Possible values are: "ires" (individual residuals, default), "wres" (weighted population residuals), "npde" (normalised prediction distribution errors), "pd" (prediction discrepancies), "iwres" (individual weighted residuals) and "icwres" (conditional individual weighted residuals). See user guide for details.
#' @param ... further arguments to be passed to or from other methods
#' 
#' @return By default, individual residuals are extracted from the model object
#' 
#' @docType methods
#' @keywords methods
#' @export 

resid.SaemixRes<-function (object, type = c("ires", "wres", "npde", "pd", "iwres", "icwres"), ...) {
            type <- match.arg(type)
            res <- switch(type, ires=object@ires, wres=object@wres, npde=object@npde, iwres=object@iwres, icwres=object@icwres, pd=object@pd)
            if(length(res)>0) return(res) else message(paste("No residuals of type",type,"available\n"))
}


#' Extract Model Predictions
#' 
#' fitted is a generic function which extracts model predictions from objects returned by modelling functions
#' 
#' @name fitted.saemix
#' @aliases fitted fitted.SaemixRes
#' 
#' @param object an object of type SaemixRes or SaemixObject
#' @param type string determining which predictions are extracted. Possible values are: "ipred" (individual predictions obtained using the mode of the individual distribution for each subject, default), "ppred" (population predictions obtained using the population parameters f(E(theta))), "ypred" (mean of the population predictions (E(f(theta)))) and "icpred" (individual predictions obtained using the conditional mean of the individual distribution). See user guide for details.
#' @param ... further arguments to be passed to or from other methods
#' 
#' @return Model predictions
#' 
#' @docType methods
#' @keywords methods
#' @export

fitted.SaemixRes<-function (object, type = c("ipred", "ypred", "ppred", "icpred"), ...) {
  type <- match.arg(type)
  pred <- switch(type, ipred=object@ipred, ypred=object@ypred, ppred=object@ppred, icpred=object@icpred)
  if(length(pred)>0) return(pred) else message(paste("No fitted values of type",type,"available\n"))
}

####################################################################################
####            SaemixRes class - extract variance-covariance matrix            ####
####################################################################################

#' Extracts the Variance-Covariance Matrix for a Fitted Model Object
#' 
#' Returns the variance-covariance matrix of the main parameters of a fitted model object
#' 
#' @name vcov
#' @aliases vcov vcov.SaemixRes vcov.SaemixObject
#' 
#' @param object a fitted object from a call to saemix
#' @param ...	further arguments to be passed to or from other methods
#' 
#' @return A matrix of the estimated covariances between the parameter estimates in model. In saemix, this matrix is obtained as the inverse of the Fisher Information Matrix computed by linearisation
#' 
#' @export

vcov.SaemixRes<-function(object, ...) {
  if(object@status!="empty" && min(dim(object@fim))>0)
      try(solve(object@fim)) else NULL
}

#' @rdname vcov
#' @export

vcov.SaemixObject<-function(object, ...) {
  vcov(object@results)
}

# setMethod("vcov","SaemixRes",
#           function (object, ...) 
#           {
#             try(solve(object@fim))
#           }
# )
# 
# setMethod("vcov","SaemixObject",
#           function (object, ...) 
#           {
#             vcov(object@results)
#           }
# )


####################################################################################
####			SaemixRes class - method to plot			####
####################################################################################


####################################################################################
####				SaemixRes class - summary method			####
####################################################################################
#' @rdname summary-methods
#' 
#' @exportMethod summary

setMethod("summary","SaemixRes",
          function(object, print=TRUE, ...) {
            if(length(object@fixed.effects)==0) {
              cat("Object of class SaemixRes, no fit performed yet.\n")
              return()
            }
            digits<-2;nsmall<-2
            if(print) {
              cat("----------------------------------------------------\n")
              cat("-----------------  Fixed effects  ------------------\n")
              cat("----------------------------------------------------\n")
            }
            #    browser()
            if(length(object@se.fixed)==0) {
              if(length(grep("structural",object@modeltype))>0) {
                tab<-data.frame(c(object@name.fixed, object@name.sigma[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]))
              }else{
                tab<-data.frame(c(object@name.fixed), c(object@fixed.effects))
              }
              colnames(tab)<-c("Parameter","Estimate")
            } else {
              if(length(grep("structural",object@modeltype))>0) {
                tab<-data.frame(c(object@name.fixed, object@name.sigma[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]),c(object@se.fixed,object@se.respar[object@indx.res]), stringsAsFactors=FALSE)
              }else{
                tab<-data.frame(c(object@name.fixed), c(object@fixed.effects),c(object@se.fixed), stringsAsFactors=FALSE)
              }
              tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
              colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
              if(length(object@indx.cov)>0) {
                wstat<-as.double(tab[,2])/as.double(tab[,3])
                pval<-rep("-",length(wstat))
                pval[object@indx.cov]<-2*(1-normcdf(abs(wstat[object@indx.cov])))
                tab<-cbind(tab,"p-value"=pval,stringsAsFactors=FALSE)
              }
            }
            tab.fix<-tab
            for(i in 2:dim(tab)[2]) {
              xcol<-as.double(as.character(tab[,i]))
              idx<-which(!is.na(xcol) & xcol!="-")
              tab[idx,i]<-format(xcol[idx],digits=digits,nsmall=nsmall)
            }
            if(print) {
              print(tab,quote=FALSE)
              cat("----------------------------------------------------\n")
              cat("-----------  Variance of random effects  -----------\n")
              cat("----------------------------------------------------\n")
              #  cat("   ECO TODO: check if Omega or Omega2 (SD or variances) and can we choose ?\n")
            }
            if(length(object@se.omega)==0) {
              tab<-data.frame(object@name.random, diag(object@omega)[object@indx.omega])
              colnames(tab)<-c("Parameter","Estimate")
            } else {
              tab<-data.frame(object@name.random, diag(object@omega)[object@indx.omega], object@se.omega[object@indx.omega])
              tab<-cbind(tab,100*as.double(tab[,3])/as.double(tab[,2]))
              nampar<-unlist(strsplit(object@name.random,"omega2."))
              nampar<-nampar[nampar!=""]
              for(i in 1:(length(nampar)-1)) {
                for(j in (i+1):length(nampar)) {
                  if(object@omega[i,j]!=0) {
                    nam1<-paste("cov",nampar[i],nampar[j],sep=".")
                    covar<-data.frame(nam1, object@omega[i,j],object@se.cov[i,j],object@se.cov[i,j]/object@omega[i,j]*100)
                    names(covar)<-names(tab)
                    tab<-rbind(tab,covar)
                    tab[dim(tab),1]<-paste("cov",nampar[i],nampar[j],sep=".")
                  }
                }
              }
              colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
            }
            tab.random<-tab
            for(i in 2:dim(tab)[2]) 
              tab[,i]<-format(as.double(as.character(tab[,i])),digits=digits,nsmall=nsmall)
            if(print) {
              print(tab,quote=FALSE)
              cat("----------------------------------------------------\n")
              cat("------  Correlation matrix of random effects  ------\n")
              cat("----------------------------------------------------\n")
            }
            tab<-cov2cor(object@omega[object@indx.omega, object@indx.omega,drop=FALSE])
            tab.corr<-tab
            for(i in 1:dim(tab)[2])
              tab[,i]<-format(as.double(as.character(tab[,i])),digits=digits,nsmall=nsmall)
            try(colnames(tab)<-rownames(tab)<-object@name.random)
            if(print) print(tab,quote=FALSE)
            l1<-rep(NA,3)
            tab.ll<-data.frame(Method=c("Linearisation","Importance Sampling","Gaussian Quadrature"),"-2xLL"=l1,AIC=l1,BIC=l1)
            if(length(object@ll.lin)>0 | length(object@ll.is)>0 | length(object@ll.gq)>0) {
              if(print) {
                cat("----------------------------------------------------\n")
                cat("---------------  Statistical criteria  -------------\n")
                cat("----------------------------------------------------\n")
              }
              if(length(object@ll.lin)>0) {
                if(print) {
                  cat("Likelihood computed by linearisation\n")
                  cat("      -2LL=",(-2*object@ll.lin),"\n")
                  cat("      AIC =",object@aic.lin,"\n")
                  cat("      BIC =",object@bic.lin,"\n")
                }
                tab.ll[1,2:4]<-c((-2*object@ll.lin),object@aic.lin, object@bic.lin)
              }
              if(length(object@ll.is)>0) {
                if(print) {
                  cat("\nLikelihood computed by importance sampling\n")
                  cat("      -2LL=",(-2*object@ll.is),"\n")
                  cat("      AIC =",object@aic.is,"\n")
                  cat("      BIC =",object@bic.is,"\n")
                }
                tab.ll[2,2:4]<-c((-2*object@ll.is),object@aic.is, object@bic.is)
              }  
              if(length(object@ll.gq)>0) {
                if(print) {
                  cat("\nLikelihood computed by Gaussian quadrature\n")
                  cat("      -2LL=",(-2*object@ll.gq),"\n")
                  cat("      AIC =",object@aic.gq,"\n")
                  cat("      BIC =",object@bic.gq,"\n")
                }
                tab.ll[3,2:4]<-c((-2*object@ll.gq),object@aic.gq, object@bic.gq)
              }
              if(print) cat("----------------------------------------------------\n")
            }
            tab<-data.frame(Id=1:dim(object@cond.mean.phi)[1], object@cond.mean.psi)
            namparfix<-object@name.fixed
            try(colnames(tab)[-c(1)]<-namparfix(-grep("beta",object@name.fixed)),silent=TRUE)
            npar<-length(object@name.fixed)
            coef<-list(fixed=tab.fix[1:npar,2],random=list(map.psi=object@map.psi, cond.mean.psi=tab))
            sigma<-tab.fix[-c(1:npar),2]
            
            res<-list(modeltype=object@modeltype,fixed.effects=tab.fix,sigma=sigma,random.effects=tab.random, correlation.matrix=tab.corr,logLik=tab.ll,coefficients=coef)
            if(length(object@fim)>0) res$FIM<-object@fim
            if(length(object@ypred)>0 | length(object@ipred)>0  | length(object@ppred)>0 | length(object@icpred)>0) {
              res$fitted<-list(population=list(pop.param=object@ppred, pop.mean=object@ppred),individual=list(map.ipred=object@ipred, cond.ipred=object@icpred))
            }
            if(length(object@wres)>0 | length(object@iwres)>0  | length(object@icwres)>0 | length(object@pd)>0) {
              res$residuals<-list(population=list(wres=object@wres), individual=list(map.iwres=object@iwres,cond.iwres=object@icwres, pd=object@pd, npde=object@npde))
            }
            
            invisible(res)
          }
)

####################################################################################
