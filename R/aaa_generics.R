#######################################################
# Documentation for generic methods created by saemix
#######################################################

#' Methods for Function read
#' 
#' Reads longitudinal data to create a SaemixData object (internal)
#' 
#' @name read-methods
#' @aliases read
#' @docType methods
#' @keywords internal
#' @export
# Importing packages grDevices, stats and utils
#' @import stats
#' @import grDevices
#' @import utils

setGeneric(name="read",
           def=function(object, verbose=TRUE){standardGeneric("read")}
)

#' Methods for Function showall
#' 
#' This function is used to visualise the majority of the elements of an object
#' 
#' @name showall-methods
#' @docType methods
#' @param object showall methods are available for objects of type SaemixData,
#' SaemixModel and SaemixObject
#' @return None
#' @section Methods: 
#' \describe{ 
#' \item{list("signature(x = \"SaemixData\")")}{
#' Prints a extensive summary of a SaemixData object }
#' \item{list("signature(x = \"SaemixModel\")")}{ Prints a extensive summary of
#' a SaemixModel object }
#' \item{list("signature(x = \"SaemixObject\")")}{ Prints a extensive summary
#' of the results from a SAEMIX fit }
#' \item{list("signature(x = \"SaemixRes\")")}{ Not user-level } 
#' }
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}}, \code{\link{SaemixObject}}
#' @examples
#' # A SaemixData object
#' data(theo.saemix)
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' showall(saemix.data)
#' 
#' # A SaemixModel object
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
#' saemix.model<-saemixModel(model=model1cpt,
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' showall(saemix.model)
#' @keywords methods
#' @export

setGeneric(name="showall",
           def=function(object){standardGeneric("showall")}
)



#' Functions to extract the individual estimates of the parameters and random
#' effects
#' 
#' These three functions are used to access the estimates of individual
#' parameters and random effects.
#' 
#' The psi_i represent the individual parameter estimates. In the SAEM
#' algorithm, these parameters are assumed to be a transformation of a Gaussian
#' random vector phi_i, where the phi_i can be written as a function of the
#' individual random effects (eta_i), the covariate matrix (C_i) and the vector
#' of fixed effects (mu):
#' 
#' phi_i = C_i mu + eta_i
#' 
#' More details can be found in the PDF documentation.
#' 
#' @name psi-methods
#' 
#' @aliases psi-methods phi-methods eta-methods
#' @aliases phi,SaemixObject-method eta,SaemixObject-method psi,SaemixObject-method 
#' @aliases psi phi eta 
#' @aliases psi.SaemixObject psi.saemix phi.SaemixObject eta.SaemixObject phi.saemix eta.saemix
#' 
#' @param object an SaemixObject object returned by the \code{\link{saemix}} function
#' @param type a string specifying whether to use the MAP (type="mode") or the mean (type="mean") of the conditional distribution of the individual parameters. Defaults to mode
#' @return a matrix with the individual parameters (psi/phi) or the random effects (eta). 
#' These functions are used to access and output the estimates of
#' parameters and random effects. When the object passed to the function does
#' not contain these estimates, they are automatically computed. The object is
#' then returned (invisibly) with these estimates added to the results.
#' @section Methods: \describe{
#' \item{list("signature(object = \"SaemixObject\")")}{ please refer to the PDF
#' documentation for the models} }
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{SaemixObject}}, \code{\link{saemixControl}},
#' \code{\link{plot.saemix}}
#' @references Comets  E, Lavenu A, Lavielle M. Parameter estimation in nonlinear mixed effect models using saemix, an R implementation of the SAEM algorithm. Journal of Statistical Software 80, 3 (2017), 1-41.
#' 
#' Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.
#' 
#' Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm.
#' 20th meeting of the Population Approach Group in Europe, Athens, Greece
#' (2011), Abstr 2173.
#' @docType methods
#' @keywords methods
#' @examples 
#' 
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
#' saemix.options<-list(algorithm=c(1,0,0),seed=632545,save=FALSE,save.graphs=FALSE)
#' 
#' # Not run (strict time constraints for CRAN)
#' # saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' 
#' # psi(saemix.fit)
#' # phi(saemix.fit)
#' # eta(saemix.fit,type="mean")
#' @export
setGeneric(name="psi",
           def=function(object,type=c("mode","mean")) {standardGeneric("psi")}
)

#' @rdname psi-methods
#' @export
setGeneric(name="phi",
           def=function(object,type=c("mode","mean")) {standardGeneric("phi")}
)

#' @rdname psi-methods
#' @export
setGeneric(name="eta",
           def=function(object,type=c("mode","mean")) {standardGeneric("eta")}
)

#######################################################
# Importation of generic methods - R
#######################################################

#' @import methods
NULL

#' @import graphics
NULL

#' @importFrom stats predict
NULL

#' @importFrom stats logLik
NULL

#' @importFrom stats coef
NULL

# Suggestion when compiling data alone
#' @importFrom utils head read.table
NULL

#######################################################
# Documentation for generic methods - R
#######################################################

#' Methods for Function initialize
#' 
#' Constructor functions for Classes in the saemix package
#' 
#' @name initialize-methods
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(.Object = \"SaemixData\")")}{ create a SaemixData
#' object. Please use the \code{\link{saemixData}} function.}
#' 
#' \item{list("signature(.Object = \"SaemixModel\")")}{ create a SaemixModel
#' object Please use the \code{\link{saemixModel}} function.}
#' 
#' \item{list("signature(.Object = \"SaemixObject\")")}{ create a SaemixObject
#' object. This object is obtained after a successful call to
#' \code{\link{saemix}}}
#' 
#' \item{list("signature(.Object = \"SaemixRepData\")")}{ create a
#' SaemixRepData object}
#' 
#' \item{list("signature(.Object = \"SaemixRes\")")}{ create a SaemixRes
#' object}
#' 
#' \item{list("signature(.Object = \"SaemixSimData\")")}{ create a
#' SaemixSimData object} }
#' @keywords methods
NULL


#' Methods for Function print
#' 
#' Prints a summary of an object
#' 
#' 
#' @name print-methods
#' @aliases print.saemix print-methods print,ANY-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"ANY\")")}{ Default print function }
#' 
#' \item{list("signature(x = \"SaemixData\")")}{ Prints a summary of a
#' SaemixData object }
#' 
#' \item{list("signature(x = \"SaemixModel\")")}{ Prints a summary of a
#' SaemixModel object }
#' 
#' \item{list("signature(x = \"SaemixObject\")")}{ Prints a summary of the
#' results from a SAEMIX fit }
#' 
#' \item{list("signature(x = \"SaemixRes\")")}{ Not user-level } }
#' @keywords methods
NULL


#' Methods for Function show
#' 
#' Prints a short summary of an object
#' 
#' 
#' @name show-methods
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"ANY\")")}{ Default show function }
#' 
#' \item{list("signature(x = \"SaemixData\")")}{ Prints a short summary of a
#' SaemixData object }
#' 
#' \item{list("signature(x = \"SaemixModel\")")}{ Prints a short summary of a
#' SaemixModel object }
#' 
#' \item{list("signature(x = \"SaemixObject\")")}{ Prints a short summary of
#' the results from a SAEMIX fit }
#' 
#' \item{list("signature(x = \"SaemixRes\")")}{ Not user-level }
#' 
#' \item{list("signature(object = \"SaemixRepData\")")}{ Prints a short summary
#' of a SaemixRepData object }
#' 
#' \item{list("signature(object = \"SaemixSimData\")")}{ Prints a short summary
#' of a SaemixSimData object } }
#' @keywords methods
NULL


#' Methods for Function summary
#' 
#' Methods for function \code{summary}
#' 
#' @name summary-methods
#' @aliases summary-methods summary,ANY-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"ANY\")")}{ default summary function ?}
#' 
#' \item{list("signature(x = \"SaemixData\")")}{ summary of the data }
#' 
#' \item{list("signature(x = \"SaemixModel\")")}{ summary of the model }
#' 
#' \item{list("signature(x = \"SaemixObject\")")}{ summary of an SaemixObject}
#' 
#' }
#' @keywords methods
NULL


#' Methods for Function predict
#' 
#' Methods for function \code{predict}
#' 
#' 
#' @name predict-methods
#' @aliases predict-methods predict,ANY-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(object = \"ANY\")")}{Default predict functions}
#' 
#' \item{list("signature(object = \"SaemixObject\")")}{Computes predictions
#' using the results of an SAEM fit} }
#' @keywords methods
NULL


#' Methods for Function plot
#' 
#' Methods for function \code{plot}
#' 
#' @name plot-methods
#' @aliases plot-methods plot,ANY-method
#' @docType methods
#' @section Methods: 
#' \describe{
#' \item{list("signature(x = \"ANY\")")}{ default plot function ?}
#' \item{list("signature(x = \"SaemixData\")")}{ Plots the data. Defaults to a
#' spaghetti plot of response versus predictor, with lines joining the data for
#' one individual.}
#' \item{list("signature(x = \"SaemixModel\")")}{ Plots prediction of the model
#' }
#' \item{list("signature(x = \"SaemixObject\")")}{ This method gives access to
#' a number of plots that can be performed on a SaemixObject}
#' \item{list("signature(x = \"SaemixSimData\")")}{ Plots simulated datasets} }
#' @keywords methods
NULL


#######################################################
# Internal functions (undocumented, link only)
#######################################################

#' Internal saemix objects
#' 
#' Internal saemix objects.
#' 
#' These are not to be called by the user.
#' 
#' @name saemix.internal
#' @aliases .First.lib plotnpde ssq error.typ
#' @keywords internal
NULL
