#######################################################
# Documentation for generic methods created by saemix
#######################################################

#' Methods for Function read
#' 
#' Reads longitudinal data to create a SaemixData object (internal)
#' 
#' @name readSaemix-methods
#' @aliases readSaemix
#' @docType methods
#' @keywords internal
#' @export
# Importing packages grDevices, stats and utils
#' @import stats
#' @import grDevices
#' @import utils
## #' @importFrom utils head read.table modifyList ## needed if import all of utils ?

setGeneric(name="readSaemix",
           def=function(object, dat=NULL) standardGeneric("readSaemix")
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
           def=function(object) standardGeneric("showall")
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
#' @param type a string specifying whether to use the MAP (type="mode") or the mean (type="mean") 
#' of the conditional distribution of the individual parameters. Defaults to mode
#' @return a matrix with the individual parameters (psi/phi) or the random effects (eta). 
#' These functions are used to access and output the estimates of
#' parameters and random effects. When the object passed to the function does
#' not contain these estimates, they are automatically computed. The object is
#' then returned (invisibly) with these estimates added to the results.
#' @section Methods: \describe{
#' \item{list("signature(object = \"SaemixObject\")")}{ please refer to the PDF
#' documentation for the models} 
#' }
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{SaemixObject}}, \code{\link{saemixControl}},
#' \code{\link{plot.saemix}}
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
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
#' saemix.options<-list(algorithm=c(1,0,0),seed=632545,save=FALSE,save.graphs=FALSE, 
#' displayProgress=FALSE)
#' 
#' # Not run (strict time constraints for CRAN)
#' \donttest{
#' saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' psi(saemix.fit)
#' phi(saemix.fit)
#' eta(saemix.fit,type="mean")
#' }
#' 
#' @export
setGeneric(name="psi",
           def=function(object,type=c("mode","mean")) standardGeneric("psi")
)

#' @rdname psi-methods
#' @export
setGeneric(name="phi",
           def=function(object,type=c("mode","mean")) standardGeneric("phi")
)

#' @rdname psi-methods
#' @export
setGeneric(name="eta",
           def=function(object,type=c("mode","mean")) standardGeneric("eta")
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
#' @aliases .First.lib
#' @keywords internal
NULL

############################# help files for datasets included with the package

#' Pharmacokinetics of theophylline
#'
#' The \code{theo.saemix} data frame has 132 rows and 6 columns of data from
#' an experiment on the pharmacokinetics of theophylline. A column with gender was
#' added to the original data for demo purposes, and contains simulated data.
#'
#' Boeckmann, Sheiner and Beal (1994) report data from a study by Dr. Robert
#' Upton of the kinetics of the anti-asthmatic drug theophylline.  Twelve
#' subjects were given oral doses of theophylline then serum concentrations
#' were measured at 11 time points over the next 25 hours.
#'
#' These data are analyzed in Davidian and Giltinan (1995) and Pinheiro and
#' Bates (2000) using a two-compartment open pharmacokinetic model.
#'
#' These data are also available in the library \code{datasets} under the name
#' \code{Theoph} in a slightly modified format and including the data at time
#' 0.  Here, we use the file in the format provided in the \emph{NONMEM}
#' installation path (see the User Guide for that software for details).
#'
#' @docType data
#' @name theo.saemix
#' @usage theo.saemix
#' @format This data frame contains the following columns: 
#' \describe{
#' \item{Id}{an ordered factor with levels \code{1}, \dots{}, \code{12}
#' identifying the subject on whom the observation was made.  The ordering is
#' by time at which the observation was made.  } 
#' \item{Dose}{dose of theophylline administered orally to the subject (mg/kg).  } 
#' \item{Time}{time since drug administration when the sample was drawn (hr).  }
#' \item{Concentration}{theophylline concentration in the sample (mg/L).  } 
#' \item{Weight}{ weight of the subject (kg).  } 
#' \item{Sex}{ gender (simulated, 0=male, 1=female} }
#' 
#' @source AJ Boeckmann, LB Sheiner, SL Beal (1994),
#' \emph{NONMEM Users Guide: Part V}, NONMEM Project Group, University of
#' California, San Francisco.
#'
#' @references M Davidian, DM Giltinan (1995) \emph{Nonlinear Models for Repeated
#' Measurement Data}, Chapman & Hall (section 5.5, p. 145 and section 6.6, p.
#' 176)
#'
#' JC Pinheiro, DM Bates (2000) \emph{Mixed-effects Models in S and
#' S-PLUS}, Springer (Appendix A.29)
#' 
#' @examples
#' data(theo.saemix)
#'
#' #Plotting the theophylline data
#' plot(Concentration~Time,data=theo.saemix,xlab="Time after dose (hr)",
#' ylab="Theophylline concentration (mg/L)")
#' 
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#'   model1cpt<-function(psi,id,xidep) { 
#'     dose<-xidep[,1]
#'     tim<-xidep[,2]  
#'     ka<-psi[id,1]
#'     V<-psi[id,2]
#'     CL<-psi[id,3]
#'     k<-CL/V
#'     ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#'     return(ypred)
#'     }
#' # Default model, no covariate
#' saemix.model<-saemixModel(model=model1cpt,
#'        description="One-compartment model with first-order absorption",
#'        psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, 
#'        dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))
#'  # Note: remove the options save=FALSE and save.graphs=FALSE 
#'  # to save the results and graphs
#'  saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
#'  \donttest{
#'  # Not run (strict time constraints for CRAN)
#'  saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#'  }
#'  # Model with covariates
#'  saemix.model<-saemixModel(model=model1cpt,
#'       description="One-compartment model with first-order absorption",
#'       psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, 
#'       dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), 
#'       covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'       covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),
#'       omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")
#'       
#'  saemix.options<-list(seed=39546,save=FALSE,save.graphs=FALSE,displayProgress=FALSE)
#'  \donttest{
#'  # Not run (strict time constraints for CRAN)
#'  saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#'  }
#'  
#' @keywords datasets

NULL

#' Data simulated according to an Emax response model, in SAEM format
#'
#' The \code{PD1.saemix} and \code{PD2.saemix} data frames were simulated
#' according to an Emax dose-response model.
#'
#' @docType data
#' @name PD1.saemix
#' @aliases PD2.saemix
#' 
#' @usage PD1.saemix
#' @usage PD2.saemix
#' 
#' @format This data frame contains the following columns: 
#' \describe{
#' \item{subject}{an variable identifying the subject on whom the observation was made.  
#' The ordering is by the dose at which the observation was made.  } 
#' \item{dose}{ simulated dose.  } 
#' \item{response}{simulated response}
#' \item{gender}{ gender (0 for male, 1 for female) } }
#' @details These examples were used by P. Girard and F. Mentre for the symposium dedicated to Comparison of Algorithms Using Simulated Data Sets and Blind Analysis, that took place in Lyon, France, September 2004.
#' The datasets contain 100 individuals, each receiving 3 different doses:(0, 10, 90), (5, 25, 65) or (0, 20, 30). 
#' It was assumed that doses were given in a cross-over study with sufficient wash out period to avoid carry over. 
#' Responses (y_ij) were simulated with the following pharmacodynamic model:
#' y_ij = E0_i + D_ij Emax_i/(D_ij + ED50_i) +epsilon_ij 
#' The individual parameters were simulated according to
#'   log (E0_i) = log (E0) + eta_i1
#'   log (Emax_i) = log (Emax) + eta_i2
#'   log (E50_i) = log (E50) +  beta w_i + eta_i3
#'   
#' @details PD1.saemix contains the data simulated with a gender effect, beta=0.3.
#' PD2.saemix contains the data simulated without a gender effect, beta=0.
#' 
#' @references P Girard, F Mentre (2004). Comparison of Algorithms Using Simulated Data Sets and Blind Analysis workshop, Lyon, France. 
#' 
#' @examples
#' data(PD1.saemix)
#' saemix.data<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
#'       name.predictors=c("dose"),name.response=c("response"),
#'       name.covariates=c("gender"), units=list(x="mg",y="-",covariates=c("-")))
#' 
#' modelemax<-function(psi,id,xidep) {
#' # input:
#' #   psi : matrix of parameters (3 columns, E0, Emax, EC50)
#' #   id : vector of indices 
#' #   xidep : dependent variables (same nb of rows as length of id)
#' # returns:
#' #   a vector of predictions of length equal to length of id
#'   dose<-xidep[,1]
#'   e0<-psi[id,1]
#'   emax<-psi[id,2]
#'   e50<-psi[id,3]
#'   f<-e0+emax*dose/(e50+dose)
#'   return(f)
#' }
#' 
#' # Plotting the data
#' plot(saemix.data,main="Simulated data PD1")
#' \donttest{
#' # Compare models with and without covariates with LL by Importance Sampling
#' model1<-saemixModel(model=modelemax,description="Emax growth model", 
#'        psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,
#'        c("E0","Emax","EC50"))), transform.par=c(1,1,1),
#'        covariate.model=matrix(c(0,0,0), ncol=3,byrow=TRUE),fixed.estim=c(1,1,1))
#' 
#' model2<-saemixModel(model=modelemax,description="Emax growth model", 
#'        psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL, 
#'        c("E0","Emax","EC50"))), transform.par=c(1,1,1),
#'        covariate.model=matrix(c(0,0,1), ncol=3,byrow=TRUE),fixed.estim=c(1,1,1))
#' 
#' # SE not computed as not needed for the test
#' saemix.options<-list(algorithms=c(0,1,1),nb.chains=3,seed=765754, 
#'        nbiter.saemix=c(500,300),save=FALSE,save.graphs=FALSE,displayProgress=FALSE)
#' 
#' fit1<-saemix(model1,saemix.data,saemix.options)
#' fit2<-saemix(model2,saemix.data,saemix.options)
#' wstat<-(-2)*(fit1["results"]["ll.is"]-fit2["results"]["ll.is"])
#' 
#' cat("LRT test for covariate effect on EC50: p-value=",1-pchisq(wstat,1),"\n")
#' }
#' @keywords datasets
NULL
  
#' Heights of Boys in Oxford
#'
#' The \code{oxboys.saemix} data collects the height and age of students in Oxford measured over time.
#' The data frame has 234 rows and 4 columns.
#' 
#' @docType data
#' @name oxboys.saemix
#' 
#' @usage oxboys.saemix
#' 
#' @format This data frame contains the following columns: 
#' \describe{
#' \item{Subject}{an ordered factor giving a unique identifier for each boy in the experiment } 
#' \item{age}{a numeric vector giving the standardized age (dimensionless)} 
#' \item{height}{a numeric vector giving the height of the boy (cm)}
#' \item{Occasion}{ an ordered factor - the result of converting 'age' from a continuous variable to a count so these slightly unbalanced data can be analyzed as balanced } }
#'   
#' @details These data are described in Goldstein (1987) as data on the height of a selection of boys from Oxford, England 
#' versus a standardized age. The dataset can be found in the package \code{nlme}.
#'  We use an linear model for this data:
#'  y_ij = Base_i + slope_i x_ij +epsilon_ij 
#' 
#' @references JC Pinheiro, DM Bates (2000) \emph{Mixed-effects Models in S and S-PLUS}, Springer, New York (Appendix A.19)
#' 
#' @examples
#' data(oxboys.saemix)
#' saemix.data<-saemixData(name.data=oxboys.saemix,header=TRUE,
#'       name.group=c("Subject"),name.predictors=c("age"),name.response=c("height"),
#'       units=list(x="yr",y="cm"))
#' 
#' # plot the data
#' plot(saemix.data)
#' 
#' growth.linear<-function(psi,id,xidep) {
#'   x<-xidep[,1]
#'   base<-psi[id,1]
#'   slope<-psi[id,2]
#'   f<-base+slope*x
#'   return(f)
#' }
#' saemix.model<-saemixModel(model=growth.linear,description="Linear model",
#'       psi0=matrix(c(140,1),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
#'       transform.par=c(1,0),covariance.model=matrix(c(1,1,1,1),ncol=2,byrow=TRUE), 
#'       error.model="constant")
#' 
#' saemix.options<-list(algorithms=c(1,1,1),nb.chains=1,seed=201004,
#'       save=FALSE,save.graphs=FALSE,displayProgress=FALSE)
#' \donttest{
#' saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' }
#' 
#' @keywords datasets
NULL


#' Evolution of the weight of 560 cows, in SAEM format
#'
#' The \code{cow.saemix} data contains records of the weight of 560 cows on 9 or 10 occasions. 
#' 
#' @docType data
#' @name cow.saemix
#' 
#' @usage cow.saemix
#' 
#' @format This data frame contains the following columns: 
#' \describe{
#' \item{cow}{the unique identifier for each cow } 
#' \item{time}{time (days)} 
#' \item{weight}{a numeric vector giving the weight of the cow (kg)}
#' \item{birthyear}{year of birth (between 1988 and 1998)}
#' \item{twin}{existence of a twin (no=1, yes=2)}
#' \item{birthrank}{the rank of birth (beetween 3 and 7)}
#' }
#'   
#' @details An exponential model was assumed to describe the weight gain with time:
#'  y_ij = A_i (1- B_i exp( - K_i t_ij)) +epsilon_ij 
#' 
#' @references JC Pinheiro, DM Bates (2000) \emph{Mixed-effects Models in S and S-PLUS}, Springer, New York (Appendix A.19)
#' 
#' @examples
#' data(cow.saemix)
#' saemix.data<-saemixData(name.data=cow.saemix,header=TRUE,name.group=c("cow"), 
#'       name.predictors=c("time"),name.response=c("weight"), 
#'       name.covariates=c("birthyear","twin","birthrank"), 
#'       units=list(x="days",y="kg",covariates=c("yr","-","-")))
#' 
#' growthcow<-function(psi,id,xidep) {
#'   x<-xidep[,1]
#'   a<-psi[id,1]
#'   b<-psi[id,2]
#'   k<-psi[id,3]
#'   f<-a*(1-b*exp(-k*x))
#'   return(f)
#' }
#' saemix.model<-saemixModel(model=growthcow,
#'       description="Exponential growth model", 
#'       psi0=matrix(c(700,0.9,0.02,0,0,0),ncol=3,byrow=TRUE, 
#'         dimnames=list(NULL,c("A","B","k"))),transform.par=c(1,1,1),fixed.estim=c(1,1,1), 
#'       covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
#'       covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), 
#'       omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' saemix.options<-list(algorithms=c(1,1,1),nb.chains=1,nbiter.saemix=c(200,100), 
#'              seed=4526,save=FALSE,save.graphs=FALSE,displayProgress=FALSE)
#' 
#' # Plotting the data
#' plot(saemix.data,xlab="Time (day)",ylab="Weight of the cow (kg)")
#' \donttest{
#' saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' }
#' 
#' @keywords datasets
NULL


#' Wheat yield in crops treated with fertiliser, in SAEM format
#'
#' The\code{yield.saemix} contains data from winter wheat experiments.
#' 
#' @docType data
#' @name yield.saemix
#' 
#' @usage yield.saemix
#' 
#' @format This data frame contains the following columns: 
#' \describe{
#' \item{site}{the site number}
#' \item{dose}{dose of nitrogen fertiliser (kg/ha)} 
#' \item{yield}{grain yield (kg/ha)} 
#' \item{soil.nitrogen}{end-of-winter mineral soil nitrogen (NO3- plus NH4+) in the 0 to 90 cm layer was measured on each site/year (kg/ha)}
#' }
#'   
#' @details The data in the \code{yield.saemix} comes from 37 winter wheat experiments carried out between 1990 and 1996 
#' on commercial farms near Paris, France. Each experiment was from a different site. 
#' Two soil types were represented, a loam soil and a chalky soil. Common winter wheat varieties were used. 
#' Each experiment consisted of five to eight different nitrogen fertiliser rates, for a total of 224 nitrogen treatments. 
#' Nitrogen fertilizer was applied in two applications during the growing season. For each nitrogen treatment, 
#' grain yield (adjusted to 150 g.kg-1 grain moisture content) was measured. In addition, 
#' end-of-winter mineral soil nitrogen (NO3- plus NH4+) in the 0 to 90 cm layer was measured on each site-year 
#' during February when the crops were tillering. Yield and end-of-winter mineral soil nitrogen measurements 
#' were in the ranges 3.44-11.54 t.ha-1 , and 40-180 kg.ha-1 respectively.
#' 
#' @source Makowski, D., Wallach, D., and Meynard, J.-M (1999). Models of yield, grain protein, and residual mineral 
#' nitrogen responses to applied nitrogen for winter wheat. Agronomy Journal 91: 377-385. 
#' 
#' @examples
#' data(yield.saemix)
#' saemix.data<-saemixData(name.data=yield.saemix,header=TRUE,name.group=c("site"),
#'       name.predictors=c("dose"),name.response=c("yield"),
#'       name.covariates=c("soil.nitrogen"),units=list(x="kg/ha",y="t/ha",covariates=c("kg/ha")))
#' 
#' #  Model: linear + plateau
#' yield.LP<-function(psi,id,xidep) {
#'   x<-xidep[,1]
#'   ymax<-psi[id,1]
#'   xmax<-psi[id,2]
#'   slope<-psi[id,3]
#'   f<-ymax+slope*(x-xmax)
#'   #'  cat(length(f),"  ",length(ymax),"\n")
#'   f[x>xmax]<-ymax[x>xmax]
#'   return(f)
#' }
#' saemix.model<-saemixModel(model=yield.LP,description="Linear plus plateau model",   
#'         psi0=matrix(c(8,100,0.2,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,   
#'             c("Ymax","Xmax","slope"))),covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
#'         transform.par=c(0,0,0),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
#'             byrow=TRUE),error.model="constant")
#' 
#' saemix.options<-list(algorithms=c(1,1,1),nb.chains=1,seed=666, 
#'        save=FALSE,save.graphs=FALSE,displayProgress=FALSE)
#' 
#' # Plotting the data
#' plot(saemix.data,xlab="Fertiliser dose (kg/ha)", ylab="Wheat yield (t/ha)")
#' 
#' \donttest{
#' saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' 
#' # Comparing the likelihoods obtained by linearisation and importance sampling 
#' # to the likelihood obtained by Gaussian Quadrature
#' saemix.fit<-llgq.saemix(saemix.fit)
#' {
#'    cat("LL by Importance sampling, LL_IS=",saemix.fit["results"]["ll.is"],"\n")
#'    cat("LL by linearisation, LL_lin=",saemix.fit["results"]["ll.lin"],"\n")
#'    cat("LL by Gaussian Quadrature, LL_GQ=",saemix.fit["results"]["ll.gq"],"\n")
#' }
#' 
#' # Testing for an effect of covariate soil.nitrogen on Xmax
#' saemix.model2<-saemixModel(model=yield.LP,description="Linear plus plateau model", 
#'          psi0=matrix(c(8,100,0.2,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL, 
#'             c("Ymax","Xmax","slope"))),covariate.model=matrix(c(0,1,0),ncol=3,byrow=TRUE), 
#'          transform.par=c(0,0,0),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
#'              byrow=TRUE),error.model="constant")
#' 
#' saemix.fit2<-saemix(saemix.model2,saemix.data,saemix.options)
#' # BIC for the two models
#' {
#'   cat("Model without covariate, BIC=",saemix.fit["results"]["bic.is"],"\n")
#'   cat("Model with covariate, BIC=",saemix.fit2["results"]["bic.is"],"\n")
#'   pval<-1-pchisq(-2*saemix.fit["results"]["ll.is"]+2*saemix.fit2["results"]["ll.is"],1)
#'   cat("        LRT: p=",pval,"\n")
#' }
#' }
#' #' @keywords datasets
NULL

#' Toenail data
#'
#' The \code{toenail.saemix} data are from a multicenter study comparing two oral treatments for toe-nail infection, including information for 294 patients measured at 7 weeks, comprising a total of 1908 measurements. The outcome binary variable "onycholysis" indicates the degree of separation of the nail plate from the nail-bed (none or mild versus moderate or severe). Patients were evaluated at baseline (week 0) and at weeks 4, 8, 12, 24, 36, and 48 thereafter.
#' 
#' @docType data
#' @name toenail.saemix
#' 
#' @usage toenail.saemix
#' 
#' @format This data frame contains the following columns: 
#' \describe{
#' \item{id}{subject index in file}
#' \item{time}{time of measurement (in months)} 
#' \item{y}{degree of onycholysis (0 if none or mild, 1 if moderate or severe)} 
#' \item{treatment}{treatment indicator (1=Treatment A, 0= Treatment B)}
#' \item{visit}{visit number (visit numbers 1-7 correspond to scheduled visits at 0, 4, 8, 12, 24, 36, and 48 weeks)}
#' }
#' #'   
#' @details The data in the \code{toenail.saemix} was copied from the Toenail dataset provided by the prLogistic package. Different models
#' and analyses have been performed to describe this dataset in Molenberg and Verbeke (2000). 
#' Please refer to the PDF documentation (chapter 4, section 4.6) for more details on the analysis, including how to obtain diagnostic plots.
#' 
#' @source prLogistic package in R
#' 
#' @references M De Backer, C De Vroey, E Lesaffre, I Scheys, P De Keyser (1998). 
#' Twelve weeks of continuous oral therapy for toenail onychomycosis caused by dermatophytes: 
#' A double-blind comparative trial of terbinafine 250 mg/day versus itraconazole 200 mg/day. 
#' Journal of the American Academy of Dermatology, 38:57-63.
#' 
#' E Lesaffre, B Spiessens (2001). On the effect of the number of quadrature points in a logistic random-effects model: An example. 
#' Journal of the Royal Statistical Society, Series C, 50:325-335.
#'
#' G Verbeke, G Molenberghs (2000). Linear mixed models for longitudinal data, Springer, New York.
#' 
#' S Rabe-Hesketh, A Skrondal (2008). Multilevel and Longitudinal Modeling Using Stata. Mahwah, NJ: Lawrence Erlbaum Associates. Second Edition.
#' 
#' @examples
#' 
#' data(toenail.saemix)
#' saemix.data<-saemixData(name.data=toenail.saemix,name.group=c("id"), name.predictors=c("time","y"), 
#'  name.response="y", name.covariates=c("treatment"),name.X=c("time"))
#' binary.model<-function(psi,id,xidep) {
#'   tim<-xidep[,1]
#'   y<-xidep[,2]
#'   inter<-psi[id,1]
#'   slope<-psi[id,2]
#'   logit<-inter+slope*tim
#'   pevent<-exp(logit)/(1+exp(logit))
#'   pobs = (y==0)*(1-pevent)+(y==1)*pevent
#'   logpdf <- log(pobs)
#'   return(logpdf)
#' }
#' 
#' saemix.model<-saemixModel(model=binary.model,description="Binary model",
#'      modeltype="likelihood",
#'      psi0=matrix(c(-5,-.1,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
#'      transform.par=c(0,0), covariate.model=c(0,1),
#'      covariance.model=matrix(c(1,0,0,1),ncol=2))
#' \donttest{
#' saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, 
#'    nb.chains=10, fim=FALSE)
#' binary.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' plot(binary.fit, plot.type="convergence")
#' }
#'  
#' @keywords datasets
NULL


#' Knee pain data
#' 
#' The \code{knee.saemix} data represents pain scores recorded in a clinical study 
#' in 127 patients with sport related injuries treated with two different therapies. 
#' After 3,7 and 10 days of treatment the pain occuring during knee movement was observed.
#' 
#' @docType data
#' @name knee.saemix
#' 
#' @usage knee.saemix
#' 
#' @format This data frame contains the following columns: 
#' \describe{
#' \item{id}{subject index in file}
#' \item{time}{time of measurement (in days)} 
#' \item{y}{knee pain (0=none to 4=severe)} 
#' \item{Age}{patient age (scaled and centered)} 
#' \item{Sex}{patient gender (0=male, 1=female)}
#' \item{RD}{moderate knee pain (defined as pain score 2 or more)} 
#' \item{treatment}{treatment indicator (0=placebo, 1=treatment)}
#' \item{Age2}{patient age, squared (Age^2)} 
#' }
#' #'   
#' @details The data in the \code{knee.saemix} was reformatted from the knee dataset provided by the
#' catdata package (see data(knee, package="catdata")). A time column was added representing the day 
#' of the measurement (with 0 being the baseline value) and each observation corresponds to a different
#' line in the dataset. Treatment was recoded as 0/1 (placebo/treatment), gender as 0/1 (male/female)
#' and Age2 represents the squared of centered Age.
#' 
#' Please refer to the PDF documentation (chapter 4, section 4.6) for more details on the analysis, 
#' including examples of diagnostic plots.
#' 
#' @source catdata package in R
#' 
#' @references G Tutz (2012), Regression for Categorical Data, Cambridge University Press.
#' 
#' #' @examples
#' data(knee.saemix)
#' 
#' #' @keywords datasets
NULL


#' NCCTG Lung Cancer Data, in SAEM format
#'
#' The \code{lung.saemix} contains survival data in patients with advanced lung cancer from the North Central Cancer Treatment Group. 
#' Performance scores rate how well the patient can perform usual daily activities. This data is available in the survival library for R
#' and has been reformatted here for use in saemix (see details).
#' 
#' @docType data
#' @name lung.saemix
#' 
#' @usage lung.saemix
#' 
#' @format This data frame contains the following columns: 
#' \describe{
#' \item{id}{subject index in file}
#' \item{inst}{institution code} 
#' \item{time}{observation time since the beginning of follow-up} 
#' \item{status}{0=alive, 1=dead} 
#' \item{cens}{0=observed, 1=censored}
#' \item{sex}{patient gender (0=male, 1=female)} 
#' \item{age}{age in years} 
#' \item{ph.ecog}{ECOG performance score as rated by the physician. 0=asymptomatic, 
#' 1= symptomatic but completely ambulatory, 2= in bed <50% of the day, 3= in bed > 50% of the day but not bedbound, 
#' 4 = bedbound} 
#' \item{ph.karno}{Karnofsky performance score (bad=0-good=100) rated by physician (%)} 
#' \item{pat.karno}{Karnofsky performance score (bad=0-good=100) rated by patient (%)} 
#' \item{meal.cal}{calories consumed at meals (cal)} 
#' \item{wt.loss}{weight loss in last six months (pounds)} 
#' }
#'   
#' @details The data in the \code{lung.saemix} was reformatted from the lung cancer dataset (see data(cancer, package="survival")). 
#' Patients with missing age, sex, institution or physician assessments were removed from the dataset. Status was recoded as 1 for death 
#' and 0 for a censored event, and a censoring column was added to denote whether the patient was dead or alive at the time of 
#' the last observation. For saemix, a line at time=0 was added for all subjects. Finally, subjects were numbered consecutively from 0 to 1.
#' 
#' @source Terry Therneau from the survival package in R
#' 
#' @references CL Loprinzi, JA Laurie, HS Wieand, JE Krook, PJ Novotny, JW Kugler, et al. (1994).
#' Prospective evaluation of prognostic variables from patient-completed questionnaires. 
#' North Central Cancer Treatment Group. Journal of Clinical Oncology. 12(3):601-7.
#' 
#' @examples
#' data(lung.saemix)
#' 
#' saemix.data<-saemixData(name.data=lung.saemix,header=TRUE,name.group=c("id"),
#' name.predictors=c("time","status","cens"),name.response=c("status"),
#' name.covariates=c("age", "sex", "ph.ecog", "ph.karno", "pat.karno", "wt.loss","meal.cal"),
#' units=list(x="days",y="",covariates=c("yr","","-","%","%","cal","pounds")))
#' weibulltte.model<-function(psi,id,xidep) {
#'   T<-xidep[,1]
#'   y<-xidep[,2] # events (1=event, 0=no event)
#'   cens<-which(xidep[,3]==1) # censoring times (subject specific)
#'   init <- which(T==0)
#'   lambda <- psi[id,1] # Parameters of the Weibull model
#'   beta <- psi[id,2]
#'   Nj <- length(T)
#'   ind <- setdiff(1:Nj, append(init,cens)) # indices of events
#'   hazard <- (beta/lambda)*(T/lambda)^(beta-1) # H'
#'   H <- (T/lambda)^beta # H
#'   logpdf <- rep(0,Nj) # ln(l(T=0))=0
#'   logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))
#'   logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))
#'   return(logpdf)
#' }
#' saemix.model<-saemixModel(model=weibulltte.model,description="time model",modeltype="likelihood",
#'                 psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("lambda","beta"))),
#'                 transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE))
#' saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
#' \donttest{
#' tte.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' }
#' # The fit from saemix using the above Weibull model may be compared 
#' # to the non-parametric KM estimate
#' \dontrun{
#' library(survival)
#'   lung.surv<-lung.saemix[lung.saemix$time>0,]
#'   lung.surv$status<-lung.surv$status+1
#'   Surv(lung.surv$time, lung.surv$status) # 1=censored, 2=dead
#'   f1 <- survfit(Surv(time, status) ~ 1, data = lung.surv)
#'   xtim<-seq(0,max(lung.saemix$time), length.out=200)
#'   estpar<-tte.fit@results@fixed.effects
#'   ypred<-exp(-(xtim/estpar[1])^(estpar[2]))
#'   plot(f1, xlab = "Days", ylab = "Overall survival probability")
#'   lines(xtim,ypred, col="red",lwd=2)
#' }
#' @keywords datasets
NULL

  
#' Epilepsy count data
#' 
#' The epilepsy data from Thall and Vail (1990), available from the MASS package, records two-week seizure counts for 59 epileptics. 
#' The number of seizures was recorded for a baseline period of 8 weeks, and then patients were randomly assigned to a treatment group 
#' or a control group. Counts were then recorded for four successive two-week periods. The subject's age is the only covariate. See the
#' documentation for epil in the MASS package for details on the dataset.
#' 
#' @docType data
#' @name epilepsy.saemix 
#' 
#' @source MASS package in R
#' 
#' @references P Thall, S Vail (1990). Some covariance models for longitudinal count data with overdispersion. Biometrics 46(3):657-71.
#' 
#' @examples
#' # You need to have MASS installed to successfully run this example
#' if (requireNamespace("MASS")) {
#'   
#'   epilepsy<-MASS::epil
#'   saemix.data<-saemixData(name.data=epilepsy, name.group=c("subject"),
#'      name.predictors=c("period","y"),name.response=c("y"),
#'      name.covariates=c("trt","base", "age"), units=list(x="2-week",y="",covariates=c("","","yr")))
#'   ## Poisson model with one parameter
#'   countPoi<-function(psi,id,xidep) { 
#'     y<-xidep[,2]
#'     lambda<-psi[id,1]
#'     logp <- -lambda + y*log(lambda) - log(factorial(y))
#'     return(logp)
#'     }
#'  saemix.model<-saemixModel(model=countPoi,description="Count model Poisson",modeltype="likelihood", 
#'    psi0=matrix(c(0.5),ncol=1,byrow=TRUE,dimnames=list(NULL, c("lambda"))), transform.par=c(1))
#'  \donttest{
#'   saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
#'   poisson.fit<-saemix(saemix.model,saemix.data,saemix.options)
#'   }
#' }
#'   
#' @keywords datasets
NULL


#' Rutgers Alcohol Problem Index 
#' 
#' The RAPI data studies gender differences across two years in alcohol-related problems, as measured by the Rutgers
#' Alcohol Problem Index (RAPI; White & Labouvie, 1989). The dataset includes 3,616 repeated measures of counts 
#' representing the number of alcohol problems reported over six months period, across five time points from 818 individuals. 
#' 
#' @format This data frame contains the following columns: 
#' \describe{
#' \item{id}{subject identification number} 
#' \item{time}{time since the beginning of the study (months)}
#' \item{rapi}{the number of reported alcohol problems during a six-months period} 
#' \item{gender}{ gender (0=Men, 1=Women} }
#' 
#' @docType data
#' @name rapi.saemix
#' 
#' @source David Atkins, University of Washington
#' 
#' @references D Atkins, S Baldwin, C Zheng, R Gallop, C Neighbors C (2013). 
#' A tutorial on count regression and zero-altered count models
#' for longitudinal substance use data. Psychology of Addictive Behaviors, 27(1):166–177. 
#' 
#' C Neighbors, Lewis MA, Atkins D, Jensen MM, Walter T, Fossos N, Lee C, Larimer M (2010). 
#' Efficacy of web-based personalized normative feedback: A two-year randomized controlled trial.
#' Journal of Consulting and Clinical Psychology 78(6):898-911.
#' 
#' C Neighbors, N Barnett, D Rohsenow, S Colby, P Monti (2010). Cost-Effectiveness of a Motivational lntervention for 
#' Alcohol-Involved Youth in a Hospital Emergency Department. Journal of Studies on Alcohol and Drugs 71(3):384-394.
#' 
#' @examples
#'  data(rapi.saemix)
#'  saemix.data<-saemixData(name.data=rapi.saemix, name.group=c("id"),
#'                      name.predictors=c("time","rapi"),name.response=c("rapi"),
#'                      name.covariates=c("gender"),units=list(x="months",y="",covariates=c("")))
#' hist(rapi.saemix$rapi, main="", xlab="RAPI score", breaks=30)
#' 
#' \donttest{
#' # Fitting a Poisson model
#' count.poisson<-function(psi,id,xidep) { 
#'   time<-xidep[,1]
#'   y<-xidep[,2]
#'   intercept<-psi[id,1]
#'   slope<-psi[id,2]
#'   lambda<- exp(intercept + slope*time)
#'   logp <- -lambda + y*log(lambda) - log(factorial(y))
#'   return(logp)
#' }
#' # Gender effect on intercept and slope
#' rapimod.poisson<-saemixModel(model=count.poisson,
#'    description="Count model Poisson",modeltype="likelihood",   
#'    psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
#'    transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)),
#'     covariance.model =matrix(data=1, ncol=2, nrow=2),
#'     covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE))
#' saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, fim=FALSE)
#' poisson.fit<-saemix(rapimod.poisson,saemix.data,saemix.options)
#' 
#' # Fitting a ZIP model
#' count.poissonzip<-function(psi,id,xidep) {
#'   time<-xidep[,1]
#'   y<-xidep[,2]
#'   intercept<-psi[id,1]
#'   slope<-psi[id,2]
#'   p0<-psi[id,3] # Probability of zero's
#'   lambda<- exp(intercept + slope*time)
#'   logp <- log(1-p0) -lambda + y*log(lambda) - log(factorial(y)) # Poisson
#'   logp0 <- log(p0+(1-p0)*exp(-lambda)) # Zeroes
#'   logp[y==0]<-logp0[y==0]
#'   return(logp)
#' }
#' rapimod.zip<-saemixModel(model=count.poissonzip,
#'    description="count model ZIP",modeltype="likelihood",   
#'    psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,
#'    dimnames=list(NULL, c("intercept", "slope","p0"))), 
#'    transform.par=c(0,0,3), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.3,0)),
#'    covariate.model = matrix(c(1,1,0),ncol=3, byrow=TRUE))
#' zippoisson.fit<-saemix(rapimod.zip,saemix.data,saemix.options)
#' 
#' # Simulation functions to simulate from the models
#' saemix.simulatePoisson<-function(psi, id, xidep) {
#'   time<-xidep[,1]
#'   y<-xidep[,2]
#'   intercept<-psi[id,1]
#'   slope<-psi[id,2]
#'   lambda<- exp(intercept + slope*time)
#'   y<-rpois(length(time), lambda=lambda)
#'   return(y)
#' }
#' saemix.simulatePoissonZIP<-function(psi, id, xidep) {
#'   time<-xidep[,1]
#'   y<-xidep[,2]
#'   intercept<-psi[id,1]
#'   slope<-psi[id,2]
#'   p0<-psi[id,3] 
#'   lambda<- exp(intercept + slope*time)
#'   prob0<-rbinom(length(time), size=1, prob=p0)
#'   y<-rpois(length(time), lambda=lambda)
#'   y[prob0==1]<-0
#'   return(y)
#' }
#' 
#' # Using simulations to compare the predicted proportion of 0's in the two models
#' nsim<-100
#' yfit1<-simulateDiscreteSaemix(poisson.fit, saemix.simulatePoisson, nsim=nsim)
#' yfit2<-simulateDiscreteSaemix(zippoisson.fit, saemix.simulatePoissonZIP, 100)
#' {
#' nobssim<-length(yfit1@sim.data@datasim$ysim)
#' cat("Observed proportion of 0's", 
#'    length(yfit1@data@data$rapi[yfit1@data@data$rapi==0])/yfit1@data@ntot.obs,"\n")
#' cat("      Poisson model, p=",
#'    length(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim==0])/nobssim,"\n")
#' cat("          ZIP model, p=",
#'    length(yfit2@sim.data@datasim$ysim[yfit2@sim.data@datasim$ysim==0])/nobssim,"\n")
#' }
#'   }
#'   
#' @keywords datasets
NULL


