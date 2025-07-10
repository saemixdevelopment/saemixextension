####################################################################################
####			SaemixData class - definition				####
####################################################################################

#' @include aaa_generics.R
NULL

###############################
# ECO TODO: check on name validity

#' Class "SaemixData"
#' 
#' An object of the SaemixData class, representing a longitudinal data structure, used by the SAEM algorithm.
#' 
#' @name SaemixData-class 
#' @docType class
#' @aliases SaemixData SaemixData-class 
#' @aliases print,SaemixData showall,SaemixData show,SaemixData
#' @aliases SaemixRepData-class SaemixRepData 
#' @aliases SaemixSimData-class SaemixSimData 
#' 
#' @section Objects from the Class: 
#' An object of the SaemixData class can be created by using the function \code{\link{saemixData}} and contain the following slots:
#' @slot name.data Object of class \code{"character"}: name of the dataset
#' @slot header Object of class \code{"logical"}: whether the dataset/file contains a header. Defaults to TRUE 
#' @slot sep Object of class \code{"character"}: the field separator character
#' @slot na Object of class \code{"character"}: a character vector of the strings which are to be interpreted as NA values
#' @slot messages Object of class \code{"logical"}: if TRUE, the program will display information about the creation of the data object
#' @slot automatic Object of class \code{"logical"}: if TRUE, automatic name recognition is on (used at the creation of the object)
#' @slot name.group Object of class \code{"character"}: name of the column containing the subject id
#' @slot name.predictors Object of class \code{"character"}: name of the column(s) containing the predictors
#' @slot name.response Object of class \code{"character"}: name of the column containing the response variable y modelled by predictor(s) x
#' @slot name.covariates Object of class \code{"character"}: name of the column(s) containing the covariates, if present (otherwise empty)
#' @slot name.X Object of class \code{"character"}: name of the column containing the regression variable to be used on the X axis in the plots
#' @slot name.mdv Object of class \code{"character"}: name of the column containing the indicator variable denoting missing data
#' @slot name.cens Object of class \code{"character"}: name of the column containing the indicator variable denoting censored data (the value in the name.response column will be taken as the censoring value)
#' @slot name.occ Object of class \code{"character"}: name of the column containing the value of the occasion
#' @slot name.ytype Object of class \code{"character"}: name of the column containing the response number
#' @slot trans.cov Object of class \code{"list"}: the list of transformations leading to the new covariate # note: there could be a covariate class
#' @slot units Object of class \code{"list"}: list with up to three elements, x, y and optionally covariates, containing the units for the X and Y variables respectively, as well as the units for the different covariates
#' @slot data Object of class \code{"data.frame"}: dataframe containing the data, with columns for id (name.group), predictors (name.predictors), response (name.response), and covariates if present in the dataset (name.covariates). A column "index" contains the subject index (used to map the subject id). The column names, except for the additional column index, correspond to the names in the original dataset.
#' @slot N Object of class \code{"numeric"}: number of subjects
#' @slot yorig Object of class \code{"numeric"}: response data, on the original scale. Used when the error model is exponential
#' @slot ocov Object of class \code{"data.frame"}: original covariate data (before transformation in the algorithm)
#' @slot ind.gen Object of class \code{"logical"}: indicator for genetic covariates (internal)
#' @slot ntot.obs Object of class \code{"numeric"}: total number of observations
#' @slot nind.obs Object of class \code{"numeric"}: vector containing the number of observations for each subject
#' @section Methods:
#'   \describe{
#'     \item{[<-}{\code{signature(x = "SaemixData")}: replace elements of object}
#'     \item{[}{\code{signature(x = "SaemixData")}: access elements of object}
#'     \item{initialize}{\code{signature(.Object = "SaemixData")}: internal function to initialise object, not to be used}
#'     \item{plot}{\code{signature(x = "SaemixData")}: plot the data}
#'     \item{print}{\code{signature(x = "SaemixData")}: prints details about the object (more extensive than show)}
#'     \item{read}{\code{signature(object = "SaemixData")}: internal function, not to be used }
#'     \item{showall}{\code{signature(object = "SaemixData")}: shows all the elements in the object}
#'     \item{show}{\code{signature(object = "SaemixData")}: prints details about the object}
#'     \item{summary}{\code{signature(object = "SaemixData")}: summary of the data. Returns a list with a number of elements extracted from the dataset (N: the number of subjects; nobs: the total number of observations; nind.obs: a vector giving the number of observations for each subject; id: subject ID; x: predictors; y: response, and, if present in the data, covariates: the covariates (as many lines as observations) and ind.covariates: the individual covariates (one line per individual).}
#'     \item{subset}{\code{signature(object = "SaemixData")}: extract part of the data; this function will operate on the rows of the dataset (it can be used for instance to extract the data corresponding to the first ten subjects)}
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
#' showClass("SaemixData")
#' 
#' # Specifying column names
#' data(theo.saemix)
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' 
#' # Specifying column numbers
#' data(theo.saemix)
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
#'   name.group=1,name.predictors=c(2,3),name.response=c(4), name.covariates=5:6, 
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' 
#' # No column names specified, using automatic recognition of column names
#' data(PD1.saemix)
#' saemix.data<-saemixData(name.data=PD1.saemix,header=TRUE, 
#'   name.covariates=c("gender"),units=list(x="mg",y="-",covariates=c("-")))
#' 
#' @keywords classes
#' @exportClass SaemixData

setClass(
  Class="SaemixData",
  representation=representation(
    name.data="character",	# name of dataset
    header="logical",		# for file, whether has header
    messages="logical",		# whether to print messages when creating the object
    automatic="logical",		# whether to use automatic name recognition when creating the object
    sep="character",		# if file, separator
    na="character",		# if file, NA symbol(s)
    name.group="character",	# name of column with ID
    name.predictors="character",# name of column(s) with predictors 
    name.response="character",	# name of column with response
    name.covariates="character",# name of column(s) with covariates
    name.X="character",		# name of predictor used on X axis for graphs
    name.mdv="character", # name of column indicating a missing response
    name.cens="character", # name of column indicating a censored response
    name.occ="character", # name of column with the occasion
    name.ytype="character", # name of column with type of response (1,2,...)
    trans.cov="list",	# a list of transformations applied to the covariates
    units="list",		# units (list with components for x, y, and cov)
    data="data.frame",		# the data (data frame with columns name.group (subject id), index (id renamed to 1:N), name.predictors (predictors), name.response (possibly transformed during fit), name.covariates), mdv (missing data), cens (censored data, 1=censored & value in column response is the LOQ, ytype (type of the response), occ (occasion)); binary covariates are modified to 0/1
    ocov="data.frame",		# original scale for the covariates
    N="numeric",		# number of subjects
    yorig="numeric",		# vector of responses in original dataset
    ind.gen="logical",	# vector of booleans (same size as name.covariates); TRUE=genetic covariate, FALSE=non-genetic covariates
    ntot.obs="numeric",		# total number of observations (=dim(tab)[1])
    nind.obs="numeric"		# number of observations for each subject
  ),
  validity=function(object){
    #    cat ("--- Checking SaemixData object ---\n")
    if (length(object@name.data)==0) {
      if(object@messages) message("[ SaemixData : validation ] Please provide a name for the data (dataset or datafile on disk).")
      return("No dataset provided")
    }
    # Ici ou a la creation, detection automatique ?
    misVar<-0
    errors<-c()
    if ( length(object@name.group)==0) {
      if(object@messages) message("Missing Id column and automatic recognition is off\n")
      errors<-c(errors,"Missing Id column")
      misVar<-1
    }
    if (length(object@name.predictors)==0) {
      if(object@messages) message("No predictors found and automatic recognition is off\n")
      errors<-c(errors,"Missing predictors")
      misVar<-1
    }
    if (length(object@name.response)==0) {
      if(object@messages) message("No response found and automatic recognition is off\n")
      errors<-c(errors,"Missing response")
      misVar<-1
    }
    if(misVar==1) {
      if(object@messages) message("[ SaemixData : validation ] At least one of group, predictor or response variable name is missing, and automatic recognition is off.")
    }
    if(length(errors)==0) {
    N<-object@N
    if(length(object@data)>0) {
     if(N<2) {
          if (object@messages) message("Warning: There is only",N,"subject in the dataset, the SAEM algorithm is a population algorithm designed to analyse longitudinal data from non-linear mixed effect models and may not work with too few subjects.\n")
            }
          }
    }

    if(length(errors)==0) return(TRUE) else return(errors)
  }
)


#' @rdname SaemixData-class
#' @exportClass SaemixSimData

setClass(
  Class="SaemixRepData", # Saemix data, replicated for different chains
  representation=representation(
    N="numeric",		# number of subjects
    NM="numeric",		# number of subjects, replicated
    dataM="data.frame",		# replicated data with columns IdM, xM, yM
#    IdM="numeric",		# subject id
#    XM="data.frame",		# matrix of predictors
#    yM="numeric",		# vector of responses 
    nrep="numeric"		# number of replicates
  ),
  validity=function(object){
#    cat ("--- Checking SaemixData object ---\n")
    return(TRUE)
  }
)

#' @rdname SaemixData-class
#' @exportClass SaemixRepData

setClass(
  Class="SaemixSimData", # Saemix predicted and simulated data
  representation=representation(
    N="numeric",		# number of subjects
    name.group="character", # name of column with ID element
    name.response="character",	# name of column with response
    name.predictors="character",# name of column(s) with predictors 
    name.X="character",		# name of predictor used on X axis for graphs
    units="list",		# units (list with components for x, y, and cov)
    data="data.frame",		# ECO TODO: do we need to keep it here ?
    nsim="numeric",		# number of simulations
    datasim="data.frame",	# simulated data with columns idsim (id in replications), irep (replication number), ypred (simulated predictions, without error), ysim (simulated data, with error)
    sim.psi="data.frame"	# simulated parameters
  ),
  validity=function(object){
#    cat ("--- Checking saemixSimData object ---\n")
    return(TRUE)
  }
)

###############################
# ECO validity ne semble pas etre appele automatiquement quand on cree un objet => il faut l'appeler dans initialize

#' @rdname initialize-methods
#' 
#' @param name.data name of the dataset (can be a character string giving the name of a file on disk or of a dataset in the R session, or the name of a dataset
#' @param header whether the dataset/file contains a header. Defaults to TRUE
#' @param sep the field separator character. Defaults to any number of blank spaces ("")
#' @param na a character vector of the strings which are to be interpreted as NA values. Defaults to c(NA)
#' @param name.group name (or number) of the column containing the subject id
#' @param name.predictors name (or number) of the column(s) containing the predictors (the algorithm requires at least one predictor x)
#' @param name.response name (or number) of the column containing the response variable y modelled by predictor(s) x
#' @param name.covariates name (or number) of the column(s) containing the covariates, if present (otherwise missing)
#' @param name.X name of the column containing the regression variable to be used on the X axis in the plots (defaults to the first predictor)
#' @param units list with up to three elements, x, y and optionally covariates, containing the units for the X and Y variables respectively, as well as the units for the different covariates (defaults to empty)
#' @param name.mdv name of the column containing the indicator for missing variable
#' @param name.cens name of the column containing the indicator for censoring
#' @param name.occ name of the column containing the occasion
#' @param name.ytype name of the column containing the index of the response
#' @param verbose a boolean indicating whether messages should be printed out during the creation of the object
#' @param automatic a boolean indicating whether to attempt automatic name recognition when some colum names are missing or wrong (defaults to TRUE)
#' 
#' @exportMethod initialize

setMethod(
  f="initialize",
  signature="SaemixData",
  definition= function (.Object,name.data,header,sep,na,name.group, name.predictors, name.response, name.covariates, name.X, units, name.mdv, name.cens, name.occ, name.ytype, verbose=TRUE, automatic=TRUE){
#    cat ("--- initialising SaemixData Object --- \n")
    if(missing(name.data)) stop ("Please provide a name for the data (dataset or datafile on disk).")
    .Object@name.data<-name.data
    if(missing(header)) header<-TRUE
    .Object@header<-header
    if(missing(verbose)) verbose<-TRUE
    .Object@messages<-verbose
    .Object@automatic<-automatic
    if(missing(sep)) sep<-""
    .Object@sep<-sep
    if(missing(na)) na<-"NA"
    .Object@na<-na
    if(missing(name.group)) {
      if(automatic) {
      if(verbose) cat("   Missing ID identifier, assuming the ID is in column 1 of the dataset.\n")
      name.group<-"1"
      } else name.group<-character()
    }
# ECO TODO: reconnaissance automatique (avant affectation a la valeur 2) ?
    if(missing(name.predictors)) {
      if(automatic) {
      name.predictors<-"2"      
      if(verbose) cat("   Missing predictors identifier, assuming there is one predictor in column 2 of the dataset.\n")
      } else name.predictors<-character()
    }
    if(missing(name.response)) {
      if(automatic) {
    	if(verbose) cat("   Missing response identifier, assuming the response is in column 3 of the dataset.\n")
      name.response<-"3"
      } else name.response<-character()
    }
    if(missing(name.covariates)) name.covariates<-character()
		if(missing(name.mdv)) name.mdv<-character()
		if(missing(name.cens)) name.cens<-character()
		if(missing(name.occ)) name.occ<-character()
		if(missing(name.ytype)) name.ytype<-character()
    if(missing(name.X)) name.X<-character()
		.Object@name.group<-name.group
    .Object@name.predictors<-name.predictors
    .Object@name.response<-name.response
    .Object@name.covariates<-name.covariates
		.Object@name.mdv<-name.mdv
		.Object@name.cens<-name.cens
		.Object@name.occ<-name.occ
		.Object@name.ytype<-name.ytype
		.Object@trans.cov<-list()
		if(missing(units)) units<-list(x="-",y="-")
    if(is.null(units$x)) units$x<-"-"
    if(is.null(units$y)) units$y<-"-"
    ncov<-length(name.covariates)
    if(ncov>0) {
      nunit<-length(units$covariates)
      if(nunit==0) units$covariates<-rep("-",ncov)
      if(nunit>ncov) units$covariates<-units$covariates[1:ncov]
      if(nunit<ncov) {
        length(units$covariates)<-ncov
        units$covariates[(nunit+1):ncov]<-"-"
      }
    }
    .Object@units<-units
    .Object@name.X<-name.X
# Object validation
#    validObject(.Object)
    return (.Object )
  }
)

# Initialize method for saemixRepData and saemixSimData
#' @rdname initialize-methods
#' 
#' @param nb.chains number of chains used in the algorithm
#' 
#' @exportMethod initialize

setMethod(
  f="initialize",
  signature="SaemixRepData",
  definition= function (.Object,data=NULL,nb.chains=1){
#    cat ("--- initialising SaemixData Object --- \n")
    if(is.null(data)) {
      .Object@N<-.Object@NM<-numeric(0)
      .Object@dataM<-data.frame()
    } else {
    N<-data@N
    .Object@N<-N
    .Object@nrep<-nb.chains
    .Object@NM<-N*nb.chains
    IdM<-kronecker(c(0:(nb.chains-1)),rep(N,data@ntot.obs))+rep(data@data[,"index"], nb.chains)
    yM<-rep(data@data[,data@name.response],nb.chains)
    XM<-do.call(rbind,rep(list(data@data[,c(data@name.predictors,data@name.cens,data@name.mdv,data@name.ytype),drop=FALSE]), nb.chains))
    .Object@dataM<-data.frame(IdM=c(IdM),XM,yM=yM)
   }
# Object validation
#    validObject(.Object)
    return (.Object )
  }
)

#' @rdname initialize-methods
#' 
#' @param datasim dataframe containing the simulated data
#' 
#' @exportMethod initialize
#' 

## Warning: currently won't work with RTTE (not the same dimension for data and for datasim, there can be different number of events...)
setMethod(
  f="initialize",
  signature="SaemixSimData",
  definition= function (.Object,data=NULL,datasim=NULL) {
#    cat ("--- initialising SaemixData Object --- \n")
    if(!is.null(data)) {
      .Object@N<-data@N
      .Object@name.group<-data@name.group
      .Object@name.response<-data@name.response
      .Object@name.predictors<-data@name.predictors
      .Object@name.X<-data@name.X
      .Object@units<-data@units
      .Object@data<-data@data
    }
    if(is.null(data) || is.null(datasim) || dim(datasim)[1]==0) {
      .Object@datasim<-data.frame()
      .Object@nsim<-0
    } else {
      .Object@datasim<-datasim
      .Object@nsim<-dim(datasim)[1]/dim(data@data)[1] # not for RTTE, will need to change this by hand in the simulations
    }
# Object validation
#    validObject(.Object)
    return (.Object )
  }
)

####################################################################################
####			SaemixData class - accesseur				####
####################################################################################

#' Get/set methods for SaemixData object
#' 
#' Access slots of an SaemixData object using the object\["slot"\] format
#' 
#' @name [
#' @aliases [<-,SaemixData-method [,SaemixData-method
#' 
#' @param x object
#' @param i element to be replaced
#' @param j element to replace with
#' @param value value to replace with
#' @param drop whether to drop unused dimensions
#' @docType methods
#' @keywords methods
#' @exportMethod [
#' @exportMethod [<-
#' @exportPattern "^[[:alpha:]]+"
#' @rdname extract-methods


# Getteur
setMethod(
  f ="[",
  signature = "SaemixData" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "name.data"={return(x@name.data)},
    "header"={return(x@header)},
    "messages"={return(x@messages)},
    "automatic"={return(x@automatic)},
    "sep"={return(x@sep)},
    "na"={return(x@na)},
    "name.group"={return(x@name.group)},
    "name.predictors"={return(x@name.predictors)},
    "name.response"={return(x@name.response)},
    "name.covariates"={return(x@name.covariates)},
    "name.X"={return(x@name.X)},
    "name.mdv"={return(x@name.mdv)},
    "name.cens"={return(x@name.cens)},
    "name.occ"={return(x@name.occ)},
    "name.ytype"={return(x@name.ytype)},
    "trans.cov"={return(x@trans.cov)},    
    "units"={return(x@units)},
    "data"={return(x@data)},
    "ocov"={return(x@ocov)},
    "N"={return(x@N)},
    "yorig"={return(x@yorig)},
    "ind.gen"={return(x@ind.gen)},
    "ntot.obs"={return(x@ntot.obs)},
    "nind.obs"={return(x@nind.obs)},
    stop("No such attribute\n")
   )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixData" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "name.data"={x@name.data<-value},
    "header"={x@header<-value},
    "messages"={x@messages<-value},
    "automatic"={x@automatic<-value},
    "sep"={x@sep<-value},
    "na"={x@na<-value},
    "name.group"={x@name.group<-value},
    "name.predictors"={x@name.predictors<-value},
    "name.response"={x@name.response<-value},
    "name.covariates"={x@name.covariates<-value},
    "name.X"={x@name.X<-value},
    "name.mdv"={x@name.mdv<-value},
    "name.cens"={x@name.cens<-value},
    "name.occ"={x@name.occ<-value},
    "name.ytype"={x@name.ytype<-value},
    "trans.cov"={x@trans.cov<-value},
    "units"={x@units<-value},
    "data"={x@data<-value},
    "ocov"={x@ocov<-value},
    "N"={x@N<-value},
    "ind.gen"={x@ind.gen<-value},
    "yorig"={x@yorig<-value},
    "ntot.obs"={x@ntot.obs<-value},
    "nind.obs"={x@nind.obs<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

# For saemixRepData

#' extract parts of SaemixRepData
#'
#' @name [
#' @aliases [,SaemixRepData-method
#' @docType methods
#' @rdname extract-methods
#' 
setMethod(
  f ="[",
  signature = "SaemixRepData" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "N"={return(x@N)},
    "NM"={return(x@NM)},
    "dataM"={return(x@dataM)},
    "nrep"={return(x@nrep)},
    stop("No such attribute\n")
   )
  }
)

#' replace names of SaemixRepData
#'
#' @name [
#' @aliases [<-,SaemixRepData-method
#' @docType methods
#' @rdname extract-methods

setReplaceMethod(
  f ="[",
  signature = "SaemixRepData" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "N"={x@N<-value},
    "NM"={x@NM<-value},
    "dataM"={x@dataM<-value},
    "nrep"={x@nrep<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

# For saemixSimData

#' extract parts of SaemixSimData
#'
#' @name [
#' @aliases [,SaemixSimData-method
#' @docType methods
#' @rdname extract-methods
#' 

setMethod(
  f ="[",
  signature = "SaemixSimData" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "N"={return(x@N)},
    "name.group"={return(x@name.group)},
    "name.response"={return(x@name.response)},
    "name.predictors"={return(x@name.predictors)},
    "name.X"={return(x@name.X)},
    "units"={return(x@units)},
    "data"={return(x@data)},
    "nsim"={return(x@nsim)},
    "sim.psi"={return(x@sim.psi)},
    "datasim"={return(x@datasim)},
    stop("No such attribute\n")
   )
  }
)

#' replace names of SaemixSimData
#'
#' @name [
#' @aliases [<-,SaemixSimData-method
#' @docType methods
#' @rdname extract-methods

setReplaceMethod(
  f ="[",
  signature = "SaemixSimData" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "N"={x@N<-value},
    "name.group"={x@name.group<-value},
    "name.response"={x@name.response<-value},
    "name.predictors"={x@name.predictors<-value},
    "name.X"={x@name.X<-value},
    "units"={x@units<-value},
    "data"={x@data<-value},
    "sim.psi"={x@sim.psi<-value},
    "datasim"={x@datasim<-value},
    "nsim"={x@nsim<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

####################################################################################
