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
#' @slot trans.cov Object of class \code{"list"}: the list of transformation applied to the covariates (currently unused, TODO)
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
      .Object@nsim<-dim(datasim)[1]/dim(data@data)[1]
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
#' Access slots of a SaemixData object using the object\["slot"\] format
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
####			SaemixData class - method to read data			####
####################################################################################

#' Name validation (## )Helper function not intended to be called by the user)
#' 
#' Helper function, checks if the names given by the user match to the names in the dataset. If not, automatic recognition is attempted when automatic=TRUE.
#'
#' @name validate.names
#' 
#' @param usernames vector of strings
#' @param datanames vector of strings
#' @param recognisednames vector of strings, values for automatic recognition
#' @param verbose logical, whether to print warning messages
#' @param automatic a boolean indicating whether to attempt automatic name recognition when some colum names are missing or wrong (defaults to TRUE)
#' @return a vector with valid names
#' @examples 
#' # TODO
#' @keywords methods
#' @export 
NULL

validate.names<-function(usernames,datanames,recognisednames=c(),verbose=TRUE, automatic=TRUE) {
  valnames<-usernames
  remcol<-c() # keep track of columns to remove
  # Detect names given as column numbers
  convnames<-suppressWarnings(as.integer(usernames))
  icol.int<-which(!is.na(convnames))
  if(length(icol.int)>0) {
    namcol.int<-suppressWarnings(as.integer(usernames[icol.int]))
    ioutrange<-namcol.int[namcol.int>length(datanames) | namcol.int<0]
    if(length(ioutrange)>0) {
      if(verbose) cat("Column number(s)",ioutrange,"do(es) not exist in the dataset, please check\n")
      remcol<-icol.int[namcol.int %in% ioutrange]
      icol.int<-icol.int[!namcol.int %in% ioutrange]
      namcol.int<-namcol.int[!namcol.int %in% ioutrange]
    }
    valnames[icol.int]<-datanames[namcol.int]
  }
  # Validate name columns given as strings
  icol.str<-which(is.na(convnames))
  if(length(icol.str)>0) {
    namcol.str<-usernames[icol.str]
    inotfound<-which(is.na(match(namcol.str,datanames)))
    if(length(inotfound)>0) {
      remcol<-c(remcol,icol.str[inotfound])
      if(verbose) cat("Column name(s)",namcol.str[inotfound],"do(es) not exist in the dataset, please check\n")
    }
  }
  if(length(remcol)>0) {
    if(verbose) cat("Remove columns",remcol,"(",usernames[remcol],")","\n")
    valnames<-valnames[-remcol]
  }
  # Automatic column detection (only if no name given or none recognised)
  if(length(valnames)==0) {
    if(length(recognisednames)==0) {
      if(verbose) cat("Please check input\n")
      return(valnames)
    } 
    if(automatic) {
    if(verbose) cat("No valid name given, attempting automatic recognition\n")
    valnames<-datanames[match(recognisednames,tolower(datanames),nomatch=0)]
    if(length(valnames)==0) {
      if(verbose) cat("Automatic recognition failed, please check input\n")} else {
        if(verbose) cat("Automatic recognition of columns",valnames,"successful \n")
      }
    } else return(NULL)
  }
  # Remove columns
  return(valnames)
}

## Create a longitudinal data structure from a file or a dataframe
## Helper function not intended to be called by the user
## @exportMethod read

#' Create a longitudinal data structure from a file or a dataframe
#'  
#' Helper function not intended to be called by the user
#'  
#' @param object an SaemixData object
#' @param dat the name of a dataframe in the R environment, defaults to NULL; if NULL, the function will
#' attempt to read the file defined by the slot name.data.
#'  
#' @rdname read-methods
#' @aliases readSaemix,SaemixData readSaemix,SaemixData-method
#'  
#' @exportMethod readSaemix

setMethod("readSaemix",
          signature="SaemixData",
  function(object, dat = NULL) {
    ow <- options("warn")
    options("warn"=-1)
# ce test devrait aller dans la definition de la classe
    if(!is(object@name.data,"character")) {
      if(object@messages) message("Please provide the name of the data (data.frame or path to file on disk) as a character string.\n")
    return("Creation of SaemixData object failed")
  }
    if(is.null(dat)) {
      if(object@messages) cat("Reading data from file",object@name.data,"\n")
      header<-object@header
      if(is.null(header)) header<-TRUE
      sep<-object@sep
      if(is.null(sep)) sep<-""
      na.strings<-object@na
      if(is.null(na.strings)) na.strings<-"NA"
      dat<-try(read.table(object@name.data,header=header,sep=sep,na.strings=na.strings))
      if(inherits(dat,"try-error")) stop("The file ",object@name.data," does not exist. Please check the name and path.\n")      
      if(object@messages) {
        cat("These are the first lines of the dataset as read into R. Please check the format of the data is appropriate, if not, modify the na and/or sep items and retry:\n")
        print(head(dat))
      }
    }
    if(dim(dat)[2]<3) {
      if(object@messages) message("The dataset does not contain enough data. The non-linear mixed effect model requires at least 3 columns, with subject ID, predictor (at least one) and response. \nPlease check the field separator, currently given as:", paste("sep=\"",object@sep,"\"",sep=""),"\n")
      return("Creation of SaemixData object failed")
    }
# Automatic recognition of columns 
#    ID (one of id, subject, sujet, group, groupe regardless of case)
#    response (one of Y, conc, concentration, resp, response, y, dv regardless of case)
#    predictors (time and/or dose, x, regardless of case)
# ECO TODO: improve automatic recognition ?
# check that we have at least a column id, response and X

    vnames<-validate.names(object@name.group,colnames(dat),recognisednames=c("id","subject","sujet","group","groupe"),verbose = object@messages, automatic=object@automatic)
    if(length(vnames)==0) {
      if(object@messages) message("Please provide a valid name for the ID column.\n")
      return("Creation of SaemixData object failed")
    }
    object@name.group<-vnames
    
    vnames<-validate.names(object@name.predictors,colnames(dat),recognisednames=c("time","temps","tps","tim","x","dose"),verbose = object@messages, automatic=object@automatic)
    if(length(vnames)==0) {
      if(object@messages) message("Please provide a valid name for the predictor(s).\n")
      return("Creation of SaemixData object failed")
    }
    object@name.predictors<-vnames
    
    vnames<-validate.names(object@name.response,colnames(dat),recognisednames=c("response","resp","conc","concentration","y","dv"),verbose = object@messages, automatic=object@automatic)
    if(length(vnames)==0) {
      if(object@messages) message("Please provide a valid name for the response.\n")
      return("Creation of SaemixData object failed")
    }
    object@name.response<-vnames[1]
    if(length(vnames)>1 & object@messages) cat("Using the response",object@name.response,"as dependent variable.\n")
    
		if(length(object@name.covariates)>0) {
    	if(object@name.covariates[1]!="") {
  		i1<-suppressWarnings(as.integer(object@name.covariates[!is.na(suppressWarnings(as.integer(object@name.covariates)))]))
      object@name.covariates[!is.na(suppressWarnings(as.integer(object@name.covariates)))]<- colnames(dat)[i1]
    	}
  		idx<-object@name.covariates[!(object@name.covariates %in% colnames(dat))]
  		if(length(idx)>0) {
  		  if(object@messages) cat("Covariates",object@name.covariates[idx],"not found.\n") 
  			object@units$covariates<-object@units$covariates[object@name.covariates %in% colnames(dat)]
  			object@name.covariates<-object@name.covariates[object@name.covariates %in% colnames(dat)]
  		}
    }
    if(nchar(object@name.group)*length(object@name.predictors)*nchar(object@name.response)<=0) {
      stop("Please check the structure of the data file and provide information concerning which columns specify the group structure (ID), the predictors (eg dose, time) and the response (eg Y, conc). See documentation for automatic recognition of column names for these elements.\n")
    }
	if(nchar(object@name.X)==0)
      object@name.X<-object@name.predictors[1]
    if(!is.na(suppressWarnings(as.integer(object@name.X)))) {
      if(dim(dat)[2]<suppressWarnings(as.integer(object@name.X))) {
        if(object@messages) cat("Attribute name.X",object@name.X,"does not correspond to a valid column in the dataset, setting the X axis for graphs to",object@name.predictors[1],".\n")
	object@name.X<-object@name.predictors[1]
      } else object@name.X<-colnames(dat)[suppressWarnings(as.integer(object@name.X))]
    } 
    if(match(object@name.X,object@name.predictors,nomatch=0)==0) {
      if(object@messages) cat("Attribute name.X",object@name.X,"does not correspond to a valid column in the dataset, setting the X axis for graphs to",object@name.predictors[1],".\n")
      object@name.X<-object@name.predictors[1]
    }
		if(nchar(object@name.mdv)==0) mdv<-rep(0,dim(dat)[1]) else {mdv<-dat[,object@name.mdv]}
		mdv[is.na(dat[,object@name.response])]<-1
		#    if(sum(mdv)>0) 
		object@name.mdv<-"mdv"
		if(nchar(object@name.cens)==0) cens<-rep(0,dim(dat)[1]) else cens<-dat[,object@name.cens]
		object@name.cens <-"cens"
		if(nchar(object@name.occ)==0) occ<-rep(1,dim(dat)[1]) else occ<-dat[,object@name.occ]
		object@name.occ<-"occ"
		if(nchar(object@name.ytype)==0) ytype<-rep(1,dim(dat)[1]) else ytype<-dat[,object@name.ytype]
		object@name.ytype<-"ytype"
		all.names<-c(object@name.group,object@name.predictors, object@name.response,object@name.covariates)
		
    dat<-dat[,all.names,drop=FALSE]
		dat<-cbind(dat,mdv=mdv,cens=cens,occ=occ,ytype=ytype)

    if(!is(dat,"data.frame")) dat<-as.data.frame(dat)
# Saving covariates in the original format in ocov, transforming binary covariates in dat to factors
    object@ocov<-dat[,object@name.covariates,drop=FALSE]
    for(icov in object@name.covariates) {
      if(length(unique(dat[,icov]))==2) dat[,icov]<-suppressWarnings(as.integer(factor(dat[,icov])))-1
    }
# Removing missing values in predictor columns
# dat<-dat[!is.na(dat[,object@name.response]),]
	for(i in object@name.predictors) {
		if(sum(is.na(dat[,i]))>0) {
		  if(object@messages) cat("Removing missing values for predictor",i,"\n")
		  if(!is.null(dim(object@ocov)[1])) object@ocov<-object@ocov[!is.na(dat[,i]),,drop=FALSE]
  		dat<-dat[!is.na(dat[,i]),]
	}
	}
# Removing subjects with only MDV in responses
	idx<-c();inull<-c()
	for(isuj in unique(dat[,object@name.group])) {
		if(sum(1-dat$mdv[dat[,object@name.group]==isuj])==0) {
			inull<-c(inull,isuj)
			idx<-c(idx,which(dat[,object@name.group]==isuj))
		}
	}
#  if(object@messages) print(idx)
  if(length(inull)>0) {
    if(object@messages) cat("Some subjects have no observations, removing them:",inull,"\n")
  	dat<-dat[-idx,]
  	if(!is.null(dim(object@ocov)[1])) object@ocov<-object@ocov[-idx,,drop=FALSE]
  }
	if(!is.null(dim(object@ocov)[1])) object@ocov <- object@ocov[order(dat[,object@name.group], dat[,object@name.X]),,drop=FALSE]
	
# ECO TODO: missing data in covariates kept for the moment, only excluded depending on the model
#    for(i in object@name.covariates) dat<-dat[!is.na(dat[,i]),]
  	object@ntot.obs<-dim(dat)[1] # total number of observations
    dat <- dat[order(dat[,object@name.group], dat[,object@name.X]),]
    id<-dat[,object@name.group]
    object@N<-length(unique(id))
    nind.obs<-tapply(id,id,length) # individual numbers of observations (1xN)
    nind.obs<-nind.obs[match(unique(id),names(nind.obs))]
    object@nind.obs<-c(nind.obs)
    dat<-cbind(index=rep(1:object@N,times=nind.obs),dat)
    object@data<-dat
    
    options(ow) # reset
    validObject(object)
    return(object)
  }
)

####################################################################################
####			SaemixData class - method to print/show data		####
####################################################################################

#' @rdname print-methods
#' 
#' @param x an object of type SaemixData, SaemixModel, SaemixRes or SaemixObject
#' @param nlines maximum number of lines of data to print (defaults to 10)
#' @param ... additional arguments passed on the print function
#' 
#' @exportMethod print

setMethod("print","SaemixData",
  function(x,nlines=10,...) {
    digits<-2;nsmall<-2
    cat("Object of class SaemixData\n")
    cat("    longitudinal data for use with the SAEM algorithm\n")
    cat("Dataset",x@name.data,"\n")
    st1<-paste(x@name.response," ~ ",paste(x@name.predictors,collapse=" + ")," | ", x@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    if(length(x@name.predictors)>1) {
      cat("    X variable for graphs:",x@name.X,paste("(",x@units$x,")",sep=""),"\n")
    } else  cat("    Predictor:",x@name.X,paste("(",x@units$x,")",sep=""),"\n")
    ncov<-length(x@name.covariates)
    if(ncov>0) {
      cat("    covariates:",paste(paste(x@name.covariates," (",x@units$covariates,")",sep=""),collapse=", "),"\n")
      if(length(x@ocov)>0) {
      for(icov in 1:ncov) {
      if(is.factor(x@ocov[,icov]) | length(unique(x@ocov[,icov]))==2) cat("      reference class for covariate",x@name.covariates[icov],": ",levels(as.factor(x@ocov[,icov]))[1],"\n")
      }
      }
    }
    if(FALSE) {
      cat("    Group column:",x@name.group,"\n")
      cat("    Predictors:",x@name.predictors,"\n")
      cat("    X variable for graphs:",x@name.X,paste("(",x@units$x,")",sep=""),"\n")
      cat("    Response column:",x@name.response, paste("(",x@units$y,")",sep=""),"\n")
      cat("    Covariates:",x@name.covariates,"\n")
    }
    if(length(x@data)>0) {
      if(nlines==0) return()
      cat("Dataset characteristics:\n")
      cat("    number of subjects:    ",x@N,"\n")
      if(x@N>0) {
        cat("    number of observations:",x@ntot.obs,"\n")
        cat("    average/min/max nb obs:",format(mean(x@nind.obs),digits=digits, nsmall=nsmall), " / ", min(x@nind.obs)," / ",max(x@nind.obs),"\n")
#    if(length(x@data)>0) print(x@data)
      }
      xdat<-x@data
      if(length(x@ocov)>0) try(xdat[,x@name.covariates]<-x@ocov)
      if(nlines==(-1)) {
        cat("Data:\n")
        print(xdat)
      } else {
        cat("First",nlines,"lines of data:\n")
        nrowShow <- min (nlines , nrow(xdat ))
        print(xdat[1:nrowShow,-c(1)])
      }
    } else cat("No data.\n")
  }
)

#' @rdname show-methods
#' 
#' @param object an object of type SaemixData, SaemixModel, SaemixRes or SaemixObject
#' 
#' @exportMethod show

setMethod("show","SaemixData",
  function(object) {
    cat("Object of class SaemixData\n")
    cat("    longitudinal data for use with the SAEM algorithm\n")
    cat("Dataset",object@name.data,"\n")
    st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + ")," | ", object@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    if(length(object@name.predictors)>1) {
      cat("    X variable for graphs:",object@name.X, paste("(",object@units$x,")",sep=""),"\n")
    }
    ncov<-length(object@name.covariates)
    if(ncov>0) {
      cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
      if(length(object@ocov)>0) {
      for(icov in 1:ncov) {
      if(is.factor(object@ocov[,icov])) cat("      reference class for covariate",object@name.covariates[icov],": ",levels(object@ocov[,icov])[1],"\n")
      }
      if(length(object@data)>0) try(object@data[,object@name.covariates]<-object@ocov)
      }
    }
    if(length(object@data)>0) {
      if(object@N>0) cat(object@ntot.obs,"    observations in",object@N,"subjects\n")
      cat("First lines of data:\n")
      nrowShow <- min (10 , nrow(object@data ))
      print(object@data[1:nrowShow,-c(1)])
    } else cat("No data.\n")
  }
)

#' @rdname showall-methods
#' @aliases showall
#' @exportMethod showall

# Could be print, with only head of data
setMethod("showall",signature="SaemixData",
  function(object) {
    digits<-2;nsmall<-2
    cat("Object of class SaemixData\n")
    cat("    longitudinal data for use with the SAEM algorithm\n")
    cat("Dataset",object@name.data,"\n")
    cat("    header:",object@header,"\n")
    cat("    sep:",object@sep,"\n")
    cat("    na:",object@na,"\n")
    st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + ")," | ", object@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    cat("    subject identifier:    ",object@name.group,"\n")
    cat("    predictors:       ",object@name.predictors,"\n")
    cat("    response:         ",object@name.response,paste("(",object@units$y,")",sep=""),"\n")
    cat("    X variable for graphs:",object@name.X,paste("(",object@units$x,")",sep=""),"\n")
    ncov<-length(object@name.covariates)
    if(ncov>0) {
      cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
      if(length(object@ocov)>0) {
      for(icov in 1:ncov) {
      if(is.factor(object@ocov[,icov])) cat("      reference class for covariate",object@name.covariates[icov],": ",levels(object@ocov[,icov])[1],"\n")
      }
      if(length(object@data)>0) try(object@data[,object@name.covariates]<-object@ocov)
      }
    }
    cat("Dataset characteristics:\n")
    cat("    number of subjects:    ",object@N,"\n")
    if(object@N>0) {
      cat("    number of observations:",object@ntot.obs,"\n")
      cat("    average/min/max nb obs:",format(mean(object@nind.obs),digits=digits, nsmall=nsmall), " / ", min(object@nind.obs)," / ",max(object@nind.obs),"\n")
#    if(length(object@data)>0) print(object@data)
    }
    if(length(object@data)>0) {
      cat("First lines of data:\n")
      nrowShow <- min (10 , nrow(object@data ))
      ncolShow <- min (10 , ncol(object@data))
      print(object@data[1:nrowShow,-c(1)])
    } else cat("No data.\n")
  }
)

# SaemixRepData
#' @rdname show-methods
#' @exportMethod show

setMethod("show","SaemixRepData",
  function(object) {
    cat("Object of class saemixRepData\n")
    if(length(object@N)>0) {
	    cat("    replicated data used in the SAEM algorithm\n")
	    cat("    number of subjects in initial dataset",object@N,"\n")
	    cat("    number of replications",object@nrep,"\n")
	    cat("    number of subjects in replicated dataset",object@NM,"\n")
    } else cat("Empty object \n")
    } 
)
     
# SaemixSimData
#' @rdname show-methods
#' @exportMethod show

setMethod("show","SaemixSimData",
  function(object) {
    cat("Object of class SaemixSimData\n")
    cat("    data simulated according to a non-linear mixed effect model\n")
    if(length(object@N)>0) {
    cat("Characteristics of original data\n")
    cat("    number of subjects:",object@N,"\n")
    cat("    summary of response:\n")
    print(summary(object@data[,object@name.response]))
    cat("Characteristics of simulated data\n")
    if(dim(object@datasim)[1]>0) {
      cat("    number of simulated datasets:",object@nsim,"\n")
      cat("    summary of simulated response\n")
      print(summary(object@datasim$ysim))
    } else cat("    no simulations performed yet\n")
  }}
)

####################################################################################
####				Summary method for SaemixData			####
####################################################################################

#' Summarising longitudinal data
#' 
#' summary method for class SaemixData
#' @name summary-methods
#' @aliases summary summary,SaemixData summary,SaemixData-method
#' 
#' @param object an object of class SaemixData
#' @param print a boolean controlling whether to print the output or return it silently
#' @param ... additional arguments (ignored)
#'     
#' @return a list with a number of elements extracted from the dataset
#' \describe{
#' \item{N}{ number of subjects}
#' \item{nobs}{ the total number of observations} 
#' \item{nind.obs}{a vector giving the number of observations for each subject}
#' \item{id}{subject ID; x: predictors; y: response, and, if present in the data, covariates: the covariates (as many lines as observations) and ind.covariates: the individual covariates (one}
#' }
#' @exportMethod summary

setMethod("summary","SaemixData",
  function(object, print=TRUE, ...) {
    digits<-2;nsmall<-2
    if(print) {
    	cat("Object of class SaemixData\n")
      cat("    longitudinal data for use with the SAEM algorithm\n")
      cat("Dataset",object@name.data,"\n")
      st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + ")," | ", object@name.group,sep="")
    	cat("    Structured data:",st1,"\n")
      if(length(object@name.predictors)>1) cat("    X variable for graphs:",object@name.X,paste("(",object@units$x,")",sep=""),"\n")
      if(length(object@name.covariates)>0) {
        cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
      }
      cat("Dataset characteristics:\n")
      cat("    number of subjects:    ",object@N,"\n")
      if(object@N>0) {
        cat("    number of observations:",object@ntot.obs,"\n")
        cat("    average/min/max nb obs:",format(mean(object@nind.obs),digits=digits, nsmall=nsmall), " / ", min(object@nind.obs)," / ",max(object@nind.obs),"\n")
#    if(length(object@data)>0) print(object@data)
      }
    }
    res<-list(N=object@N,nobs=list(ntot=object@ntot.obs,nind=object@nind.obs), id=object@data[,object@name.group],x=object@data[,object@name.predictors,drop=FALSE], y=object@data[,object@name.response])
    if(length(object@name.covariates)>0) {
      res$covariates<-object@ocov
      if(dim(object@data)[1]==dim(object@ocov)[1]) {
        ucov<-cbind(object@data[,object@name.group],object@ocov)
        colnames(ucov)[1]<-object@name.group
        ucov<-ucov[match(unique(object@data$index),object@data$index),]
        res$ind.covariates<-ucov
      }
    }
    invisible(res)
 }
)

####################################################################################
####			SaemixData class - method to plot			####
####################################################################################
### #' @export
### #' @docType methods
### #' @rdname plot-methods
### #' @aliases plot,SaemixData 


saemix.data.setoptions<-function(saemix.data) {
# setting default plot options
  plot.opt<-list(
# General graphical options
    new=TRUE,				# whether a new page should be called
    ask=FALSE,				# whether the program should ask before creating a new page
    ilist=c(1:saemix.data["N"]),
    separate=FALSE,	# if TRUE, plots individual subjects (a la nlme), if FALSE plots a single plot with all the subjects
# Options for individual plots
    nmax=12,					  # maximum number of subjects
    limit=TRUE,					# limit to nmax plots
    sample=FALSE,				# if FALSE=use the (nmax) first subjects; TRUE=randomly sample (nmax) subjects from the dataset
    interactive=FALSE,  # whether the program should prompt the user for the number of subjects to plot in the individual plots if this number exceeds nmax
    which.cov="none",  # whether to split over covariates
# Layout and plots options
    mfrow=c(),				# page layout (if empty, defaults to the default layout for each graph type)
    main=" ",				# title
    xlab=" ",
    ylab=" ",
    col="black",
    pch=20,
    lty=1,
    lwd=1,
    xlim=c(),
    ylim=c(),
    xlog=FALSE,
    ylog=FALSE,
    type="b",
    cex=1,
    cex.axis=1,
    cex.lab=1,
    cex.main=1)
        
     if(is.null(plot.opt$name.X))
        plot.opt$name.X<-saemix.data["name.predictors"][1]
    plot.opt$xlab<-paste(plot.opt$name.X," (",plot.opt$units$x,")", sep="")
     if(length(saemix.data["name.response"])>0)
    plot.opt$ylab<-paste(saemix.data["name.response"]," (",saemix.data["units"]$y,")", sep="")
   return(plot.opt)
}

replace.data.options<-function(plot.opt,...) {
  args1<-match.call(expand.dots=TRUE)
  # These arguments are used by other functions and may be passed on via "...", so we want to ignore them. Other arguments not in list will raise warnings
  legacy<-c("plot.type","individual")
  if(length(args1)>2) {
# Other arguments
    for(i in 3:length(args1)) {
      if(match(names(args1)[i],names(plot.opt),nomatch=0)>0) {
        #    plot.opt[[names(args1)[i]]]<-args1[[i]] else {
        if(!is.null(eval(args1[[i]]))) plot.opt[[names(args1)[i]]]<-eval(args1[[i]])
      } else {
        if(is.na(match(names(args1)[i],legacy))) message(paste("Argument",names(args1)[i],"not available, check spelling"))
      }
    }
  }
  return(plot.opt)
}

#' Plot of longitudinal data 
#' 
#' This function will plot a longitudinal dataframe contained in an SaemixData object. By default it produces a spaghetti plot, but arguments can be passed on to modify this behaviour. 
#' 
## #' @name plot-SaemixData
#' 
#' @param x an SaemixData object or an SaemixSimData object
#' @param y unused, present for compatibility with base plot function
#' @param ... additional arguments to be passed on to plot (titles, legends, ...)
#' 
#' @aliases plot,SaemixData-methods 
#' @aliases plot-SaemixData
#' @aliases plot,SaemixData plot,SaemixData,ANY-method
#' @keywords plot
### #' @docType methods
#' @rdname plot-SaemixData
#' 
#' @import ggplot2 grid gridExtra
#' @importFrom graphics plot
#' @importFrom rlang is_missing
#' @method plot SaemixData
#' @export 


# Plot the data, either as points or as lines grouped by x@name.group
plot.SaemixData<-function(x,y,...) {
  if(length(x@data)==0) {
    message("No data to plot.\n")
    return("Missing data")
  }
  # Eco: commented, otherwise was resetting the graphical layout on exit and preventing the graphs to be set on the same page
  #    oldpar <- par(no.readonly = TRUE)    # code line i
  #    on.exit(par(oldpar))            # code line i + 1
  
  # User-defined options
  userPlotOptions<-list(...)
  if(!is_missing(y) && is.list(y)) {
    userPlotOptions<-c(y,userPlotOptions)
  }
  i1<-match("individual",names(userPlotOptions))
  if(!is.na(i1)) {
    individual<-as.logical(eval(userPlotOptions[[i1]]))
  } else individual<-FALSE
  i1<-match("type",names(userPlotOptions))
  if(!is.na(i1)) {
    plot.type<-as.character(userPlotOptions[[i1]])
    plot.type<-plot.type[plot.type!="c"]
  } else plot.type<-c()
  if(length(plot.type)==0) plot.type<-ifelse(individual,"b","l")
  
  # Default options for data plot
  plot.opt <- saemix.data.setoptions(x)
  plot.opt$xlab<-paste(x@name.X," (",x@units$x,")",sep="")
  plot.opt$ylab<-paste(x@name.response," (",x@units$y,")",sep="")
  plot.opt$type<-ifelse(individual,"b","l")
  
  # Replace default options by options passed explicitly
  if(length(userPlotOptions)>0)
    plot.opt <- modifyList(plot.opt, userPlotOptions[intersect(names(userPlotOptions), names(plot.opt))])
  
  #plot.opt<-replace.data.options(plot.opt,...)
  logtyp<-paste(ifelse(plot.opt$xlog,"x",""),ifelse(plot.opt$ylog,"y",""),sep="")
  if(individual) { # separate plots subject per subject
    if(length(plot.opt$ilist)>plot.opt$nmax & plot.opt$limit) {
      if(plot.opt$interactive) {
        x1<-readline(prompt=paste("The number of subjects may be too large to be plotted. Should I plot only",plot.opt$nmax,"subjects ? (Y/n) \n"))
        if(tolower(x1)=="y") {
          plot.opt$limit<-TRUE
          plot.opt$ilist<-plot.opt$ilist[1:plot.opt$nmax]
          if(plot.opt$sample) plot.opt$ilist<-sort(sample(plot.opt$ilist, plot.opt$nmax)) else plot.opt$ilist<-plot.opt$ilist[1:plot.opt$nmax]
          if(!plot.opt$ask) {
            x1<-readline(prompt="Stop after each page of plot ? (Y/n) \n")
            if(tolower(x1)=="y") plot.opt$ask<-TRUE
          }
        }
      } else {
        if(plot.opt$interactive) {
          cat("The number of subjects is too large, I will plot only")
          if(plot.opt$sample) cat(" the data for",plot.opt$nmax,"subjects sampled randomly;") else cat(" only the data for the first",plot.opt$nmax,"subjects;")
          cat(" use limit=FALSE in the call to plot to force plotting all the subjects.\n")
        }
        if(plot.opt$sample) plot.opt$ilist<-sort(sample(plot.opt$ilist, plot.opt$nmax)) else plot.opt$ilist<-plot.opt$ilist[1:plot.opt$nmax]
      }
    } # end of test on length(ilist)
    if(plot.opt$new) {
      if(length(plot.opt$mfrow)==0) {
        np<-length(plot.opt$ilist)
        if(np>12) np<-12
        n1<-round(sqrt(np))
        n2<-ceiling(np/n1)
        par(mfrow=c(n1,n2),ask=plot.opt$ask)
      } else par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
    }
    xind<-x["data"][,x["name.predictors"], drop=FALSE]
    id<-x["data"][,"index"]
    yobs<-x["data"][,x["name.response"]]
    for(isuj in plot.opt$ilist) {
      if(plot.opt$main=="") main<-paste("Subject",isuj) else main<-plot.opt$main
      plot(xind[id==isuj,x@name.X],yobs[id==isuj],type=plot.type, xlab=plot.opt$xlab,ylab=plot.opt$ylab,col=plot.opt$col,pch=plot.opt$pch,log=logtyp, xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=main,cex=plot.opt$cex, cex.axis=plot.opt$cex.axis,cex.lab=plot.opt$cex.lab,lty=plot.opt$lty, lwd=plot.opt$lwd)
    }
  } else {	# One plot for all the data
    if(plot.opt$new) par(mfrow=c(1,1))
    if(plot.type=="p" | plot.type=="b") {
      plot(x@data[,x@name.X],x@data[,x@name.response],xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,pch=plot.opt$pch,log=logtyp,xlim=plot.opt$xlim, ylim=plot.opt$ylim,main=plot.opt$main,cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab) }
    if(plot.type=="l") {
      plot(x@data[,x@name.X],x@data[,x@name.response],xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,lty=plot.opt$lty,lwd=plot.opt$lwd,type="n", log=logtyp,xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=plot.opt$main, cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab)
    }
    if(plot.type=="l" | plot.type=="b") {
      for(isuj in unique(x@data[,x@name.group])) {
        lines(x@data[x@data[,x@name.group]==isuj,x@name.X], x@data[x@data[,x@name.group]==isuj,x@name.response],col=plot.opt$col, lty=plot.opt$lty,lwd=plot.opt$lwd)
      }
    }
  }
}


#' When applied to an SaemixSimData object, mirror plots are produced which help assess whether the simulated data has similar features when compared to the original data.
#'
#' @name plot-SaemixData
#' 
#' @param irep which replicate datasets to use in the mirror plot (defaults to -1, causing a random simulated dataset to be sampled from the nsim
#' simulated datasets)
#' 
#' @aliases plot,SaemixSimData-method plot,SaemixSimData plot,SaemixSimData,ANY-method
### #' @docType methods
#' @rdname plot-SaemixData
#' @importFrom graphics plot
#' @method plot SaemixSimData
#' @export 

# Check simulations using mirror plots
plot.SaemixSimData<-function(x,y,irep=-1,...) {
  #    oldpar <- par(no.readonly = TRUE)    # code line i
  #    on.exit(par(oldpar))            # code line i + 1
  # User-defined options
  userPlotOptions<-list(...)
  if(!is_missing(y) && is.list(y)) {
    userPlotOptions<-c(y,userPlotOptions)
  }
  i1<-match("interactive",names(userPlotOptions))
  if(!is.na(i1)) interactive<-as.logical(userPlotOptions[[i1]]) else interactive<-FALSE
  i1<-match("warnings",names(userPlotOptions))
  if(!is.na(i1)) printwarnings<-as.logical(userPlotOptions[[i1]]) else printwarnings<-FALSE
  if(dim(x@datasim)[1]==0) {
    if(interactive | printwarnings) message("No simulated data.\n")} else {  
      i1<-match("type",names(userPlotOptions))
      if(!is.na(i1)) {
        plot.type<-as.character(userPlotOptions[[i1]])
        plot.type<-plot.type[plot.type!="c"]
      } else plot.type<-"l"
      # Default options for data plot
      plot.opt <- saemix.data.setoptions(x)
      plot.opt$new<-FALSE
      plot.opt$plot.type<-"b"
      plot.opt$xlab<-paste(x@name.X," (",x@units$x,")",sep="")
      plot.opt$ylab<-paste(x@name.response," (",x@units$y,")",sep="")
      # Replace default options by options passed explicitly
      if(length(userPlotOptions)>0)
        plot.opt <- modifyList(plot.opt, userPlotOptions[intersect(names(userPlotOptions), names(plot.opt))])
      logtyp<-paste(ifelse(plot.opt$xlog,"x",""),ifelse(plot.opt$ylog,"y",""),sep="")
      
      if(irep[1]<0) irep<-sample(unique(x@nsim),1)
      for(irep1 in irep) {
        if(plot.opt$main==" ") tit<-paste("Mirror plot (replication ",irep1,")",sep="") else tit<-plot.opt$main
        tab<-data.frame(id=x@data[,x@name.group],x=x@data[,x@name.X], y=x@datasim$ysim[x@datasim$irep==irep1])
        if(plot.type=="p" | plot.type=="b") {
          plot(tab[,"x"],tab[,"y"],xlab=plot.opt$xlab, ylab=plot.opt$ylab, col=plot.opt$col,pch=plot.opt$pch,log=logtyp,xlim=plot.opt$xlim, ylim=plot.opt$ylim,main=tit,cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab) }
        if(plot.type=="l") {
          plot(tab[,"x"],tab[,"y"],type="n",xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,lty=plot.opt$lty,lwd=plot.opt$lwd, log=logtyp,xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=tit, cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab)
        }
        if(plot.type=="l" | plot.type=="b") {
          for(isuj in unique(tab[,"id"])) {
            lines(tab[tab[,"id"]==isuj,"x"],tab[tab[,"id"]==isuj,"y"])
          }
        }
        
      }
    }
}

####################################################################################
####		Creating an object of SaemixData class - User-level function	####
####################################################################################

#' Function to create a SaemixData object
#' 
#' This function creates a SaemixData object. The only mandatory argument is
#' the name of the dataset. If the dataset has a header (or named columns), the
#' program will attempt to detect which column correspond to ID, predictor(s)
#' and response. Warning messages will be printed during the object creation
#' and should be read for details.
#' 
#' This function is the user-friendly constructor for the SaemixData object
#' class. The read.saemixData is a helper function, used to read the dataset,
#' and is not intended to be called directly.
#' 
#' @name saemixData
#' 
#' @param name.data name of the dataset (can be a character string giving the name of a file on disk or of a dataset in the R session, or the name of a dataset
#' @param header whether the dataset/file contains a header. Defaults to TRUE
#' @param sep the field separator character. Defaults to any number of blank spaces ("")
#' @param na a character vector of the strings which are to be interpreted as NA values. Defaults to c(NA)
#' @param name.group name (or number) of the column containing the subject id
#' @param name.predictors name (or number) of the column(s) containing the predictors (the algorithm requires at least one predictor x)
#' @param name.response name (or number) of the column containing the response variable y modelled by predictor(s) x
#' @param name.covariates name (or number) of the column(s) containing the covariates, if present (otherwise missing)
#' @param name.genetic.covariates name (or number) of the column(s) containing the covariates, if present (otherwise missing)
#' @param name.mdv name of the column containing the indicator for missing variable
#' @param name.cens name of the column containing the indicator for censoring
#' @param name.occ name of the column containing the occasion
#' @param name.ytype name of the column containing the index of the response
#' @param name.X name of the column containing the regression variable to be used on the X axis in the plots (defaults to the first predictor)
#' @param units list with up to three elements, x, y and optionally covariates, containing the units for the X and Y variables respectively, as well as the units for the different covariates (defaults to empty)
#' @param verbose a boolean indicating whether messages should be printed out during the creation of the object
#' @param automatic a boolean indicating whether to attempt automatic name recognition when some colum names are missing or wrong (defaults to TRUE)
#' @details This function is the user-friendly constructor for the SaemixData object class. The read is a helper function, used to read the dataset, and is not intended to be called directly.
#' @return A SaemixData object (see \code{\link{saemixData}}).
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}, Audrey Lavenu, Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}}, \code{\link{saemixControl}},\code{\link{saemix}}
#' @examples
#' 
#' data(theo.saemix)
#' 
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' 
#' print(saemix.data)
#' 
#' plot(saemix.data)
#' @export 

saemixData<-function(name.data,header,sep,na,name.group,name.predictors, name.response,name.X, name.covariates=c(), name.genetic.covariates=c(), name.mdv="", name.cens="",name.occ="",name.ytype="", units=list(x="",y="",covariates=c()), verbose=TRUE, automatic=TRUE) {
# setting proper types for the SaemixData class
  if(missing(name.data)) {
    if(verbose) cat("Error in saemixData: please provide the name of the datafile or dataframe (between quotes)\n")
    return("Creation of SaemixData object failed")
  }
  if(is.data.frame(name.data)) {
    data_from_name.data <- TRUE
    dat <- name.data
    name.data<-deparse(substitute(name.data))
  } else {
    data_from_name.data <- FALSE
  }
  if(missing(header)) header<-TRUE
  if(missing(sep)) sep<-""
  if(missing(na)) na<-"NA" else {na<-as.character(na);na[is.na(na)]<-"NA"}
  if(missing(name.group)) name.group<-"" else name.group<-as.character(name.group)
  if(missing(name.predictors)) name.predictors<-"" else name.predictors<-as.character(name.predictors)
  if(missing(name.response)) name.response<-"" else  name.response<-as.character(name.response)
  if(missing(name.mdv)) name.mdv<-"" else  name.mdv<-as.character(name.mdv)
  if(missing(name.cens)) name.cens<-"" else  name.cens<-as.character(name.cens)
  if(missing(name.occ)) name.occ<-"" else  name.occ<-as.character(name.occ)
  if(missing(name.ytype)) name.ytype<-"" else  name.ytype<-as.character(name.ytype)
  if(missing(name.X)) name.X<-"" else name.X<-as.character(name.X)
  name.covariates<-c(as.character(name.covariates),as.character(name.genetic.covariates))
  x<-new(Class="SaemixData",name.data=name.data,header=header,sep=sep,na=na, name.group=name.group,name.predictors=name.predictors,name.X=name.X, name.response=name.response,name.covariates=name.covariates,units=units, name.mdv=name.mdv, name.cens=name.cens, name.occ=name.occ, name.ytype=name.ytype, verbose=verbose, automatic=automatic)
#  showall(x)
  if(data_from_name.data) {
    x1<-readSaemix(x, dat)
  } else
  x1<-readSaemix(x)
  if(is(x1,"SaemixData")) {
  	igen<-rep(FALSE,length(name.covariates))
  	igen[match(name.genetic.covariates,name.covariates)]<-TRUE
  	x1@ind.gen<-igen
  	if(verbose) cat("\n\nThe following SaemixData object was successfully created:\n\n")
  }
  if(verbose) print(x1,nlines=0)
  return(x1)
}

####################################################################################
####		Covariate transformation																								####
####################################################################################

# Apply transform to the data element of an object
# BUT: very dangerous !!! (eg mdv can be transformed to other than 0/1, NA values can be added, etc...)

#' @export 
transform.SaemixData<-function(`_data`, ...) {
	`_data`@data <- data.frame(transform(`_data`@data,...))
	`_data`
}

# > transform.data.frame
# function (`_data`, ...) 
# {
# 	e <- eval(substitute(list(...)), `_data`, parent.frame())
# 	tags <- names(e)
# 	inx <- match(tags, names(`_data`))
# 	matched <- !is.na(inx)
# 	if (any(matched)) {
# 		`_data`[inx[matched]] <- e[matched]
# 		`_data` <- data.frame(`_data`)
# 	}
# 	if (!all(matched)) 
# 		do.call("data.frame", c(list(`_data`), e[!matched]))
# 	else `_data`
# }

#' Transform covariates
#' 
#' Transform and/or center a vector
#'
#' @name transform
#' @aliases transform.numeric
#' 
#' @param _data a vector with values of type numeric
#' @param transformation transformation function. Defaults to no transformation
#' @param centering string, giving the value used to center the covariate; can be "mean" or "median", in which case this value will be computed from the data, 'none' or 0 for no centering, or a value given by the user. Defaults to the median value over the dataset.
#' @param verbose a boolean, prints messages during the execution of the function if TRUE. Defaults to FALSE.
#' @param \dots unused, for consistency with the generic method
#' @examples 
#' # TODO
#' @return a vector
#' @keywords data
#' @export

transform.numeric<-function(`_data`,transformation=function(x) x, centering="median",verbose=FALSE, ...) {
  x <- `_data`
  if(!(centering %in% c('mean','median')) & is.na(as.double(centering))) {
    if(verbose) cat("Need a proper value to center. Please specify mean, median or a numerical value\n")
    return(x)
  }
  if(tolower(centering)=="none") centering<-0
  if(centering %in% c('mean','median')) {
    f1<-match.fun(centering)
    xcent<-f1(x)
  } else xcent<-as.double(centering)
  if(verbose) cat("Data centered with respect to the value:",xcent,"\n")
  xt<-transformation(x-xcent)
  return(xt)
}


#' Transform covariates
#' 
#' Transform and/or center continuous covariates
#'
#' @name transformContCov
#' @aliases transform.SaemixData
#' 
#' @param object saemixData object
#' @param covariate name of the covariate
#' @param transformation transformation function. Defaults to no transformation
#' @param centering string, giving the value used to center the covariate; can be "mean" or "median", in which case this value will be computed from the data, 'none' or 0 for no centering, or a value given by the user. Defaults to the median value over the dataset.
#' @param verbose a boolean, prints messages during the execution of the function if TRUE. Defaults to FALSE.
#' @examples 
#' # TODO
#' @return an object of class \code{"\linkS4class{SaemixData}"}
#' @keywords data
#' @export
#' 
transformContCov<-function(object, covariate, transformation=function(x) x, centering="median",verbose=FALSE) {
 	covariate<-deparse(substitute(covariate))
 	name.trans<-deparse(substitute(transformation))
	if(!(covariate %in% object@name.covariates)) {
		if(verbose) cat("Covariate",covariate,"not found\n")
		return(object)
	}
	if(!(centering %in% c('mean','median')) & is.na(as.double(centering))) {
	  if(verbose) cat("Need a proper value to center. Please specify mean, median or a numerical value\n")
		return(object)
	}
	if(tolower(centering)=="none") centering<-0
	if(centering %in% c('mean','median')) {
		f1<-match.fun(centering)
		labl<-paste(object@data[,object@name.group],object@data$occ,sep="-")
		covar<-object@data[!duplicated(labl),covariate]
		xcent<-f1(covar)
	} else xcent<-as.double(centering)
	if(verbose) cat(covariate,"centered with respect to the value:",xcent,"\n")
	xcov<-object@data[,covariate]
	xcov<-transformation(xcov-xcent)
	object@data[,covariate]<-xcov
	object@trans.cov[[covariate]]<-list(type="continuous",transformation=transformation,centering=centering)
	return(object)
}

#' Transform covariates
#' 
#' Regroup categorical covariates
#'
#' @name transformCatCov
#' 
#' @param object saemixData object
#' @param covariate name of the covariate
#' @param group a vector giving the categories to which the initial values of the covariates should be mapped. If the resulting covariate is binary, it will be stored as 0/1. If it has more than 2 categories, dummy covariates will be created for the analysis.
#' @param reference the reference group
#' @param verbose a boolean, prints messages during the execution of the function if TRUE. Defaults to FALSE.
#' @return an object of class \code{"\linkS4class{SaemixData}"}
#' @examples 
#' # TODO
#' @keywords data
#' @export 
 
transformCatCov<-function(object, covariate, group, reference, verbose=FALSE) {
	covariate<-deparse(substitute(covariate))
	if(!(covariate %in% object@name.covariates)) {
	  if(verbose) cat("Covariate",covariate,"not found\n")
		return(object)
	}
	if(length(object@ocov)>0) xcov<-object@ocov[,covariate] else xcov<-object@data[,covariate]
	if(missing(reference)) reference<-sort(group)[1]
	if(!(reference %in% group) & verbose) {
	  if(verbose) cat("Reference category not in group\n")
		reference<-sort(group)[1]
	}
	if(length(group)>10 & verbose) {
	  if(verbose) cat("Warning: more than 10 categories\n")
	}
	ifac<-(is.factor(xcov))
	gr<-as.character(xcov)
	ugr<-sort(unique(gr))
	if(length(ugr)!=length(group)) {
	  if(verbose) cat("The argument group must be the same size as the initial number of categories\n")
		return(object)
	}
	for(i in 1:length(ugr)) gr[as.character(xcov)==ugr[i]]<-group[i]
	if(ifac) gr<-as.factor(gr)
	uugr<-unique(gr)
	uugr<-c(reference,uugr[uugr!=reference])
	ngr<-length(uugr)
	if(ngr>2) { # remove initial covariate from data object
		object@data[,covariate]<-NULL
		tdum<-NULL # generate dummy covariates
		for(i in 2:ngr) {
			dum<-ifelse(gr==uugr[i],1,0)
			tdum<-cbind(tdum,dum)
			colnames(tdum)[(i-1)]<-paste(covariate,".G",i,sep="")
		}
		object@data<-cbind(object@data,tdum)
	} else { # 2 categories, remapping to 0/1
		xgr<-ifelse(unclass(gr)==reference,0,1)
		object@data[,covariate]<-xgr
	}
	object@trans.cov[[covariate]]<-list(type="cat",group=group,reference=reference)
	return(object)
}

####################################################################################
####				saemixObject class - S3 methods			####
####################################################################################

#' Data subsetting
#' 
#' Return an SaemixData object containing the subset of data which meets conditions.
#'
#' @name subset
#' @aliases subset-methods subset.SaemixData
#' 
#' @param x saemixData object
#' @param subset logical expression indicating elements or rows to keep: missing values are taken as false
#' @param ... additional parameters (ignored)
#' @return an object of class \code{"\linkS4class{SaemixData}"}
#' @examples 
#' # TODO
#' @keywords methods
#' @export 

subset.SaemixData<-function (x, subset, ...) {
    if (missing(subset)) 
        return(x)
    else {
        e <- substitute(subset)
	      xdat<-x["data"]
        r <- eval(e, xdat, parent.frame())
        if (!is.logical(r)) 
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
    }
    x1<-x
    x1["data"]<-x["data"][r,,drop=FALSE]
    if(length(x["yorig"])>0) x1["yorig"]<-x["yorig"][r]
    if(length(x["ocov"])>0) x1["ocov"]<-x["ocov"][r,,drop=FALSE]
    id<-x1["data"][,x1["name.group"]]
    x1["N"]<-length(unique(id))
    nind.obs<-tapply(id,id,length) # individual numbers of observations (1xN)
    nind.obs<-c(nind.obs[match(unique(id),names(nind.obs))])
    x1["nind.obs"]<-nind.obs
    x1["ntot.obs"]<-length(id)
    x1["data"]$index<-rep(1:x1["N"],times=nind.obs)
    return(x1)
}

####################################################################################
