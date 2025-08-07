#' @include aaa_generics.R
#' @include SaemixParameter.R
NULL

################################################################################
# Computational functions
## TODO remove from func_aux.R

#' @rdname saemix.internal
#' 
#' @aliases normcdf norminv
#' 
#' @keywords internal

# Inverse of the normal cumulative distribution fct: using erfcinv from ?pnorm
norminv<-function(x,mu=0,sigma=1)  mu-sigma*qnorm(x,lower.tail=FALSE)

# Truncated gaussian distribution (verifie par rapport a definition de erf/matlab)
normcdf<-function(x,mu=0,sigma=1)
  cutoff(pnorm(-(x-mu)/sigma,lower.tail=FALSE),1e-30)

################################################################################
# Creator function 

#' Function to create an SaemixParam object
#' 
#' This function creates a parameter object (class: SaemixParam). 
#' This function is the user-friendly constructor for the SaemixParam object class, and is the preferred way to specify parameters starting from saemix 4.0. 
#' 
#' @name saemixParam
#' 
#' @param name name of the parameter
#' @param verbose a boolean indicating whether messages should be printed out during the creation of the object
#' 
#' @details Starting from saemix 4.0, parameters are specified through their distribution.
#' The default value, when called without any argument, is to return a parameter "theta" with a lognormal distribution, one level of variability and initial values of 1 for both the parameter and its standard deviations.
#' This new definition allows for setting initial values, additional levels of variability, correlation structures (when combined with other parameters), as well as specifying which of the associated population parameters are fixed or estimated.
#' 
#' TODO: list of parameters
#' 
#' 
#' @return An SaemixParam object (see \code{\link{SaemixParam}}).
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' DDMORE
#' 
#' MLXTRAN
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}}, \code{\link{saemixControl}},\code{\link{saemix}}
#' @examples
#' # Defining a single parameter, with default values
#' ka<-saemixParam()
#' print(ka)
#'   
#' # Defining the list of parameters for a model
#' saemix.parameters <- list(ka=saemixParam(sd.init=c(0.8,0.5), varlevel=c("iiv","iov")),
#'      vd=saemixParam(sd.init=0.7),  
#'      cl=saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(0.6,0.3), corr = list(iiv=c("ka","vd"),iov=c("vd")), covariate=c("wt","sex","age"), covariate.init=c(0.75,0,0), covariate.estim=c("fixed","estimated","estimated"), corr.init=list(iiv=c(-0.5,0.7), iov=0.7)))
#' print(saemix.parameters)
#' @export 


## ToDo: secure code against type mismatches

saemixParam <- function(name,distribution="lognormal", transform="",invtransform="", mu.estim='estimated', mu.init=1, varlevel="iiv", sd.init=c(), sd.estim=character(), corr=list(), corr.init=list(),  corr.estim=list(), covariate=c(), covariate.estim=c(),covariate.init=c(), covariate.varlevel=c(), verbose=FALSE) {
  if(missing(name)) name<-"theta"
  theta<-new(Class="SaemixParameter", name=name, distribution=distribution, transform=transform, invtransform=invtransform, 
             mu.init=mu.init, varlevel=varlevel, covariate=covariate, 
             covariate.estim=covariate.estim, covariate.init=covariate.init, covariate.varlevel=covariate.varlevel, verbose=verbose)
  # add sd.init and sd.estim, matching the varlevels
  names(varlevel)<-varlevel
  nvar<-length(varlevel)
  if(nvar>0) {
    if(length(sd.init)==0) sd.init<-rep(1,length(varlevel))
    if(length(sd.init)<nvar) sd.init<-c(sd.init, rep(1,nvar-length(sd.init)))
    if(length(sd.estim)==0) sd.estim<-rep("estimated",length(varlevel))
    if(length(sd.estim)<nvar) sd.estim<-c(sd.init, rep("estimated", nvar-length(sd.init)))
  } else {
    sd.init<-numeric()
    sd.estim<-character()
  }
  theta@sd.init<-sd.init
  theta@sd.estim<-sd.estim
  logmsg<-""
  
  # add correlation and CI for corr to the object
  # Option 1: replicate corr structure => no
  # Option 2: use a named list for corr to match the levels in varlevel
  if(length(corr)>0 & length(varlevel)==0) {
    msg<-"corr given but no associated variability, ignoring\n"
    if(verbose) cat(msg)
    logmsg<-paste0(logmsg, msg)
    corr<-list()
  }
  ## match names and structure in corr, corr.init, corr.estim
  if(length(corr)>0) {
    xcheck<-checkMatchingVectorList(varlevel, corr, type="character", defaultValue = c())
    logmsg<-paste0(logmsg, xcheck$logmsg)
    corr<-xcheck$list
    
    xcheck <- checkMatchingList(corr, corr.init, name.list1="corr",name.list2="corr.init")
    corr.init <- xcheck$list
    ilog<-0
    for(ilev in 1:length(corr.init)) {
      xinit<-corr.init[[ilev]]
      idx1<-sum(xinit>1)+sum(xinit<(-1))
      if(idx1>0) {
        xinit[xinit>1]<-0.5 
        xinit[xinit<(-1)]<-(-0.5)
        if(ilog==0) {
          msg<-"Replacing CI for correlation not in [-1,1] by +/-0.5\n"
          if(verbose) cat(msg)
          theta@log<-c(theta@log,msg)
          ilog<-1
        }
      }
    }
    logmsg<-paste0(logmsg, xcheck$logmsg)
    
    xcheck <- checkMatchingList(corr, corr.estim, name.list1="corr",name.list2="corr.estim", type="character", defaultValue='estimated')
    corr.estim <- xcheck$list
    logmsg<-paste0(logmsg, xcheck$logmsg)
  } else {
    corr.init<-list()
    corr.estim<-list()
  }
  theta@corr<-corr
  theta@corr.init <- corr.init
  theta@corr.estim <- corr.estim
  return(theta)
}

################################################################################

## Helper function, checks if the names given by the user match to the names in the dataset. If not, automatic recognition is attempted when automatic=TRUE.

#' @rdname internals
#' @alias checkParameters
#' 
#' @param parameters a list of parameters
#' @param varlevel a list of variability levels. If not given, this function will attempt to derive the variability levels and their nesting order through the information in parameters
#' 
#' @return if parameters is not a list, returns an error message, otherwise returns a list with three elements: 
#' parameters: a list of named parameters with valid names
#' varlevel.order: a list of variability levels, sorted by nesting order
#' logmsg: a string recording the warnings during the check
#' 
#' @examples 
#' # TODO
#' @keywords methods
#' @export 
NULL


# Check parameters, assign names and return the order of variability levels
checkParameters <- function(parameters, varlevel=c()) {
  if(inherits(parameters,"SaemixParameter")) {
    parameters<-list(parameters)
    name(parameters)<-parameters[[1]]@name
  } else {
    icheck<-0
    for(i in parameters)
      if(!inherits(i,"SaemixParameter")) icheck<-1
    if(icheck==1) {
      return("Extraction failed: parameters should be a list of SaemixParameter objects\n")
    }
  }
  npar<-length(parameters)
  logmsg<-""
  if(is.null(names(parameters))) {
    names(parameters)<-paste0("theta",1:npar)
    for(i in 1:npar)
      if(parameters[[i]]@name=="theta") parameters[[i]]@name<-names(parameters)[i] else names(parameters)[i]<-parameters[[i]]@name
  } else {
    for(i in 1:length(parameters)) 
      if(names(parameters)[i]!="") parameters[[i]]@name<-names(parameters)[i] else {
        if(!is.null(parameters[[i]]@name)) names(parameters)[i]<-parameters[[i]]@name
      }
  }
  nampar <- names(parameters)
  if(length(varlevel)==0) {
    varlevel.order<-c()
    for(i in 1:npar) {
      if(length(parameters[[i]]@varlevel)>length(varlevel.order)) varlevel.order<-parameters[[i]]@varlevel
      varlevel<-c(varlevel,parameters[[i]]@varlevel)
    }
    varlevel<-unique(varlevel)
    if(length(varlevel)>length(varlevel.order)) {
      msg<-paste0("None of the parameters have all the levels of variability, assuming the order is:",varlevel,"\n    If this is not correct, please specify the order of the variability levels using the varlevel argument\n")
      varlevel.order<-varlevel
      logmsg<-paste0(logmsg,msg)
      if(verbose) cat(msg)
    }
  } else {
    varlevel.order<-varlevel
    # add checks
    ## ignore additional levels of variability in the parameters
    ## if varlevel.order is different for some parameters, print warning message
  }
  return(list(parameters=parameters,varlevel.order=varlevel.order,logmsg=logmsg))
}

################################################################################
# Utility functions to automatically fill in omega.estim and omega.init 

#' @rdname validate.names
#' @aliases checkMatchingList
#' 
#' @param list1 a list to match against
#' @param list2 a list to match
#' @param name.list1 name of list1 (for informative error messages)
#' @param name.list2 name of list2
#' @param checkContent a boolean. If TRUE (default), the structure of list2 will be 
#' @param defaultValue vector of strings, values for automatic recognition
#' @param type type of the elements in list2 
#' @param verbose whether to print messages
#' 
#' @details
#' checkMatchingList verifies and changes the content of list2 to match the names, type and structure of list1.
#' checkMatchingListVector verifies and changes the content of list2 to match the names and type of vec1.
NULL


checkMatchingList <- function(list1, list2, checkContent=TRUE, name.list1="",name.list2="", defaultValue=0, type="numeric", verbose=FALSE) {
  # match names in list2 to names in list1, potentially reordering
  # make sure both lists have the same subelement
  # check type is type
  logmsg<-""
  if(length(list2)>0) { 
    if(is.null(names(list2))) { # names not given, assume same order and levels as list1
      if(length(list2)>length(list1)) {
        msg<-paste0("More elements in ",name.list2," than in ",name.list1," ignoring additional levels\n")
        if(verbose) cat(msg)
        logmsg<-paste0(logmsg, msg)
        list2<-list2[1:length(list1)]
      }
      names(list2)<-names(list1)
    } else {# match names and potentially re-order
      temp <- list()
      for(ilev in names(list1)) {
        idx<-grep(ilev, names(list2))
        if(length(idx)>0) temp[ilev]<-list2[idx] else temp[ilev]<-rep(defaultValue, length(list1[ilev]))
      }
      if(!identical(temp,list2)) {
        msg<-paste0("Matching elements in ",name.list2," to ",name.list1,"\n")
        if(verbose) cat(msg)
        logmsg<-paste0(logmsg, msg)
      }
      list2<-temp
    }
    # check type and length of each element
    if(checkContent) {
      ilog.type<-ilog.ign<-ilog.comp<-0
      for(ilev in 1:length(list2)) {
        idx<-grep(names(list2)[ilev], names(list1))
        nelem <- length(list1[[idx]])
        if(!is(list2[[ilev]],type)) {
          if(ilog.type==0) {
            msg<-paste0("Elements in ",name.list2," should be of type ",type,", setting to default value\n")
            if(verbose) cat(msg)
            logmsg<-paste0(logmsg, msg)
            ilog.type<-1
          }
          list2[[ilev]]<-rep(defaultValue, nelem)
        }
        if(length(list2[[ilev]])>nelem) {
          if(ilog.ign==0) {
            msg<-paste0("Too many elements in ",name.list2,", ignoring extra elements\n")
            if(verbose) cat(msg)
            logmsg<-paste0(logmsg, msg)
            ilog.ign<-1
          }
          list2[[ilev]]<-list2[[ilev]][1:nelem]
        }
        if(length(list2[[ilev]])<nelem) {
          if(ilog.comp==0) {
            msg<-paste0("Too few elements in ",name.list2,", adding default values\n")
            if(verbose) cat(msg)
            logmsg<-paste0(logmsg, msg)
            ilog.comp<-1
          }
          list2[[ilev]]<-c(list2[[ilev]], rep(defaultValue,nelem))
          list2[[ilev]]<-list2[[ilev]][1:nelem]
        }
      }
    }
  } else {
    list2<-vector(mode="list", length=length(list1))
    names(list2)<-names(list1)
    if(checkContent) {
      for(ilev in 1:length(list1)) 
        list2[[ilev]]<-rep(defaultValue,length(list1[[ilev]]))
    }
  }
  return(list(list=list2, logmsg=logmsg))
}

#' @rdname validate.names
#' @aliases checkMatchingVectorList
#' 
#' @param vec1 a vector to match names against
#' @param name.vec1 name of vec1 (for informative error messages)
NULL

checkMatchingVectorList <- function(vec1, list2, checkContent=TRUE, name.vec1="",name.list2="", defaultValue=0, type="numeric", verbose=FALSE) {
  # match names in list2 to names in named vector vec1, potentially reordering
  # check type of list2 is type, attribute default value
  logmsg<-""
  if(is.null(names(vec1))) names(vec1)<-vec1
  if(length(list2)>0) { 
    if(is.null(names(list2))) { # names not given, assume same order and levels as vec1
      if(length(list2)>length(vec1)) {
        msg<-paste0("More elements in ",name.list2," than in ",name.vec1," ignoring additional levels\n")
        if(verbose) cat(msg)
        logmsg<-paste0(logmsg, msg)
        list2<-list2[1:length(vec1)]
      }
      if(length(list2)<length(vec1)) {
        msg<-paste0("Less elements in ",name.list2," than in ",name.vec1,", using only the first ",length(list2)," elements of ",name.vec1,"\n")
        if(verbose) cat(msg)
        logmsg<-paste0(logmsg, msg)
        vec1<-vec1[1:length(list2)]
      } 
      names(list2)<-names(vec1)
    } else {# match names and potentially re-order
      temp <- list()
      for(ilev in names(vec1)) {
        idx<-grep(ilev, names(list2))
        if(length(idx)>0) temp[ilev]<-list2[idx] else temp[ilev]<-defaultValue
      }
      if(!identical(temp,list2)) {
        msg<-paste0("Matching elements in ",name.list2," to ",name.vec1,"\n")
        if(verbose) cat(msg)
        logmsg<-paste0(logmsg, msg)
      }
      list2<-temp
    }
    # check type and length of each element
    if(checkContent) {
      ilog.type<-ilog.ign<-ilog.comp<-0
      for(ilev in 1:length(list2)) {
        idx<-grep(names(list2)[ilev], names(vec1))
        if(!is(list2[[ilev]],type)) {
          if(ilog.type==0) {
            msg<-paste0("Elements in ",name.list2," should be of type ",type,", setting to default value\n")
            if(verbose) cat(msg)
            logmsg<-paste0(logmsg, msg)
            ilog.type<-1
          }
          list2[[ilev]]<-rep(defaultValue, length(list2[[ilev]]))
        }
      }
    }
  } else {
    list2<-list()
  }
  return(list(list=list2, logmsg=logmsg))
}
