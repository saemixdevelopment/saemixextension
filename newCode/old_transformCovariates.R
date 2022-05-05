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
