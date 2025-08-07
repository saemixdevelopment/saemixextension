
#' @include aaa_generics.R
#' @include SaemixData.R
#' @include SaemixData-methods.R
NULL

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


#' Transform continuous covariates
#' 
#' Transform and/or center continuous covariates
#'
#' @name transformContCov
#' @aliases transform.SaemixData
#' 
#' @param object saemixData object
#' @param covariate the name of the covariate to transform (without quotes)
#' @param transformation function used to apply the transformation (can be a name, eg log or exp, or a function, see examples). Defaults to no transformation
#' @param centering string, giving the value used to center the covariate; can be "mean" or "median", in which case this value will be computed from the data, 'none' or 0 for no centering, or a value given by the user. Defaults to the median value over the entire dataset (including potential duplicates). The transformation will be applied as transformation(covariate)-transformation(centering value)
#' @param na.rm a logical evaluating to TRUE or FALSE indicating whether NA values should be stripped before computing mean or median for centering (defaults to TRUE)
#' @param newCovName the name of the transformed covariate (if not given, will be set to "the name of the covariate to transform"+".mod")
#' @param verbose a boolean, prints messages during the execution of the function if TRUE. Defaults to FALSE.
#' @details
#' Transformations can only involve one covariate. More complex transformations should be performed by the user 
#' and integrated in the data object before a call to saemixData
#' NA values are kept unchanged (values NA before will be NA after transformation)
#' Using transformContCov on an saemixData object will create or modify the trans.cov slot, to store the transformation leading
#' to the new covariate.
#' Warning: the name of the covariate must be given without quotes, to match the behaviour of the usual transform function in R.
#' 
#' @examples 
#' # Log-transform variable birthyear and center it to the median
#'   x<-saemixData(name.data=cow.saemix,name.group="cow",name.predictors=c("time"), 
#'       name.response="weight",
#'   name.covariates=c("birthyear","twin","birthrank"), 
#'       units=list(x="d",y="kg",covariates=c("yr","-","-")),verbose=FALSE)
#'   x2<-transformContCov(x,birthyear,centering="median",transformation=log,verbose=FALSE, 
#'       newCovName = "logYear")
#'   print(summary(x2@data$logYear))
#' # Should be the same as:
#'   print(summary(x@data$birthyear-median(x@data$birthyear)))
#'   
#' @return an object of class \code{"\linkS4class{SaemixData}"}
#' @keywords data
#' @export
#' 
transformContCov<-function(object, covariate, transformation=function(x) x, centering="median", na.rm=TRUE, newCovName="",verbose=FALSE) {
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
  if(newCovName=="") newCovName <- paste0(covariate,".mod")
  if(newCovName==covariate) {
    newCovName <- paste0(covariate,".mod")
    if(verbose) cat("Transformed covariates should have a different name from the original covariate, renaming to:", newCovName,"\n")
  }
  if(tolower(centering)=="none") centering<-0
  if(centering %in% c('mean','median')) {
    f1<-match.fun(centering)
#    labl<-paste(object@data[,object@name.group],object@data$occ,sep="-")
#    covar<-object@data[!duplicated(labl),covariate]
    covar<-object@data[,covariate]
    xcent<-f1(covar, na.rm=na.rm)
  } else xcent<-as.double(centering)
  if(verbose) cat(covariate,"centered with respect to the value:",xcent,"\n")
  xcov<-object@data[,covariate]
  xcov<-transformation(xcov)-transformation(xcent)
  # Add the new covariate to the data and its name to the list of covariates
  object@data[,newCovName]<-xcov
  object@name.covariates <- c(object@name.covariates, newCovName)
  # Keep the transformation from the old to the new
  object@trans.cov[[newCovName]]<-list(fromCov=covariate,toCov=newCovName,type="continuous",transformation=transformation,centering=centering, na.rm=na.rm)
  return(object)
}

## TODO
# transformContCov2Cat<-function(object, covariate, transformation=function(x) x, centering="median", newCovname="",verbose=FALSE) {
  # Eco TODO: transform continuous covariate in categories
  ## by default: <25, 25-75, >75 with categories q25 (<25th percentile) and q75 (>75th percentile), with 25-75 as the reference
#}

## TODO: do we want a cat to cont conversion ? (not necessarily)

#' Transform categorical covariates
#' 
#' This function is used to automatically map categorical covariates to dummy binary covariates. It can also be used to regroup categorical covariates before the analysis, for example mapping a covariate initially in 5 categories to only 3 categories.
#' For binary covariates, the covariate will simply be transformed to 0/1 with the reference category being 0.
#'
#' @name transformCatCov
#' 
#' @param object saemixData object
#' @param covariate name of the covariate
#' @param oldCat a vector giving the initial categories of the covariate (if NULL, defaults to the unique values of the covariate) 
#' @param newCat a vector of consecutive integers 1:N, giving the N categories to which the initial values of the covariate should be mapped, with 1 denoting the reference class. If the resulting covariate is binary, it will be stored as 0/1. If it has more than 2 categories, N dummy covariates 0/1 will be created. If not given, the initial categories of the covariate will be mapped to 1:N and transformed into N dummy covariates 0/1.
#' @param newCatName the name of the new group. If NULL (default), binary covariates will be renamed as covariate.mod, while for categorical covariates, the new covariates will be named covariate.ref, covariate.G2, covariate.G3, etc...
#' @param verbose a boolean, prints messages during the execution of the function if TRUE. Defaults to FALSE.
#' @return an object of class \code{"\linkS4class{SaemixData}"}
#' 
#' @details
#' For a binary covariate, a dummy covariate with values 0 for the reference class and 1 for the other category will be created and will replace the original dataset in the covariate
#' For covariates with 3 categories or more (categorical covariates), dummy variables with values 0/1 will be created for the reference class and for each contrast with the reference class (category 2 versus reference, category 3 versus reference, etc...)
#' If these covariates have units, the dummy covariates will have the same unit (which may not be appropriate and can then be changed by the user)
#' 
#' @examples 
#' data(cow.saemix)
#' saemix.data<-saemixData(name.data=cow.saemix,header=TRUE,name.group=c("cow"),
#'                    name.predictors=c("time"),name.response=c("weight"),
#'                    name.covariates=c("birthyear","twin","birthrank"),
#'                    units=list(x="days",y="kg",covariates=c("yr","-","-")))
#' unique(saemix.data@data$birthrank) # 5 categories, 3 4 5 6 7
#' # create 3 dummy variables regrouping 3 (reference), 4 and 5, and 6 and 7
#' cowt <- transformCatCov(saemix.data, covariate=birthrank, newCat=c(1,2,2,3,3), verbose=TRUE) 
#' head(saemix.data@data) # the original covariate is birthrank
#' head(cowt@data) 
#' # the new covariates are birthrank.ref (initially 3), birthrank.G2 (regrouping 4 and 5) and 
#' # birthrank.G3 (6 and 7)
#' # only birthrank.G2 and birthrank.G3 are included in the name.covariates slot of the object cowt
#' cowt <- transformCatCov(cowt, covariate=birthrank, newCat=c(1,2,2,3,3), 
#'    newCatName=c("ref","preg4-5","6-7"), verbose=TRUE)
#' head(cowt@data) 
#' # new names can be assigned to the dichotomised covariates through the newCatName argument
#' # the new covariates are now "ref", preg4-5" (regrouping 4 and 5) and "6-7" (6 and 7)
#' 
#' # Changing the reference for a binary variable
#' cowt<-transformCatCov(cowt, covariate=twin, newCat=c(2,1), newCatName="Singleton", verbose=TRUE)
#' 
#' @keywords data
#' @export 
#' 

transformCatCov<-function(object, covariate, oldCat=NULL, newCat=NULL, newCatName=NULL, verbose=FALSE) {
  covariate<-deparse(substitute(covariate))
  if(!(covariate %in% object@name.covariates)) {
    if(verbose) cat("Covariate",covariate,"not found\n")
    return(object)
  }
  if(sum(is.na(as.integer(newCat)))>0) {
    if(verbose) cat("newCat should be a vector of consecutive integers\n")
    return(object)
  }
  xcov<-object@data[,covariate]
  # If the covariate is binary, check if it is also in ocov and then use the value in ocov (original coding in the database)
  if(length(unique(xcov))==2) {
    if(length(object@ocov)>0 && covariate %in% colnames(object@ocov)) xcov <- object@ocov[,covariate]
  }
  # If old categories not given, use the default order for the initial categories else check the categories match
  if(is.null(oldCat)) oldCat <- sort(unique(xcov)) else {
    if(length(unique(xcov))!=length(oldCat) | !identical(sort(unique(xcov)), sort(unique(oldCat)))) {
      if(verbose) cat("The categories given in oldCat (N=",length(oldCat),") should match the categories in the original data (N=",length(unique(xcov)),")\n")
      return(object)
    }
  }
  # If new categories not given, match old categories to 1:nb.cat
  if(is.null(newCat)) {
    newCat<-1:length(oldCat)
  }
  # Check the number of new categories match the number of old ones
  if(length(newCat)!=length(oldCat)) {
    if(verbose) cat("The number of new categories should match the number of categories in the original data (",length(oldCat),")\n")
    return(object)
  }
  if(length(unique(newCat))>10 & verbose) {
    if(verbose) cat("Warning: more than 10 categories\n")
  }
  nb.cat <- length(unique(newCat))
  # Checking names
  if(is.null(newCatName)) {
    if(nb.cat==2) newCatName <- paste0(covariate,".mod") else newCatName <- covariate
  }
  if(length(newCatName)==1 & nb.cat>2) newCatName <- c(paste0(newCatName,".ref"),paste0(newCatName,".G",2:nb.cat))
  if(nb.cat==2 & length(newCatName)>1) {
    if(verbose) cat("Warning: for a binary covariate, newCatName should be of length 1, using the first element only \n")
    newCatName<-newCatName[1]
  }
  if(nb.cat>2 & length(newCatName)<nb.cat) {
    if(verbose) cat("Warning: number of newCatNames should be equal to the number of groups, changing name to",covariate,".ref,",covariate,".G2, .G3, ...\n")
    newCatName <- c(paste0(covariate,".ref"),paste0(covariate,".G",2:nb.cat))
  }
  if(length(newCatName) > nb.cat) {
    if(verbose) cat("Warning: number of newCatNames should be equal to the number of groups, truncating\n")
    newCatName<-newCatName[1:nb.cat]
  }
  oldCat <- as.character(oldCat)
  xcov <- as.character(xcov)
  uniqCat <- unique(newCat)
  newCat <- as.character(newCat)
  reference <- which(uniqCat==1)
  # Matching old categories to new ones
  xgroup<-newCat[match(xcov,oldCat)]
  # Removing old covariate from dataframe
  object@data[,covariate]<-NULL
  if(length(object@units$covariates)>0) {
    unitcov<-object@units$covariates[which(object@name.covariates==covariate)]
    object@units$covariates<-object@units$covariates[-which(object@name.covariates==covariate)]
  }
  object@name.covariates <-object@name.covariates[object@name.covariates!=covariate]
  # Creating binary variables 0/1 for the fit
  if(nb.cat>2) { # remove initial covariate from data object
    tdum<-NULL # generate dummy covariates
    for(i in 1:nb.cat) {
      dum<-ifelse(xgroup==uniqCat[i],1,0)
      tdum<-cbind(tdum,dum)
      # Add covariate name to name.covariates (except for reference class)
      if(i!=reference) {
        object@name.covariates <-c(object@name.covariates,newCatName[i])
        if(length(object@units$covariates)>0) object@units$covariates<-c(object@units$covariates, unitcov)
      }
      # Keep the transformation from the old to the new 
      object@trans.cov[[newCatName[i]]]<-list(fromCov=covariate,toCov=newCatName[i],type="cat", fromValue=oldCat[newCat==uniqCat[i]], toGroup=uniqCat[i])
    }
    colnames(tdum)<-newCatName
    object@data<-cbind(object@data,tdum)
  } else { # 2 categories, remapping to 0/1
    xgr<-ifelse(xgroup=="1",0,1)
    object@data[,newCatName]<-xgr
    object@trans.cov[[newCatName]]<-list(fromCov=covariate,toCov=newCatName,type="cat", fromValue=oldCat[newCat!=1], toGroup=1)
  }
  return(object)
}

# Eco TODO: add handling of NA
## maybe add impute.NA argument with options 'no'/none (=> stays NA), "mode" (most frequent value), "XX" (a specific value)
# Eco TODO: check if ocov not having the same names is a problem => removed ocov bits, now only shows the covariates in the model itself

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

