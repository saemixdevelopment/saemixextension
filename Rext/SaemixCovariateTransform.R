########################################################################
#' Transforming covariates
#' 
#' This method allows to apply the transformation contained in a SaemixContinuousCovariate or SaemixDiscreteCovariate
#' object to a vector of covariate values in order to include the covariate in a saemix model.
#' 
#' @name transformCovariate-methods
#' @aliases transformCovariate 
#' @aliases transformCovariate,SaemixContinuousCovariate transformCovariate,SaemixContinuousCovariate-method
#' @aliases transformCovariate,SaemixDiscreteCovariate transformCovariate,SaemixDiscreteCovariate-method
#' 
#' @param object an object of class SaemixDiscreteCovariate or SaemixContinuousCovariate
#' @param x the values of the covariate to transform
#'     
#' @return a vector or dataframe with the transformed values
#' 
#' @details 
#' For continuous covariates, the transformation contained in object is applied to the values in x.
#' For binary covariates, x is transformed to 0 for the values corresponding to the reference value and 1 otherwise.
#' For categorical covariates with ncat (=3 or more) categories, a dataframe is created with (ncat-1) dummy variables in columns. Each column is 1 if the corresponding value of x is in the corresponding group and 0 otherwise.
#' 
#' @examples 
#' # Transforming a vector of weight to log(weight/mean(weight))
#' weightCov<-new(Class="SaemixContinuousCovariate", name="Weight", transform.function=log, centering.function=median)
#' print(c(60, 70, 80))
#' transformCovariate(weightCov, c(60, 70, 80))
#' 
#' # Transforming a gender covariate given as Female/Male to 0/1 with Female as reference (0)
#' sexCov<-new(Class="SaemixDiscreteCovariate", name="gender", reference="Female")
#' transformCovariate(sexCov, c("Female","Male","Male"))
#' # Also works with factors
#' transformCovariate(sexCov, as.factor(c("Female","Male","Male")))
#' 
#' # Regrouping a covariate with 5 categories in 3 categories
#' scoreCov<-new(Class="SaemixDiscreteCovariate", name="score", type="categorical",groups=list(c(1,2),c(3,4),c(5)))
#' transformCovariate(scoreCov, c(1,2,3,4,5))
#' 
#' @exportMethod transformCovariate

# Generic, move to aaa_generics.R
setGeneric(name="transformCovariate",
           def=function(covariateModel, x=NULL) standardGeneric("transformCovariate")
)


setMethod(
  f ="transformCovariate",
  signature = "SaemixContinuousCovariate" ,
  definition = function (covariateModel, x=NULL){
    if(is.null(x)) return(NULL)
    if(length(covariateModel@centering.value)==0) 
      covariateModel@transform.function(x/covariateModel@centering.function(x)) else 
        covariateModel@transform.function(x/covariateModel@centering.value)
  }
)

setMethod(
  f ="transformCovariate",
  signature = "SaemixDiscreteCovariate" ,
  definition = function (covariateModel, x=NULL){
    if(is.null(x)) return(NULL)
    if(covariateModel@reference=="") reference<-sort(unique(as.character(x)))[1] else reference<-as.character(covariateModel@reference)
    if(covariateModel@type=="binary" & length(unique(x))<=2) {
      return(1-as.integer(as.character(x)==reference))
    } else { # type="categorical" OR x has more than 2 unique values
      ncat<-length(unique(x))
      if(ncat>length(x)/2) {
        message("x may not be a categorical covariate, it seems to have too many different values")
      }
      # Create (ncat-1) categories
      if(is.factor(x)) xcat<-levels(x) else xcat<-sort(unique(as.character(x)))
      if(length(covariateModel@groups)==0) {
        groups<-as.vector(xcat, mode="list")
        names(groups)<-xcat
      } else groups<-covariateModel@groups
      if(is.null(names(groups))) {
        l1<-c()
        for(i in 1:length(groups)) l1<-c(l1,as.character(groups[[i]][1]))
        names(groups)<-l1
      }
      idx1<-grep(reference,groups)
      if(length(idx1)==0) {
        message(paste("Can't find reference category ", reference,", using first category instead"))
        reference<-groups[[1]][1]
        idx1<-1
      }
      ncat<-length(groups)
      tabcat<-NULL
      for(i in 1:ncat) {
        if(i!=idx1) {
          icat<-as.integer(as.character(x) %in% groups[[i]])
          tabcat<-cbind(tabcat, icat)
        }
      }
      colnames(tabcat)<-paste0(covariateModel@name,".",names(groups)[-idx1])
      if(dim(tabcat)[2]==1) tabcat<-c(tabcat) # regrouped to binary covariate
      return(tabcat)
    }
  }
)
