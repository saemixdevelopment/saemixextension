#' Predictions for an SaemixModel object appiled to the predictors in a SaemixData object
#' 
#' @param object an SaemixModel object
#' @param psi a vector or a dataframe giving the parameters for which predictions are to be computed (defaults to empty). 
#' The number of columns in psi (or the number of elements of psi, if psi is given as a vector) should match the number of
#' parameters in the model, otherwise an error message will be shown and the function will return empty.
#' If psi is NA, the predictions are computed for the population parameters in the model (first line of the psi0 slot).
#' Covariates are not taken into account in the prediction. 
#' If psi is a dataframe, each line will be used for a separate 'subject' in the predictors dataframe, as 
#' indicated by the id argument; if id is not given, only the first line of psi will be used. 
#' @param predictors an SaemixData object (the predictors will then be extracted from the object using the name.predictors slot of the object)
#' @param id the vector of subjects for which individual plots will be obtained. If empty, the first 12 subjects in the dataset will be used (subject id's are taken from the name.group slot in the data object). If id is given, individual plots will be shown for the matching subjects in the dataset (eg if id=c(1:6), the first 6 subjects in the dataframe will be used for the plots, retrieving their ID from the data object)
#' 
#' @param \dots unused argument, for consistency with the generic
#' 
#' @details The function uses the model slot of the SaemixModel object to obtain predictions, using the predictors object. The
#' user is responsible for giving all the predictors needed by the model function.
#' if psi is not given, the predictions will be computed for the population parameters (first line of the psi0 slot) of the object.
#' 
#' @details The predictions correspond to the structure of the model; for models defined in terms of their likelihood, 
#' No individual graphs are currently available for discrete data models.
#' 
#' @details Warning: this function is currently under development and the output may change in future versions of the package 
#' following changes 
#' 
#' @seealso \code{\link[predict.SaemixModel]}, \code{\link[plotDiscreteData]},  \code{\link[ggplot]}
#' 
#' @return plots the data overlayed with the model predictions for each subject in id (where id is the index in the N subjects), and returns the plot invisibly (ggplot)
#' 
#' @examples 
#' data(theo.saemix)
#' xpred<-theo.saemix[,c("Dose","Time")]
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
#' saemix.model<-saemixModel(model=model1cpt,modeltype="structural",
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' checkInitialFixedEffects(saemix.model, saemix.data, id=c(1:6))
#' checkInitialFixedEffects(saemix.model, saemix.data, id=c(1:6), psi=c(0.5, 30, 2)) # better fit
#' 
#' @export 

checkInitialFixedEffects<-function(object, data, psi=c(), id=c(), ...) {
  if(!is(data,"SaemixData") ) {
    message("The data argument should be an SaemixData object to extract the predictors from.\n")
    return()
  }
  if(object@modeltype!="structural") {
    message("Individual plots are only available for models dealing with continuous outcomes.\n")
    return()
  }
  if(length(psi)==0) psi<-object["psi0"][1,,drop=FALSE]
  if(is.null(dim(psi))) psi<-as.data.frame(t(psi)) # psi given as a vector
  if(dim(psi)[2] != object@nb.parameters) {
    message(paste0("psi must have a number of columns equal to the number of parameters in the model (",object@nb.parameters,")\n")
    )
    return()
  }
  # Select subjects corresponding to number id, or use the first 12 subjects
  if(length(id)==0) id<-1:12
  idall<-data@data[,"index"]
  nind.obs <- data@nind.obs
  predictors <- data@data[,data@name.predictors]
  id<-intersect(idall, id) # unique indexes of subjects to use for plot
  if(length(id)==0) id<-1:12
  idkeep <- which(data@data$index %in% id) # retrieve data for these subjects
  xidep<-predictors[idkeep,]
  idx<-rep(c(1:length(id)),times=nind.obs[id])   # renumber id
  if(dim(psi)[1]==1 & length(unique(id))>1)
    psi<-do.call(rbind,rep(list(psi),length(unique(id))))
  colnames(psi)<-colnames(object["psi0"])
  rownames(psi)<-NULL
  # Model predictions
  ypred<-object["model"](psi, idx, xidep)
  obspl<-data.frame(id=data@data[idkeep,data@name.group], x=data@data[idkeep,data@name.X], y=data@data[idkeep,data@name.response])
  predpl<-data.frame(id=data@data[idkeep,data@name.group], x=data@data[idkeep,data@name.X], y=ypred)

  # Individual graphs
  myplot <- ggplot(data=obspl, aes(x=x, y=y, group=id)) + geom_point() + geom_line(data=predpl) + 
    xlab(paste0(data@name.X," (",data@units$x,")")) + ylab(paste0(data@name.response," (",data@units$y,")")) + 
    theme_bw() + facet_wrap(.~id, nrow=3, ncol=4)
  print(myplot)

  invisible(myplot)
}

