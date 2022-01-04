#' Predictions for a new dataset
#' 
#' @param object an SaemixModel object
#' @param psi a vector or a dataframe giving the parameters for which predictions are to be computed (defaults to NA). 
#' The number of columns in psi (or the number of elements of psi, if psi is given as a vector) should match the number of
#' parameters in the model, otherwise an error message will be shown and the function will return empty.
#' If psi is NA, the predictions are computed for the population parameters in the model (first line of the psi0 slot). 
#' If psi is a dataframe, each line will be used for a separate 'subject' in the predictors dataframe, as 
#' indicated by the id argument; if id is not given, only the first line of psi will be used. 
#' @param predictors a dataframe with the predictors for the model (must correspond to the predictors used by the model function)
#' @param id a vector of indices of length equal to the number of lines in predictors, matching each line of predictors to the 
#' corresponding line in psi, ie the parameters for this predictors (defaults to NA). If id is given, the unique values in id must be equal
#' to the number of lines in psi, otherwise id will be set to 1. If id is given and its values do not take the consecutive values 1:N, the
#' indices will be matched to 1:N to follow the lines in psi.
#' 
#' @details The function uses the model slot of the SaemixModel object to obtain predictions, using the predictors object. The
#' user is responsible for giving all the predictors needed by the model function.
#' if psi is not given, the predictions will be computed for the population parameters (first line of the psi0 slot) of the object.
#' 
#' @details The predictions correspond to the structure of the model; for models defined in terms of their likelihood, the predictions 
#' are the log-pdf of the model (see documentation for details).
#' 
#' @details Warning: this function is currently under development and the output may change in future versions of the package to conform to the usual predict functions.
#' 
#' @return a list with two components
#' @describe{
#' \item{param}{a dataframe with the estimated parameters}
#' \item{predictions}{a dataframe with the population predictions}
#' }
#' 
#' @aliases predict.SaemixModel
#' 
#' @examples 
#' # TODO
#' @export

predict.saemixmodel<-function(object, predictors, psi=NA, id=NA) {
  xidep<-predictors
  if(is.na(id) || length(id)!=dim(predictors)[1]) 
    id<-rep(1,dim(xidep)[1]) 
  idkeep<-id
  if(max(id)>length(unique(id))) { # indexes need to go from 1 to N
    id1<-1:length(unique(id))
    id2<-unique(id)
    id<-id1[match(id,id2)]
  }
  if(is.na(psi)) psi<-object["psi0"][1,,drop=FALSE]
  if(is.null(dim(psi))) psi<-as.data.frame(t(psi)) # psi given as a vector
  if(dim(psi)[2] != object@nb.parameters) {
      message(paste0("psi must have a number of columns equal to the number of parameters in the model (",object@nb.parameters,")\n")
      )
      return()
    }
  if(dim(psi)[1]==1 & length(unique(id))>1)
     psi<-do.call(rbind,rep(list(psi),length(unique(id))))
  colnames(psi)<-colnames(object["psi0"])
  rownames(psi)<-NULL
  ypred<-object["model"](psi, id, xidep)
  return(list(param=cbind(id=1:dim(psi)[1],psi), predictions=data.frame(id=idkeep, xidep, pred=unname(ypred))))
}

#' Plot model predictions for a new dataset. If the dataset is large, only the first 20 subjects (id's) will be shown.
#' 
#' @param smx.model an SaemixModel object
#' @param smx.data an SaemixData object
#' @param psi a vector or a dataframe giving the parameters for which predictions are to be computed (defaults to NA). 
#' The number of columns in psi (or the number of elements of psi, if psi is given as a vector) should match the number of
#' parameters in the model, otherwise an error message will be shown and the function will return empty.
#' If psi is NA, the predictions are computed for the population parameters in the model (first line of the psi0 slot). 
#' If psi is a dataframe, each line will be used for a separate subject of the smx.data object. Elements of psi will be recycled 
#' if psi has less lines than the number of subjects in the dataset.
#' 
#' @details The function uses the model slot of the SaemixModel object to obtain predictions, using the dataset contained in the 
#' SaemixData object. The user is responsible for making sure data and model match.
#' if psi is not given, the predictions will be computed for the population parameters (first line of the psi0 slot) of the object.
#' 
#' @details Currently this function only works for models defined as 'structural'.
#' 
#' @return a ggplot object
#' 
#' @aliases plot.SaemixModel
#' 
#' @examples 
#' # TODO
#' @export

plot.saemixModel <- function(smx.model, smx.data, psi=NA) {
  if(smx.model@modeltype!="structural") {
    message("Currently plots of the model are only available for continuous response models\n")
    return()
  }
  if(is.na(psi)) psi<-smx.model["psi0"][1,,drop=FALSE]
  if(is.null(dim(psi))) psi<-as.data.frame(t(psi)) # psi given as a vector
  if(dim(psi)[2] != smx.model@nb.parameters) {
    message(paste0("psi must have a number of columns equal to the number of parameters in the model (",smx.model@nb.parameters,")\n"))
    return()
  }
  if(dim(psi)[1]==1 || dim(psi)[1]<smx.data@N)
    psi<-do.call(rbind,rep(list(psi),length.out=smx.data@N))
  
  nvalues<-100
  xt<-seq(min(smx.data@data[,smx.data["name.X"]]), max(smx.data@data[,smx.data["name.X"]]), length.out=nvalues)
  xidep<-data.frame(x=xt)
  colnames(xidep)<-smx.data["name.X"]
  if(length(smx.data@name.predictors)>1) {
    id<-smx.data@data[,smx.data@name.group]
    otherpred<-smx.data@name.predictors[smx.data@name.predictors != smx.data["name.X"]]
    x1<-smx.data@data[match(unique(id), id), otherpred, drop=FALSE]
    dat1<-NULL
    for(i in 1:length(unique(id)))
      dat1<-rbind(dat1, 
                  do.call(rbind,rep(list(x1[i,,drop=FALSE]), nvalues)))
    xidep<-cbind(xidep, dat1)
    colnames(xidep[2:dim(xidep)[2]])<-otherpred
    xidep<-xidep[,smx.data["name.predictors"]] # Sort the predictors back in the correct order...
  }
  id<-rep(1:length(unique(id)), each=nvalues)
  y<-predict.saemixmodel(smx.model, predictors=xidep, psi=psi, id=id)
  gpred<-cbind(id=id,xidep,y=y$predictions$pred)
  colnames(gpred)[colnames(gpred)==smx.data@name.X]<-"x"
  
  gdat<-smx.data@data
  colnames(gdat)[colnames(gdat)==smx.data@name.X]<-"x"
  colnames(gdat)[colnames(gdat)==smx.data@name.response]<-"y"
  colnames(gdat)[colnames(gdat)==smx.data@name.group]<-"id"
  
  if(length(unique(gdat$id))>20) {
    nrow<-4
    ncol<-5
    zesuj<-unique(gdat$id)
    gdat1<-gdat[gdat$id %in% zesuj[1:12],]
    gpred1<-gpred[gpred$id %in% zesuj[1:12],]
  } else {
    gdat1<-gdat
    gpred1<-gpred
    nrow<-NULL
    ncol<-NULL
  }
  
  g1<-ggplot(data=gdat1, aes(x=x, y=y, group=id)) + geom_point() + geom_line(data=gpred1,aes(x=x, y=y)) + 
    facet_wrap(.~id, nrow=nrow, ncol=ncol) + labs(x=smx.data@name.X, y=smx.data@name.response) + theme_bw()
  return(g1)
} 

