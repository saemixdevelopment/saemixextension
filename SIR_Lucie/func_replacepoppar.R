####################################################################################
######      Replace the population parameters component in an saemixObject    ######
####################################################################################

#' Replace the population parameters elements in a SaemixObject object
#' 
#' Returns an SaemixObject object where the population parameters elements have been replaced by the population parameters elements provided in 3 matrix
#'
#' @name replacePopPar
#' @aliases replacePopPar-methods replacePopPar.saemixObject
#' 
#' @param object an saemixObject object
#' @param newfixed.par a vector containing new fixed population parameters
#' @param newomega a matrix containing new omega
#' @param newres.par a vector containing new residual error
#' @return an object of class \code{"\linkS4class{SaemixObject}"}. The data and model are retained but all the population parameters, predictions, individual parameters and statistical criteria are removed.
#' @keywords methods
#' @export 
#' 




replacePopPar.saemixObject<-function(saemixObject, newparam) {
  # Takes a saemixObject fit, changes the population parameters and removes all results pertaining to individual parameters and LL
  
  #Construction of pop parameters to change from a vector of parameters
  R <- estpar.fromvect.tomatrix(newparam, saemixObject)
  newfixed.par  <- R[['newfixed.effects']]
  newomega <- R[['newomega']]
  newres.par <- R[['newres']]
  
  saemix.newObj<-saemixObject

  saemix.newObj['results']['fixed.effects'] <- newfixed.par
  saemix.newObj['results']['omega'] <- newomega
  saemix.newObj['results']['respar'] <- newres.par
  
  saemix.newObj["results"]["cond.mean.phi"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["cond.mean.psi"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["cond.var.phi"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["cond.mean.eta"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["cond.shrinkage"] <-numeric(0)
  saemix.newObj["results"]["mean.phi"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["map.psi"] <-data.frame()
  saemix.newObj["results"]["map.phi"] <-data.frame()
  saemix.newObj["results"]["map.eta"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["map.shrinkage"] <-numeric(0)
  saemix.newObj["results"]["phi"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["psi.samp"] <-array(dim=c(0,0,0),data=0)
  saemix.newObj["results"]["phi.samp"] <-array(dim=c(0,0,0),data=0)
  saemix.newObj["results"]["phi.samp.var"] <-array(dim=c(0,0,0),data=0)
  saemix.newObj["results"]["ll.lin"] <-numeric(0)
  saemix.newObj["results"]["aic.lin"] <-numeric(0)
  saemix.newObj["results"]["bic.lin"] <-numeric(0)
  saemix.newObj["results"]["bic.covariate.lin"] <-numeric(0)
  saemix.newObj["results"]["ll.gq"] <-numeric(0)
  saemix.newObj["results"]["aic.gq"] <-numeric(0)
  saemix.newObj["results"]["bic.gq"] <-numeric(0)
  saemix.newObj["results"]["bic.covariate.gq"] <-numeric(0)
  
  saemix.model <- saemixObject["model"]
  saemix.newObj["results"]["nbeta.random"]<- sum(saemix.model["betaest.model"]%*%diag(saemix.model["fixed.estim"])%*%as.matrix(diag(saemix.model["covariance.model"])))
  saemix.newObj["results"]["nbeta.fixed"]<-  sum(saemix.model["betaest.model"]%*%diag(saemix.model["fixed.estim"])%*%as.matrix(-1*diag(saemix.model["covariance.model"])+1))
  
  saemix.newObj["results"]["ll.is"] <-numeric(0)
  saemix.newObj["results"]["aic.is"] <-numeric(0)
  saemix.newObj["results"]["bic.is"] <-numeric(0)
  saemix.newObj["results"]["bic.covariate.is"] <-numeric(0)
  saemix.newObj["results"]["predictions"] <-data.frame()
  saemix.newObj["results"]["ypred"] <-numeric(0)
  saemix.newObj["results"]["ppred"] <-numeric(0)
  saemix.newObj["results"]["ipred"] <-numeric(0)
  saemix.newObj["results"]["icpred"] <-numeric(0)
  saemix.newObj["results"]["ires"] <-numeric(0)
  saemix.newObj["results"]["iwres"] <-numeric(0)
  saemix.newObj["results"]["icwres"] <-numeric(0)
  saemix.newObj["results"]["wres"] <-numeric(0)
  saemix.newObj["results"]["npde"] <-numeric(0)
  saemix.newObj["results"]["pd"] <-numeric(0)
  
  saemix.newObj <- estimateMeanParametersNewdata(saemix.newObj)
  
  saemix.newObj <- try(estimateIndividualParametersNewdata(saemix.newObj))
  
  return(saemix.newObj)
  
}


####################################################################################