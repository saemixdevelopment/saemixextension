#' Predictions for a new dataset
#' 
#' @param saemixObject an SaemixObject from a fitted run
#' @param newdata a dataframe containing the new data. The dataframe must contain the same information as the original dataset (colunm names, etc...) 
#' @param type one or several of "ipred" (individual predictions using the MAP estimates), "ppred" (population predictions obtained using the population parameters f(E(theta))), "ypred" (mean of the population predictions (E(f(theta)))), "icpred"  (individual predictions using the conditional mean estimates). Defaults to "ppred".
#' @param nsamp an integer, ignored for other types than icpred; if icpred, returns both the mean of the conditional distribution and nsamp samples, with the corresponding predictions. Defaults to 1.
#' 
#' @details This function is the workhorse behind the predict method for SaemixObject. It computes predictions for a new dataframe
#' based on the results of an saemix fit. 
#' Since the predict function only returns predicted values, this function is provided so that users can access other elements, 
#' for example the different types of parameter estimates (such as individual or population parameters) associated with the predictions
#' For other purposes such as simply obtaining model predictions, we suggest using the predict() method with or without newdata.
#' 
#' @details The function uses estimateMeanParametersNewdata() to set the population estimates for the individual parameters
#' taking into account the individual covariates and doses, and estimateIndividualParametersNewdata() to derive individual estimates 
#' by computing the mean of the conditional distributions (type="icpred") or the MAP estimate (type="ipred")
#' 
#' @return a list with three components (five if type includes "icpred" and nsamp>1)
#' \describe{
#' \item{param}{a dataframe with the estimated parameters. The columns in the dataframe depend on which type of predictions were requested (argument type)}
#' \item{predictions}{a dataframe with the predictions. The columns in the dataframe depend on which type of predictions were requested (argument type)}
#' \item{saemixObject}{the SaemixObject with the data slot replaced by the new data. The elements of the results slot pertaining to individual (including population individual parameters) predictions and likelihood will have been removed, and only the elements computed within the function will have been replaced (eg individual estimated parameters and predictions for the new data)}
#' \item{parSample}{a dataframe with parameters sampled from the conditional distribution of the individual parameters (only present if type includes 'icpred' and nsamp>1)}
#' \item{predSample}{a dataframe with the predictions corresponding to the parameters sampled from the conditional distribution of the individual parameters (only present if type includes 'icpred' and nsamp>1)}
#' }
#' 
#' @aliases estimateMeanParametersNewdata estimateIndividualParametersNewdata saemixPredictNewdata
#' 
#' @examples 
#' # TODO
#' @export
 
saemixPredictNewdata<-function(saemixObject, newdata, type=c("ipred", "ypred", "ppred", "icpred"),nsamp=1) {
  # Predictions corresponding to a model fit for newdata
  ## replace the data object in saemixObject with the newdata
  ## wipe out the individual parameter and likelihood estimates associated with the initial run
  ## estimate the parameters and predictions for the new data with the results object
  # Inputs
  ## saemixObject: a fitted saemix object
  ## newdata: a dataframe containing the predictors needed by the model
  ## type: a vector of strings indicating which types of predictions should be obtained
  ## nsamp ignored for other types than icpred; for icpred, returns both the mean of the conditional distribution and samples, with the corresponding predictions
  # Returns a list with the following elements
  ## param: estimated parameters for the subjects in newdata
  ## predictions: dataframe containing the predictions
  ## object: the modified saemixObject with newdata replacing saemixObject@data
  ## if type="icpred" and nsamp>1, the list has two additional elements
  ### parSample: sampled parameters (size nsamp)
  ### predSample: predictions corresponding to the sampled parameters
  saemixObject<-replaceData.saemixObject(saemixObject,newdata)
  if(sum(is.na(saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]))>0) {
    if(saemixObject["model"]["modeltype"]=="likelihood") {
      if(saemixObject["options"]$warnings) cat("Please provide values of the response to obtain predictions for a model defined by loglikelihood\n")
      return(NULL)
    }
    type<-type[type!="ipred" & type!="icpred"]
  }
  if(length(type)==0) type<-"ppred"
  
  # Estimate population parameters (Ci*mu) for the new subjects
  saemixObject<-estimateMeanParametersNewdata(saemixObject) # updates mean.phi (normally...)
  
  # Predictions using the mean parameters ppred=f(E(theta,x))
  newdata<-saemixObject["data"]
  chdat<-saemixObject["rep.data"]
  NM<-chdat["NM"]
  IdM<-chdat["dataM"]$IdM
  yM<-chdat["dataM"]$yM
  XM<-chdat["dataM"][,c(newdata["name.predictors"],newdata["name.cens"],newdata["name.mdv"],newdata["name.ytype"]),drop=FALSE]
  mean.phi<-saemixObject["results"]["mean.phi"]
  psiM<-transphi(mean.phi,saemixObject["model"]["transform.par"])
  fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
  colnames(psiM)<-saemixObject["model"]["name.modpar"]
  predictions<-data.frame(IdM,XM,ppred=fpred)
  colnames(predictions)[1]<-newdata["name.group"]
  parameters<-list(id=unique(newdata["data"][,newdata["name.group"]]), population=psiM)
  saemixObject["results"]["ppred"]<-fpred

  # Mean predictions over the population ypred=E(f(theta,x))
  # Technically... we can obtain ypred as E(f()) by simulating etas in their distribution
  if(length(grep(c("ypred"),type))>0) {
    ind.eta<-saemixObject["model"]["indx.omega"]
    nb.etas<-length(ind.eta)
    omega<-saemixObject["results"]["omega"]
    chol.omega<-try(chol(omega[ind.eta,ind.eta]),silent=TRUE)
    
    etaM<-matrix(data=0,nrow=NM,ncol=nb.etas)
    ypred<-matrix(data=0,nrow=dim(XM)[1],ncol=saemixObject["options"]$nb.sim)
    mean.phiM<-do.call(rbind,rep(list(mean.phi),1))
    phiMc<-mean.phiM
    if(length(grep(c("ypred"),type))==1) {
      for(isim in 1:saemixObject["options"]$nb.sim) {
        etaMc<-0.5*matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
        phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
        psiMc<-transphi(phiMc,saemixObject["model"]["transform.par"])
        fpred<-saemixObject["model"]["model"](psiMc, IdM, XM)
        ypred[,isim]<-fpred
      }
      ypred<-rowMeans(ypred)
      predictions$ypred<-ypred
      saemixObject["results"]["ypred"]<-ypred
    }
  }
  if(sum(!is.na(match(c("icpred","ipred"),type)))==0) { # only population predictions
    return(list(param=parameters,predictions=predictions, object=saemixObject))
  }
  # Estimate individual parameters, if type contains ipred and/or icpred
  # Call predict instead ? (but will compute both map and cond.mean)
  ctype<-c()
  if(length(grep("ipred",type))==1) ctype<-c(ctype,"mode")
  if(length(grep("icpred",type))==1) ctype<-c(ctype,"mean")
  saemixObject<-estimateIndividualParametersNewdata(saemixObject,type=ctype,nsamp=nsamp) # updates cond.mean.psi, map.psi (normally...)
  
  if(length(grep("icpred",type))==1) {
    psiM<-parameters$cond.mean.psi<-saemixObject["results"]["cond.mean.psi"]
    parameters$cond.var.phi<-saemixObject["results"]["cond.var.phi"]
    parameters$cond.mean.phi<-saemixObject["results"]["cond.mean.phi"]
    fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
    predictions<-cbind(predictions,icpred=fpred)
    saemixObject["results"]["icpred"]<-fpred
    if(nsamp>1) {
      samp.pred<-array(dim=c(length(fpred),nsamp))
      samp.par<-array(dim=c(dim(psiM),nsamp))
      for(isamp in 1:nsamp) {
        phiM<-saemixObject["results"]["phi.samp"][,,isamp]
        psiM<-samp.par[,,isamp]<-transphi(phiM,saemixObject["model"]["transform.par"])
        fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
        samp.pred[,isamp]<-fpred
      }
    }
  }
  if(length(grep("ipred",type))==1) {
    psiM<-parameters$map.psi<-saemixObject["results"]["map.psi"]
    fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
    predictions<-cbind(predictions,ipred=fpred)
    saemixObject["results"]["ipred"]<-fpred
  }
  saemixObject["results"]["predictions"]<-predictions
  rlist<-list(param=parameters,predictions=predictions, object=saemixObject)
  if(length(grep("icpred",type))==1 & nsamp>1) 
    rlist<-list(param=parameters,predictions=predictions, object=saemixObject, parSample=samp.par, predSample=samp.pred)
  
  return(rlist)
}


estimateMeanParametersNewdata<-function(saemixObject) {
  nb.chains<-1
  saemix.newdata<-saemixObject["data"]
  chdat<-new(Class="SaemixRepData",data=saemix.newdata, nb.chains=nb.chains)
  NM<-chdat["NM"]
  IdM<-chdat["dataM"]$IdM
  yM<-chdat["dataM"]$yM
  XM<-chdat["dataM"][,c(saemix.newdata["name.predictors"],saemix.newdata["name.cens"],saemix.newdata["name.mdv"],saemix.newdata["name.ytype"]),drop=FALSE]
  io<-matrix(data=0,nrow=saemix.newdata["N"],ncol=max(saemix.newdata["nind.obs"]))
  for(i in 1:saemix.newdata["N"])
    io[i,1:saemix.newdata["nind.obs"][i]]<-1
  ioM<-do.call(rbind,rep(list(io),nb.chains))
  ind.ioM <- which(t(ioM)!=0)
  saemixObject["rep.data"]<-chdat
  
  id<-saemix.newdata["data"][,saemix.newdata["name.group"]]
  if(length(saemix.newdata["name.covariates"])==0) tab<-data.frame(id=id) else
    tab<-data.frame(id=id,saemix.newdata["data"][, saemix.newdata["name.covariates"],drop=FALSE])
  temp<-tab[!duplicated(id),,drop=FALSE]
  
  if(length(saemix.newdata["name.covariates"])>0) {
    Mcovariates<-data.frame(id=rep(1,saemix.newdata["N"]),temp[,2:dim(temp)[2],drop=F])
  } else {
    Mcovariates<-data.frame(id=rep(1,saemix.newdata["N"]))
  }
  #namcov<-c(saemix.newdata["name.group"],rownames(saemixObject["model"]["betaest.model"]))
  namcov<-c("id",rownames(saemixObject["model"]["betaest.model"]))
  namcov<-namcov[!is.na(namcov)]
  Mcovariates<-Mcovariates[,namcov,drop=F]
  nb.parameters<-saemixObject["model"]["nb.parameters"]
  nb.betas<-sum(saemixObject["model"]["betaest.model"])
  ind.eta<-saemixObject["model"]["indx.omega"]
  nb.etas<-length(ind.eta)
  pop.par<-saemixObject["results"]["fixed.effects"]
  pop.par[saemixObject["results"]["indx.fix"]]<-transpsi(t(as.matrix(pop.par[saemixObject["results"]["indx.fix"]])),saemixObject["model"]["transform.par"])
  omega<-saemixObject["results"]["omega"]
  chol.omega<-try(chol(omega[ind.eta,ind.eta]),silent=TRUE)
  
  # Initialisation of phiM for new subjects
  # mean value for all would be phiM + Mcov*beta
  
  pfix<-matrix(data=0,nrow=1,ncol=nb.parameters)
  LCOV<-MCOV<-matrix(data=0,nrow=nb.betas,ncol=nb.parameters)
  j1<-1
  COV<-matrix(nrow=dim(Mcovariates)[1],ncol=0)
  mean.phi<-matrix(data=0,nrow=saemix.newdata["N"],ncol=nb.parameters)
  fixed.ini<-saemixObject["model"]["betaest.model"]*0
  fixed.ini[saemixObject["model"]["betaest.model"]==1]<-pop.par
  for(j in 1:nb.parameters) {
    jcov<-which(saemixObject["model"]["betaest.model"][,j]==1)
    lambdaj<-fixed.ini[jcov,j]
    aj<-as.matrix(Mcovariates[,jcov])
    COV<-cbind(COV,aj)
    nlj<-length(lambdaj)
    j2<-j1+nlj-1
    LCOV[j1:j2,j]<-matrix(data=1,nrow=nlj,ncol=1)
    j1<-j2+1
    if(length(jcov)<=1) mean.phi[,j]<-aj*lambdaj else mean.phi[,j]<-aj%*%lambdaj
    pfix[j]<-length(lambdaj)
  }
  saemixObject["results"]["mean.phi"]<-mean.phi
  return(saemixObject)
  
}
estimateIndividualParametersNewdata<-function(saemixObject,type=c("mode","mean"),nsamp=1) {
  if(length(saemixObject["results"]["mean.phi"])==0) {
    if(saemixObject["options"]$warnings) cat("Population parameters (Ci*mu) will first be estimated to provide a starting point for the estimation of the individual parameters.\n")
    estimateMeanParametersNewdata(saemixObject)
  }
  
  # Initialisation as mean.phi + eta
  chdat<-saemixObject["rep.data"]
  NM<-chdat["NM"]
  IdM<-chdat["dataM"]$IdM
  yM<-chdat["dataM"]$yM
  XM<-chdat["dataM"][,c(saemixObject["data"]["name.predictors"],saemixObject["data"]["name.cens"],saemixObject["data"]["name.mdv"],saemixObject["data"]["name.ytype"]),drop=FALSE]
  nb.parameters<-saemixObject["model"]["nb.parameters"]
  nb.betas<-sum(saemixObject["model"]["betaest.model"])
  ind.eta<-saemixObject["model"]["indx.omega"]
  nb.etas<-length(ind.eta)
  pop.par<-saemixObject["results"]["fixed.effects"]
  pop.par[saemixObject["results"]["indx.fix"]]<-transpsi(t(as.matrix(pop.par[saemixObject["results"]["indx.fix"]])),saemixObject["model"]["transform.par"])
  omega<-saemixObject["results"]["omega"]
  chol.omega<-try(chol(omega[ind.eta,ind.eta]),silent=TRUE)
  
  kt<-0
  mean.phi<-saemixObject["results"]["mean.phi"]
  mean.phiM<-do.call(rbind,rep(list(mean.phi),1))
  itest.phi<-1:NM
  ltest.phi<-length(itest.phi)
  phiM<-matrix(data=0,nrow=NM,ncol=nb.parameters)
  etaM<-matrix(data=0,nrow=NM,ncol=nb.etas)
  phiMc<-mean.phiM
  
  while (ltest.phi>0) {
    kt<-kt+1
    if (kt==100) 
      stop("Failed to find a valid initial parameter guess for the new data\n")
    end   
    etaMc<-0.5*matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
    phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
    etaM[itest.phi,]<-etaMc[itest.phi,]
    phiM[itest.phi,]<-phiMc[itest.phi,]
    psiM<-transphi(phiM,saemixObject["model"]["transform.par"])
    fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
    inan<-(is.na(fpred)+is.infinite(fpred)+(Im(fpred)!=0))
    itest.phi<-unique(IdM[inan])
    ltest.phi<-length(itest.phi)
  }
  saemixObject["results"]["phi"]<-phiM
  
  if(length(grep(c("mean"),type))==1) {
    # Estimating the conditional distribution of individual parameters for the new subjects
    saemixObject["results"]["cond.mean.phi"]<-phiM
    saemixObject<-conddist.saemix(saemixObject,nsamp=nsamp)
    saemixObject["results"]["cond.mean.psi"]<-transphi(saemixObject["results"]["cond.mean.phi"],saemixObject["model"]["transform.par"])
  }
  if(length(grep(c("mode"),type))==1) {
    saemixObject<-map.saemix(saemixObject)
  }
  return(saemixObject)
}



