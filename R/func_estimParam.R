
estimateMeanParameters.newdata<-function(saemixObject) {
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
estimateIndividualParameters.newdata<-function(saemixObject,type=c("mode","mean"),nsamp=1) {
  if(length(saemixObject["results"]["mean.phi"])==0) {
    if(saemixObject["options"]$warnings) cat("Population parameters (Ci*mu) will first be estimated to provide a starting point for the estimation of the individual parameters.\n")
    estimateMeanParameters.newdata(saemixObject)
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

predict.newdata<-function(saemixObject, saemix.newdata, type=c("ipred", "ypred", "ppred", "icpred"),nsamp=1) {
  # nsamp ignored for other types than icpred; if icpred, returns both the mean of the conditional distribution and samples, with the corresponding predictions
  # Replace the data object in saemixObject with the newdata, and wipe out the parameter and likelihood estimates associated with the initial run
  saemixObject<-replaceData.saemixObject(saemixObject,saemix.newdata)
  if(sum(is.na(saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]))>0) {
    if(saemixObject["model"]["modeltype"]=="likelihood") {
      if(saemixObject["options"]$warnings) cat("Please provide values of the response to obtain predictions for a model defined by loglikelihood\n")
      return(NULL)
    }
    type<-type[type!="ipred" & type!="icpred"]
  }
  if(length(type)==0) type<-"ypred"
  
  # Estimate population parameters (Ci*mu) for the new subjects
  saemixObject<-estimateMeanParameters.newdata(saemixObject)
  
  # Predictions using the mean parameters ypred=f(E(theta,x))
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
  predictions<-data.frame(IdM,XM,ypred=fpred)
  colnames(predictions)[1]<-newdata["name.group"]
  parameters<-list(id=unique(newdata["data"][,newdata["name.group"]]), population=psiM)
  
  # Mean predictions over the population ppred=E(f(theta,x))
  # Technically... we can obtain ppred as E(f()) by simulating etas in their distribution
  if(length(grep(c("ppred"),type))>0) {
    ind.eta<-saemixObject["model"]["indx.omega"]
    nb.etas<-length(ind.eta)
    omega<-saemixObject["results"]["omega"]
    chol.omega<-try(chol(omega[ind.eta,ind.eta]),silent=TRUE)
    
    etaM<-matrix(data=0,nrow=NM,ncol=nb.etas)
    ppred<-matrix(data=0,nrow=dim(XM)[1],ncol=saemixObject["options"]$nb.sim)
    mean.phiM<-do.call(rbind,rep(list(mean.phi),1))
    phiMc<-mean.phiM
    if(length(grep(c("ppred"),type))==1) {
      for(isim in 1:saemixObject["options"]$nb.sim) {
        etaMc<-0.5*matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
        phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
        psiMc<-transphi(phiMc,saemixObject["model"]["transform.par"])
        fpred<-saemixObject["model"]["model"](psiMc, IdM, XM)
        ppred[,isim]<-fpred
      }
      ppred<-rowMeans(ppred)
      predictions$ppred<-ppred
    }
  }
  if(sum(!is.na(match(c("icpred","ipred"),type)))==0) { # only population predictions
    return(list(param=parameters,predictions=predictions))
  }
  # Estimate individual parameters, if type contains ipred and/or icpred
  # Call predict instead ? (but will compute both map and cond.mean)
  ctype<-c()
  if(length(grep("ipred",type))==1) ctype<-c(ctype,"mode")
  if(length(grep("icpred",type))==1) ctype<-c(ctype,"mean")
  saemixObject<-estimateIndividualParameters.newdata(saemixObject,type=ctype,nsamp=nsamp)
  
  if(length(grep("icpred",type))==1) {
    psiM<-parameters$cond.mean.psi<-saemixObject["results"]["cond.mean.psi"]
    parameters$cond.var.phi<-saemixObject["results"]["cond.var.phi"]
    parameters$cond.mean.phi<-saemixObject["results"]["cond.mean.phi"]
    fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
    predictions<-cbind(predictions,icpred=fpred)
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
  }
  saemixObject["results"]["predictions"]<-predictions
  rlist<-list(param=parameters,predictions=predictions, object=saemixObject)
  if(length(grep("icpred",type))==1 & nsamp>1) rlist<-list(param=parameters,predictions=predictions, object=saemixObject, parSample=samp.par, predSample=samp.pred)

  return(rlist)
}

