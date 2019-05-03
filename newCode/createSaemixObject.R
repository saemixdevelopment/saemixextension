createSaemixObject.empty<-function(model,data,control=list()) {
  if(class(model)!="SaemixModel") {
    cat("Please provide a valid model object (see the help page for SaemixModel)\n")
    return()
  }
  if(class(data)!="SaemixData") {
    cat("Please provide a valid data object (see the help page for SaemixData)\n")
    return()
  }
  
  saemixObject<-new(Class="SaemixObject",data=data,model=model,options=control)
  #  saemixObject<-new(Class="SaemixObject",data=saemix.data, model=saemix.model,options=saemix.options)
  opt.warn<-getOption("warn")
  if(!saemixObject["options"]$warnings) options(warn=-1)
  
  saemix.data<-saemixObject["data"]
  saemix.data@ocov<-saemix.data@ocov[saemix.data@data[,"mdv"]==0,,drop=FALSE]
  saemix.data@data<-saemix.data@data[saemix.data@data[,"mdv"]==0,]
  saemix.data@ntot.obs<-dim(saemix.data@data)[1]
  saemixObject["data"]<-saemix.data
  
  return(saemixObject)
}

createSaemixObject.initial<-function(model,data,control=list()) {
  # Checking validity of input
    if(class(model)!="SaemixModel") {
    cat("Please provide a valid model object (see the help page for SaemixModel)\n")
    return()
  }
  if(class(data)!="SaemixData") {
    cat("Please provide a valid data object (see the help page for SaemixData)\n")
    return()
  }
  # Creating saemixObject (empty results)
  saemixObject<-new(Class="SaemixObject",data=data,model=model,options=control)
  opt.warn<-getOption("warn")
  if(!saemixObject["options"]$warnings) options(warn=-1)
  
  saemix.options<-saemixObject["options"]
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.data@ocov<-saemix.data@ocov[saemix.data@data[,"mdv"]==0,,drop=FALSE]
  saemix.data@data<-saemix.data@data[saemix.data@data[,"mdv"]==0,]
  saemix.data@ntot.obs<-dim(saemix.data@data)[1]

  # Initialising results component to initial estimates
  xinit<-initialiseMainAlgo(saemix.data,saemix.model,saemix.options)
  saemix.model<-xinit$saemix.model
  Uargs<-xinit$Uargs
  varList<-xinit$varList

  xres1<-new(Class="SaemixRes",modeltype="structural",status="initial",
        name.fixed=saemix.model["name.fixed"], name.random=saemix.model["name.random"],name.sigma=saemix.model["name.sigma"],
        fixed.effects=saemix.model@psi0[saemix.model@betaest.model==1],
        fixed.psi=xinit$fixedpsi.ini,
        betaC=xinit$betas[xinit$Uargs$indx.betaC],betas=xinit$betas,
        omega=varList$omega,respar=varList$pres,MCOV=varList$MCOV)
  xres1@indx.cov<-saemix.model@indx.cov
  xres1@indx.res<-saemix.model@indx.res
  xres1@indx.fix<-saemix.model@indx.fix
  xres1@indx.omega<-saemix.model@indx.omega
  saemixObject["results"]<-xres1
  
  return(saemixObject)
}
