# ToDo: adjust ind.ioM (we might need indicators specific to each level of variability)

########################### Functions used to initialise the lists in the algorithm

## Creation of Dargs 
### contains data elements - passed on to functions, unchanged
### returns Dargs and DYF (rectangular matrix of size max(nobs_i)xNM with NM=N x nb.chains)
### previously also contained elements relative to the model:  
#### (etype.exp=which(errmod == "exponential"), modeltype="structural",transform.par=saemix.model["transform.par"], error.model=saemix.model["error.model"],structural.model=structural.model)
### Question: keep here, add later (after call to createStructData) or move to another list ? (error model already associated with outcome, parameter transformation has associated functions in the model object)

createStructData <- function(saemix.data, nb.chains=1) {
  # Input
  ## SaemixData object
  # Output
  ## DYF: rectangular matrix max(nobs_i) x NM to keep track of LL(y_ijk)
  ## Dargs: data elements - passed on to functions, unchanged
  ### N (nb of subjects), nobs (nb of observations), yobs (observations)
  ### multiple chains: NM (nb of subjects, replicated if nb.chains>1), XM (predictors), yM (observations), IdM (match subject-lines in the dataset)
  ## added: indicators
  ### ind.ioM: non-null elements in DYF (in columns)
  ### is.lpdf=which elements are returned as log-pdf (corresponding to discrete outcomes)
  
  # removed from Dargs, to add in Margs (model elements)
  ### (transform.par=saemix.model["transform.par"], error.model=saemix.model["error.model"],structural.model=structural.model)
  
  # using several Markov chains
  chdat<-new(Class="SaemixRepData",data=saemix.data, nb.chains=nb.chains)
  NM<-chdat["NM"]
  IdM<-chdat["dataM"]$IdM
  yM<-chdat["dataM"]$yM
  XM<-chdat["dataM"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
  N<-saemix.data["N"]
  io<-matrix(data=0,nrow=N,ncol=max(saemix.data["nind.obs"]))
  for(i in 1:N)
    io[i,1:saemix.data["nind.obs"][i]]<-1
  ioM<-do.call(rbind,rep(list(io),nb.chains))
  ind.ioM <- which(t(ioM)!=0)
  DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])
  is.lpdf<-rep(0,saemix.data["ntot.obs"])
  which.lpdf<-which(saemix.data@outcome!="continuous")
  if(length(which.lpdf)>0) 
    is.lpdf[saemix.data@data$ytype %in% which.lpdf]<-1
  is.lpdfM<-rep(is.lpdf,nb.chains)
  Dargs<-list(IdM=IdM, XM=XM, yM=yM, NM=NM, N=N, nb.chains=nb.chains,
              nobs=saemix.data["ntot.obs"],yobs=saemix.data["data"][,saemix.data["name.response"]], 
              is.lpdfM=is.lpdfM,ind.ioM=ind.ioM)
  return(list(Dargs=Dargs, DYF=DYF))
}

createStructModel <- function(saemix.model) {
  # Input
  ## SaemixModel object
  # Output
  ## Margs: model elements - passed on to functions, unchanged
  ### structural model, list of parameter transformation
  ### model.type: type of outcome (continuous vs discrete)
  ### error.model: list of error model objects (or discrete)
  model.type<-error.model <- vector(mode="list",length=saemix.model@noutcome)
  for(i in 1:saemix.model@noutcome)
    if(saemix.model@outcome[[i]]@type.outcome=="continuous") {
      error.model[[i]]<-saemix.model@outcome[[i]]@error.model
      model.type[[i]]<-"continuous"
      } else error.model[[i]]<-model.type[[i]]<-"discrete"
    # Model related list
    Margs <- list(transform = saemix.model@transform,
                  invtransform = saemix.model@invtransform,
                  structural.model=saemix.model@model, 
                  error.model = error.model,
                  model.type = model.type)
    return(list(Margs=Margs))
}

createVarlist <- function(saemix.model, saemix.data) {
  # Margs: model elements - passed on to functions, unchanged
  # Vargs: list containing parameters updated during the fit (fixed and random)
  ## question: also design matrices ?
  # need to recover the number of occasions in the dataset
  ## needs both model and data so should be done when matching data and model
  nvarlev <- length(saemix.model@varlevel)
  nocc.max <- length(unique(saemix.data@data[,saemix.data@varlevel[nvarlev]]))
  
  # parameters with at least one level of variability
  ind.eta<-c()
  for(ivar in 1:length(saemix.model@varlevel)) ind.eta<-c(ind.eta,saemix.model@varmodel[[ivar]]@idcol.eta)
  ind.eta<-sort(unique(ind.eta))
  ind.eta.ext<-c()
  for(i in 1:nocc.max) ind.eta.ext <- c(ind.eta.ext, ind.eta+(i-1)*nocc.max)

  # Variance parameters
  omega <- list(saemix.model@varmodel[[1]]@omega)
  # with 2 levels of variability
  ## extension to more levels ?
  if(nvarlev==1) omega.sim<-omega[[1]] else {
    omega[[2]] <- saemix.model@varmodel[[2]]@omega
    omega.sim <- kronecker(matrix(1,nocc.max,nocc.max), omega[[1]]) + kronecker(diag(1,nocc.max,nocc.max), omega[[2]])
  }
  # design matrices
  
  # fixed effect parameters (mu and beta)
    
  Vargs <- list(omega=omega, omega.sim=omega.sim)
  return(list(Vargs=Vargs, ind.eta=ind.eta, ind.eta.ext=ind.eta.ext))
}

########################### Functions used to initialise the individual parameters
