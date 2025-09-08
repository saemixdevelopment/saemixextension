
## Creation of Dargs 
### contains data elements - passed on to functions, unchanged
### previously also contained elements relative to the model:  
#### (etype.exp=which(errmod == "exponential"), modeltype="structural",transform.par=saemix.model["transform.par"], error.model=saemix.model["error.model"],structural.model=structural.model)
### Question: keep here, add later (after call to createStructData) or move to another list ? (error model already associated with outcome, parameter transformation has associated functions in the model object)

createStructData <- function(saemix.data, nb.chains=1) {
  ### Dargs: data elements - passed on to functions, unchanged
  # Elements of the lists
  #   Dargs<-list(IdM=IdM, XM=XM, yM=yM, NM=NM, N=N, nobs=saemix.data["ntot.obs"],
  #               yobs=saemix.data["data"][,saemix.data["name.response"]], is.lpdf=is.lpdf)
  # removed: 
  ### (transform.par=saemix.model["transform.par"], error.model=saemix.model["error.model"],structural.model=structural.model)
  ### decide where to put these
  
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
  Dargs<-list(IdM=IdM, XM=XM, yM=yM, NM=NM, N=N, nobs=saemix.data["ntot.obs"],
              yobs=saemix.data["data"][,saemix.data["name.response"]], is.lpdfM=is.lpdfM)
  return(list(Dargs=Dargs, DYF=DYF, ind.ioM=ind.ioM))
}

