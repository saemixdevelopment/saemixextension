

createStructData(saemix.data, nb.chains=saemix.options$nb.chains)

# list for Phi ?
## length varlevel
## first = N, second Sum_i nocc_i, etc..

varlevel <- c("iiv","iov") # model variability levels
varcol <- c("id","occ") # columns in data corresponding to variability levels 

xphi <- vector(mode="list", length=length(varlevel))
names(xphi)<-varlevel

ilev <- 2
lab.ilev <- saemix.data@data[,saemix.data[varcol[1]]]
for(il in 2:ilev) lab.ilev <- paste(lab.ilev,saemix.data@data[,saemix.data[varcol[il]]],sep="|")
index.phi <- tapply(lab.ilev, lab.ilev, length) # nb of observations per id at this level [used to match psi to xidep]
nphilev <- length(!duplicated(lab.ilev)) # nb of parameters at this level (nb of rows of phi)



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