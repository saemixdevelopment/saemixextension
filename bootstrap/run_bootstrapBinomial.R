saemixDir<-"/home/eco/work/monolix/rversion/newLib/saemix"
setwd(saemixDir)
source(file.path(saemixDir,"testeco","helper-source.R"))
datDir<-file.path(saemixDir,"data")
resDir<-"/home/eco/work/monolix/rversion/newLib/bootstrapSim"

sirDir<-"/home/eco/work/theses/sirM2marilou/simulations"
sirDataDir<-file.path(sirDir,"data")

source(file.path(saemixDir,"bootstrap","caseBootstrap.R"))

binary.model<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  logpdf<-rep(0,length(tim))
  P.obs = (y==0)*(1-pevent)+(y==1)*pevent
  logpdf <- log(P.obs)
  return(logpdf)
}

saemix.model<-saemixModel(model=binary.model,description="Binary model",modeltype="likelihood",
            psi0=matrix(c(0,-.5,0,0.5),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
            transform.par=c(0,0), covariate.model=c(0,1),covariance.model=matrix(c(1,0,0,1),ncol=2),verbose=FALSE)

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE,displayProgress=FALSE,warnings=FALSE)

nboot<-200 # Nb of bootstrap samples

for(ires in 1:200) {
  namfile<-file.path(sirDataDir,paste("data",ires,".txt",sep=""))
  tab1<-read.table(namfile,header=TRUE,stringsAsFactors = FALSE)
  
  saemix.data<-saemixData(name.data=namfile,name.group=c("id"),name.predictors=c("time","y"), 
                          name.covariates=c("trt"),name.X=c("time"),verbose=FALSE)
  binary.fit<-saemix(saemix.model,saemix.data,saemix.options)
  saemix.fit<-binary.fit
  
  # nboot<-100 # Nb of samples from the conditional distribution
  source(file.path(saemixDir,"bootstrap","bootstrapDistribution.R"))
  namfile<-file.path(resDir,"binomial",paste("bootstrap",ires,".txt",sep=""))
  write.table(res.boot,namfile,quote=F,col.names=TRUE)
}
