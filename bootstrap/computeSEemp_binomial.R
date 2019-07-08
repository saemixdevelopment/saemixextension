saemixDir<-"/home/eco/work/monolix/rversion/newLib/saemix"
setwd(saemixDir)
scenario<-"n200"
scenario<-"n20"

if(scenario=="n200") {
  sirDir<-"/home/eco/work/theses/sirM2marilou/simulations"
  namres<-file.path(saemixDir,"bootstrap","resEstimBinomialMarilou.tab")
}

if(scenario=="n20") {
  sirDir<-"/home/eco/work/theses/sirM2marilou/simulations/binom20"
  namres<-file.path(saemixDir,"bootstrap","resEstimBinomialMarilou_N20.tab")
}

sirDataDir<-file.path(sirDir,"data")

source(file.path(saemixDir,"testeco","helper-source.R"))

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

saemix.model<-saemixModel(model=binary.model,description="Binary model",
                          modeltype="likelihood",
                          psi0=matrix(c(0,-.5,0,0.5),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                          transform.par=c(0,0), covariate.model=c(0,1),covariance.model=matrix(c(1,0,0,1),ncol=2))

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE)

tabres<-NULL
l1<-c("Simulation","theta1","theta2","beta","omega2.theta1","omega2.theta2")
write(l1,namres,ncolumns = length(l1))

for(ires in 1:200) {
  namfile<-file.path(sirDataDir,paste("data",ires,".txt",sep=""))
  tab1<-read.table(namfile,header=TRUE,stringsAsFactors = FALSE)
  saemix.data<-saemixData(name.data=namfile,name.group=c("id"),name.predictors=c("time","y"), 
                          name.covariates=c("trt"),name.X=c("time"))
  
  binary.fit<-saemix(saemix.model,saemix.data,saemix.options)
  
  l1<-c(ires,binary.fit@results@fixed.effects,diag(binary.fit@results@omega))
  tabres<-rbind(tabres,l1)
  write(l1,namres,ncolumns = length(l1), append=TRUE)
}
