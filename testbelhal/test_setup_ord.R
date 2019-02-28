###################################################################################
cat("Running example CAT\n")

smx.ord <- read.table(file.path(datDir,"categorical1_data.txt"),header=T)
saemix.data<-saemixData(name.data=smx.ord,header=TRUE,sep=" ",na=NA, 
name.group=c("ID"),name.predictors=c("Y"), 
name.X=c("TIME"))

ordinal.model<-function(psi,id,xidep) {
  y<-xidep[,1]
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  logit1<-alp1
  logit2<-alp1+alp2
  logit3<-alp1+alp2+alp3
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  logpdf<-rep(0,length(y))
  P.obs = (y==0)*pge1+(y==1)*(pge2 - pge1)+(y==2)*(pge3 - pge2)+(y==3)*(1 - pge3)
  logpdf <- log(P.obs)

  return(logpdf)
}


saemix.model<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
psi0=matrix(c(3,1,1),ncol=3,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3"))),
omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
transform.par=c(0,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3))
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

ord.fit<-saemix.fit

###################################################################################
# New dataset
xtim<-seq(1,3,1)
nsuj<-10
ord.newdata<-data.frame(ID=rep(1:nsuj,each=length(xtim)),TIME=rep(xtim,nsuj))

saemixObject<-saemix.fit
psiM<-data.frame(alp1=seq(2,2.9,0.1),alp2 = seq(0.6,1.5,0.1),alp3 = seq(0.3,1.2,0.1))

simul.ord<-function(psi,id,xidep) {
  y<-xidep
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  logit1<-alp1
  logit2<-alp1+alp2
  logit3<-alp1+alp2+alp3
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  p0 <- pge1
  p1 <- pge2 - pge1
  p2 <- pge3 - pge2
  p3 <- 1 - pge3
  obs <-rep(0,length(y))

  for (i in (1:length(obs))){
    obs[i] <- rmultinom(1, size = 3, prob = c(p1[i],p2[i],p3[i]))[1]
  }  

  return(obs)
}

preds <- simul.ord(psiM, ord.newdata$ID, ord.newdata[,c("TIME")])
ord.newdata$Y<-preds

ord.psiM<-psiM

###################################################################################
# test dataset

# test.newdata <- read.table(file.path(datDir,"categorical2_test_data.txt"),header=T)
# # test.newdata <- test.newdata[1:240,]
# test.newdata <- test.newdata[test.newdata$TIME<4,]
# saemixObject<-saemix.fit
# psiM<-data.frame(alp1=seq(2,2.19,0.01),alp2 = seq(0.6,0.79,0.01),alp3 = seq(0.1,0.29,0.01))
# # fpred <- ordinal.model(psiM, test.newdata$ID, test.newdata[,c("TIME","Y")])
# fpred<-saemixObject["model"]["model"](psiM, test.newdata$ID, test.newdata[,c("TIME","Y")])
# test.newdata$LogProbs<-fpred

# ord.newdata<-test.newdata
# ord.psiM<-psiM

###################################################################################
