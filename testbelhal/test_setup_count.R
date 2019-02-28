###################################################################################
cat("Running example Rcount\n")

count.saemix<-read.table(file.path(datDir,"count1_data.txt"),header=T)
# count.saemix <- count.data$Y
saemix.data<-saemixData(name.data=count.saemix,header=TRUE,sep=" ",na=NA, name.group=c("ID"),
  name.predictors=c("Y"),name.response=c("Y"))

#Basic Poisson model
countmodel<-function(psi,id,xidep) { 
  y<-xidep[,1]
  lambda<-psi[id,1]
  dummy<-psi[id,2]
  logp <- -lambda + y*log(lambda) - factorial(y)
  return(logp)
}

saemix.model<-saemixModel(model=countmodel,description="count model",modeltype="likelihood",   
  psi0=matrix(c(0.5,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","dummy"))), 
  transform.par=c(1,1),omega.init=matrix(c(0.3,0,0,0.3),ncol=2,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE),fixed.estim=c(1,1))

# ##Generalized Poisson model
# countmodel<-function(psi,id,xidep) { 
#   y<-xidep[,1]
#   delta<-psi[id,1]
#   lambda<-psi[id,2]
#   logp <- -lambda
#   pos.ind <- which(y>0)
#   logp[pos.ind] <- log(lambda) + (y-1)*log(lambda+y*delta) - (lambda+y*delta) - factorial(y)
#   return(logp)
# }


# saemix.model<-saemixModel(model=countmodel,description="count model",modeltype="likelihood",   
#   psi0=matrix(c(0.5,0.5),ncol=2,byrow=TRUE,dimnames=list(NULL,   
#   c("delta","lambda"))), 
#   transform.par=c(1,1),omega.init=matrix(c(0.3,0,0,0.3),ncol=2,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,1),ncol=2, 
#   byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE,nbiter.sa=50)
saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

count.fit<-saemix.fit

###################################################################################
# New dataset
test.newdata<-read.table(file.path(datDir,"count2_data.txt"),header=T)
saemixObject<-count.fit

psiM<-data.frame(lambda=seq(0.1,0.5,length.out=length(unique(test.newdata$ID))),dummy = seq(1,3,4))
fpred<-saemixObject["model"]["model"](psiM, test.newdata$ID, test.newdata[,c("Y"),drop=FALSE])
test.newdata$LogProbs<-fpred
test.newdata$Y<-ifelse(test.newdata$TIME>0,1,0)

count.newdata<-test.newdata
count.psiM<-psiM

###################################################################################
