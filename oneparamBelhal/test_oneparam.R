# library(saemixextension)
source('R/aaa_generics.R') 
source('R/compute_LL.R') 
source('R/func_aux.R') 
source('R/func_distcond.R') 
source('R/func_FIM.R')
source('R/func_plots.R') 
source('R/func_simulations.R') 

source('R/main.R')
source('R/main_estep.R')
source('R/main_initialiseMainAlgo.R') 
source('R/main_mstep.R') 
source('R/SaemixData.R')
source('R/SaemixModel.R') 
source('R/SaemixRes.R') 
# source('R/SaemixRes_c.R') 
source('R/SaemixObject.R') 
source('R/zzz.R') 
library("mlxR")


warfa_data <- read.table("data/warfarin_data.txt", header=T)
saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")

# model1cpt<-function(psi,id,xidep) {
#   dose<-xidep[,1]
#   tim<-xidep[,2]
#   ka<-psi[id,1]
#   V<-2
#   # V<-psi[id,2]
#   k<-0.5
#   CL<-k*V
#   ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#   return(ypred)
# }


# saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",modeltype="structural",
#   psi0=matrix(c(1),ncol=1,byrow=TRUE, dimnames=list(NULL, c("ka"))),
#   transform.par=c(1),omega.init=matrix(c(1),ncol=1,byrow=TRUE),
#   covariance.model=matrix(c(1),ncol=1,byrow=TRUE))


model1cpt<-function(psi,id,xidep) {
  dose<-xidep[,1]
  tim<-xidep[,2]
  ka<-psi[id,1]
  V<-psi[id,2]
  k<-psi[id,3]
  CL<-k*V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",modeltype="structural"
  ,psi0=matrix(c(1,7,1),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V", "k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE))



##RUNS
K1 = 100
K2 = 50
iterations = 1:(K1+K2+1)
end = K1+K2

#Warfarin
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0, nbiter.map=4)
warfa<-saemix(saemix.model_warfa,saemix.data_warfa,options_warfa)
