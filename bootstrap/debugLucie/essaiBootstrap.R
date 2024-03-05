library(saemix)
workDir <- "/home/eco/work/saemix/saemixextension/bootstrap/debugLucie"
setwd(workDir)
nboot<-1000

# Initial parameters
ka0 = 1.5
Cl0 = 0.04
V0 = 0.5

corr <- matrix(c(1,0.90,0.98,
                 0.9,1,0.9,
                 0.98,0.9,1),nrow=3,ncol=3)
omega <- c(sqrt(1.2),sqrt(1.2),sqrt(1.2))
BSV <- omega * t(corr * omega)

# PK model ----
model1cpt_firstorder=function(psi,id,xidep) { 
  dose=4
  tim=xidep[,1]  
  ka=psi[id,1]
  CL=psi[id,2]
  V = psi[id,3]
  k = CL/V
  
  pred = dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(pred)
}

for(i in 9:10) {
  simu = read.table(paste0("dataset",i,".txt"),header=T)
  simu = simu[simu$Time>0,]
  saemix_simu=saemixData(name.data=simu,
                         name.group=c("Id"),
                         name.predictors=c("Time"),
                         name.response=c("Concentration"),
                         name.covariate="Tr",
                         name.X="Time",
                         units=list(x="day",y="mg/mL"))
  
  saemix_model=saemixModel(model=model1cpt_firstorder,
                           psi0=matrix(c(ka0,Cl0,V0),ncol=3, 
                                       dimnames=list(NULL, c("ka","CL","V"))),
                           transform.par=c(1,1,1),
                           covariance.model=matrix(c(1,1,1,
                                                     1,1,1,
                                                     1,1,1),ncol=3,byrow=TRUE),
                           omega.init=matrix(as.vector(BSV),ncol=3,byrow=TRUE),
                           covariate.model=c(1,1,1),
                           error.model="proportional")
  
  saemix_options=list(save=FALSE,save.graphs=FALSE,nb.chains = 10, nbiter.saemix=c(300,100))
  saemix_fit=saemix(saemix_model,saemix_simu,saemix_options)
  if(i==9) fit9<-saemix_fit else fit10<-saemix_fit
  
  boot1 <- saemix.bootstrap(saemix_fit, method="case", nboot=nboot)
  if(i==9) boot9<-boot1 else boot10<-boot1
  
  write.table(boot1, paste0("ecoboot",i,".txt"), row.names=FALSE, quote=F)
}

if(FALSE){
  boot9<-read.table(paste0("ecoboot",9,".txt"), header=T)
  boot10<-read.table(paste0("ecoboot",10,".txt"), header=T)
  summary(boot9)
  summary(boot10)
}