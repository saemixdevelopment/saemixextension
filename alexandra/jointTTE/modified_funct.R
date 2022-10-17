#################################################
############# SIMULATION DES DONNEES ############
#################################################

saemixDir <- "/home/eco/work/saemix/saemixextension"
workDir <- file.path(saemixDir, "alexandra","jointTTE")
setwd(workDir)

library(saemix)
library(ggplot2)
library("viridis")  

################### PKPD ########################

N = 1000
N<-200
# 4 doses possibles 
dose = sample(c(50,75,100,125),N,replace = T)

ka_pop = 2
V_pop = 8
#Cl_pop = 0.15
Cl_pop = 2

omega_ka = 1.4
omega_V = 0.15
omega_Cl = 0.3

Emax_pop = 100
EC50_pop = 5

omega_Emax = 0.15   # mettre omega_Emax et omega_EC50 = 0 pour que ça marche 
omega_EC50 = 0.2

sigma_pk = 0.1
sigma_pd = 1

data_pkpd = data.frame(id=rep(1:N,each=20),time = rep(1:10,2*N),obs=NA,ytype=rep(c(rep(1,10),rep(2,10)),N),dose=NA)

for (i in 1:N){
  ka = ka_pop*exp(rnorm(1,0,omega_ka))
  V = V_pop*exp(rnorm(1,0,omega_V))
  Cl = Cl_pop*exp(rnorm(1,0,omega_Cl))
  Emax = Emax_pop*exp(rnorm(1,0,omega_Emax))
  EC50 = EC50_pop*exp(rnorm(1,0,omega_EC50))
  
  D = dose[i]
  
  pk = function(x) (D/V)*(ka/(ka-Cl/V))*(exp(-(Cl/V)*x)-exp(-ka*x))
  pd = function(x) Emax*pk(x)/(pk(x)+EC50)
  
  data_pkpd$obs[data_pkpd$id==i & data_pkpd$ytype==1] = pk(1:10)+rnorm(10,0,sigma_pk)
  data_pkpd$obs[data_pkpd$id==i & data_pkpd$ytype==2] = pd(1:10)+rnorm(10,0,sigma_pd)
  
  data_pkpd$dose[data_pkpd$id==i] = D
}

# par(mfrow=c(1,2))
# for(i in 1:2) {
#   plot(data_pkpd$time[data_pkpd$ytype==i], data_pkpd$obs[data_pkpd$ytype==i], xlab="Time",ylab="PK/PD")
# }

ggplot(data_pkpd, aes(x=time, y=obs, group=id, colour=as.factor(dose))) + geom_line() + theme_minimal() + scale_color_viridis(discrete=TRUE) + facet_wrap(.~ytype, scales="free") + ylab("Outcome") +  guides(colour=guide_legend(title="Dose"))

#write.table(data_pkpd,file="C:/Users/AlexandraLAVALLEY/Documents/Code_saemix/essais/pkpd.txt",row.names = F)
if(FALSE)
  data_pkpd= read.table(file="C:/Users/AlexandraLAVALLEY/Documents/Code_saemix/essais/pkpd.txt",header = T)


LONGIjointmodel<-function(psi,id,xidep) {
  ytype<-xidep$ytype
  
  ka <- psi[id,1] 
  V <- psi[id,2]
  Cl <- psi[id,3]
  Emax <- psi[id,4]
  EC50 <- psi[id,5]
  
  T<-xidep[,1] 
  D <- xidep[,2]
  
  pk = function(t) (D/V)*(ka/(ka-Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))
  
  ypred <- pk(T)
  ypred[ytype==2] <- Emax*ypred[ytype==1]/(ypred[ytype==1]+EC50)
  
  return(ypred)
}


# Refait sans la fonction pk à l'intérieur
LONGIjointmodel<-function(psi,id,xidep) {
  ytype<-xidep$ytype
  
  ka <- psi[id,1] 
  V <- psi[id,2]
  Cl <- psi[id,3]
  Emax <- psi[id,4]
  EC50 <- psi[id,5]
  
  tim<-xidep[,1] 
  D <- xidep[,2]

  ypred <- (D/V)*(ka/(ka-Cl/V))*(exp(-(Cl/V)*tim)-exp(-ka*tim))
  ypred[ytype==2] <- Emax*ypred[ytype==2]/(ypred[ytype==2]+EC50)
  
  return(ypred)
}




saemix.data<-saemixData(name.data=data_pkpd, name.group=c("id"), name.predictors=c("time","dose"), name.response="obs",name.ytype = "ytype")

saemix.model<-saemixModel(model=LONGIjointmodel,description="joint pkpd",modeltype="structural",
                          psi0=matrix(c(2, 8, 0.15, 100, 5),ncol=5,byrow=TRUE,dimnames=list(NULL, c("ka","V","Cl","Emax","EC50"))),
                          transform.par=c(1,1,1,1,1), covariance.model=diag(c(1,1,1,1,1)),fixed.estim = c(1,1,1,1,1),
                          omega.init = diag(c(1.96,0.02,0.09,0.02,0.04)),error.model = c("constant","constant"),error.init = c(1,0,10,0),
                          name.sigma = c("a.1","b.1","a.2","b.2"))


saemix.model<-saemixModel(model=LONGIjointmodel,description="joint pkpd",modeltype="structural",
                          psi0=matrix(c(2, 8, 0.15, 100, 5),ncol=5,byrow=TRUE,dimnames=list(NULL, c("ka","V","Cl","Emax","EC50"))),
                          transform.par=c(1,1,1,1,1), covariance.model=diag(c(1,1,1,1,1)),fixed.estim = c(1,1,1,1,1),
                          omega.init = diag(c(1.96,0.02,0.09,0.02,0.04)),error.model = c("proportional","proportional"),error.init = c(0,1,0,1),
                          name.sigma = c("a.1","b.1","a.2","b.2"))

saemix.options<-list(seed=12345,save=FALSE,save.graphs=FALSE, fim=FALSE, displayProgress=FALSE)

########################################################################################################################################################

## LANCER LES FONCTIONS MODIFIEES SUIVANTES 
source("functions_alex.R")
source("main_2longi.R")


saemix_test(saemix.model,saemix.data,saemix.options)
########################################################################################################################################################