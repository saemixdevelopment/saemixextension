workDir <- "/home/eco/work/saemix/saemixextension/lucasMultiple"
datDir <- file.path(workDir, "exemple_emmanuelle")

library(tidyverse)
library(Rcpp)

for (i in 1:2){
  fichiers = sort(list.files(path = "R" , pattern = "\\." , full.names = TRUE))
  sapply(fichiers,source)
}

fichiers_multi = list.files(path="R_multi", pattern = "\\." , full.names = TRUE)
  sapply(fichiers_multi,source)


source("gridfuncd_mean_PKPD_ver.R")
sourceCpp("saemix_modelPBPK_CPP.cpp")
source("SaemixModelPBPK_PKPD_ver_CPP.R")

######################## Model PKPD ##############################

# Modèle
modelPKPD<-function(psi,id,xidep) {
  ytype<-xidep$ytype
  dose<-xidep[,1]
  tim<-xidep[,2]
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  Emax<-psi[id,4]
  EC50<-psi[id,5]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  ypredPD<-Emax*ypred/(ypred+EC50)
  ypred[ytype==2]<-ypredPD[ytype==2]
  return(ypred)
}

# Modèle
modelPKPD.cpt<-function(psi,id,xidep) {
  ytype<-xidep$ytype
  dose<-xidep[,1]
  tim<-xidep[,2]
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  ke0 <- psi[id,4]
  Emax<-psi[id,5]
  EC50<-psi[id,6]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  ypred.site <-dose*ka*ke0/V * (exp(-ka*tim)/(k-ka)/(ke0-ka) + exp(-k*tim)/(ka-k)/(ke0-k) +
                                  exp(-ke0*tim)/(ka-ke0)/(k-ke0))
  ypredPD<-Emax*ypred.site/(ypred.site+EC50)
  ypred[ytype==2]<-ypredPD[ytype==2]
  return(ypred)
}

###########################################################
          # Simulation Numéro 1
############################################################


data_sim1 = read.csv(file.path(datDir,"basicPKPD_Sim1.csv") , header = TRUE)

# paramètres mu 
psi_ini =  c(ka = 1, V = 20, CL = 2, Emax = 100, EC50 = 5)
x.grid05 <- gridfuncd(data= data_sim1, psi0=psi_ini, e.obj = 0.05 , class = "ytype" , model = modelPKPD)


saemix.model <- saemixModelPBPK(psi0=psi_ini, 
                                x.grid=x.grid05,
                                simulations=unique(data_sim1$id), 
                                output=NULL,
				covariate.model = NULL  , 
				covariance.model=diag(1,5,5), # Matrice diagonale , pas de covariance entre les paramètres
				error.model = c("proportional" , "constant"), # proportionnel pour PK et constant pour PD
				modeltype = c("structural","structural") , 
				model.interpolate = TRUE , # on utilise la  grille
				name.response = c("y1","y2") , 
				model = modelPKPD) # le model structural 


saemix.data <- saemixData(name.data=data_sim1,header=TRUE,sep=" ",                  na=NA,name.group=c("id"),name.predictors=c("dose","tim"), name.response=c("y"), name.X="time",
					units = list(x = "mg/L", y = "hr"),name.ytype = "ytype")


options <- list(seed=54333, map=F, fim=F, ll.is=F,
                nbiter.saemix = c(300,200), nbiter.sa=20,
                nmc.is = 5000, nu.is=5, r.is=1, nbdisplay=25,
                displayProgress=T, print.is=F, save.graphs=F,
		warnings=TRUE, directory = "Fit1")

time.beg1 <- Sys.time()
fit_sim1 =  saemix.multi(saemix.model, saemix.data, options)
time.end1 <- Sys.time()
time.fit1 <- time.end1-time.beg1

saemix.model.exact <- saemixModelPBPK(psi0=psi_ini, 
                                x.grid=x.grid05,
                                simulations=unique(data_sim1$id), 
                                output=NULL,
				covariate.model = NULL  , 
				covariance.model=diag(1,5,5), # Matrice diagonale , pas de covariance entre les paramètres
				error.model = c("proportional" , "constant"), # proportionnel pour PK et constant pour PD
				modeltype = c("structural","structural") , 
				model.interpolate = FALSE , # on utilise la  grille
				name.response = c("y1","y2") , 
				model = modelPKPD) # le model structural 

options <- list(seed=54333, map=F, fim=F, ll.is=F,
                nbiter.saemix = c(300,200), nbiter.sa=20,
                nmc.is = 5000, nu.is=5, r.is=1, nbdisplay=25,
                displayProgress=T, print.is=F, save.graphs=F,
		warnings=TRUE, directory = "exactFit1")

time.beg1 <- Sys.time()
fit_sim1.exact =  saemix.multi(saemix.model.exact, saemix.data, options)
time.end1 <- Sys.time()
time.exact1 <- time.end1-time.beg1

##################### MAP


##################### FIM
