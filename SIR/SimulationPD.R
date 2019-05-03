##################### Simulation PD Emax DOSNE 2016 ######################



E0 <- 10
Emax <- 100
ED50 <- 5 
IIVE0 <- 0.3
IIVED50 <- 0.3
RUVprop <- 0.1

modelPDEmax <- function(E0, Emax, ED50, D){
  return((Emax*D/(ED50+D)+E0))
}

#Fonction de simulation du modèle PD Emax
simPDEmax <- function(nsubj, doses){
  ndose<-length(doses)
  subj<-data.frame(id=1:nsubj, E0=E0*exp(rnorm(nsubj,0,IIVE0)),ED50=ED50*exp(rnorm(nsubj,0,IIVED50)),Emax=Emax)
  datasim<-do.call(rbind,rep(list(subj),ndose))
  datasim<-datasim[order(datasim$id),]
  datasim$dose<-rep(doses,nsubj)
  datasim$ypred<-modelPDEmax(datasim$E0, datasim$Emax, datasim$ED50, datasim$dose)
  datasim$ysim<-datasim$ypred*(1+rnorm(nsubj*ndose,0,RUVprop))
  par(mfrow=c(1,1))
  plot(datasim$dose,datasim$ysim,pch=20, xlab="Doses", ylab="Simulated Data")
  return(datasim)
}

#3 simulations a priori réalisées dans l'article de Dosne 2016 (à voir pour le nombre de doses)
set.seed(1234)
simPDEmax20 <- simPDEmax(20, c(0, 2.5, 5, 15))
View(simPDEmax20)
set.seed(2345)
simPDEmax50 <- simPDEmax(50, c(0, 2.5, 5, 15))
set.seed(3456)
simPDEmax200 <- simPDEmax(200, c(0, 2.5, 5, 15))


#############################################################
#############################################################
#####################      SAEMIX       #####################
#############################################################
#############################################################


#############################################################
################# PD Emax 20 ID - 4 obs/ID ##################
#############################################################

#Construction du modèle PD Emax
modelPD<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  E0<-psi[id,1]
  Emax<-psi[id,2]
  ED50<-psi[id,3]
  ypred<-E0+Emax*dose/(ED50+dose)
  return(ypred)
}

#Construction des données pour simulation 20 sujets
saemix.data.pd20<-saemixData(name.data=simPDEmax20, 
                        name.group=c("id"),name.predictors=c("dose"),
                        name.response=c("ysim"), name.X="dose") 

#Graphes des données
plot(saemix.data.pd20,type="b",col="DarkBlue",main="Simulated dataset 1 (N=20)", xlab='dose (mg)', ylab='response (-)')
#Graphes des données individuelles
plot(saemix.data.pd20,individual=TRUE)

#psi contient les parametres (3 parametres ici)
#xidep contient les predicteurs (ici on en a 1)
#Attention aux ordres
#creation du modele type saemix
#Default model, no covariate
saemix.model.pd20<-saemixModel(model=modelPD,
                          description="PD Emax model",
                          psi0=matrix(c(20,200,10),ncol=3,byrow=TRUE, 
                          dimnames=list(NULL, c("E0","Emax","ED50"))), fixed.estim = c(1,1,1),transform.par=c(1,1,1), 
                          error.model = "proportional", covariance.model = matrix(c(1,0,0,0,0,0,0,0,1),nrow=3, ncol=3, byrow=T))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
sim.pd20.saemix<-saemix(saemix.model.pd20,saemix.data.pd20) #dans saemix : model, data, options eventuellement


#recupere que les resultats (sans modele, sans data)
y<-summary(sim.pd20.saemix)
print(y$fixed.effects)
print(y$logLik$X.2xLL[2])
sim.pd20.saemix["results"]["ll.is"]

print(head(y$coefficients$random$map.psi))



# Graphes
plot(sim.pd20.saemix)
plot(sim.pd20.saemix,plot.type="observations.vs.predictions")
plot(sim.pd20.saemix,plot.type="both")
plot(sim.pd20.saemix,plot.type="vpc")
plot(sim.pd20.saemix,plot.type="residuals.scatter")
plot(sim.pd20.saemix,plot.type="residuals.scatter",which.pres="npde")
plot(sim.pd20.saemix,plot.type="residuals.distribution",which.pres="npde")
par(mfrow=c(3,2))
cat("Distribution des parametres individuels estimes\n")
plot(saemix.fit,plot.type="marginal",indiv.histo=TRUE)



#############################################################
################# PD Emax 50 ID - 4 obs/ID ##################
#############################################################
saemix.data50<-saemixData(name.data=sim.pd50, #header=TRUE,sep=" ",na=NA,
                        name.group=c("ID"),name.predictors=c("dose"),
                        name.response=c("E"),#name.covariates=c("Weight","Sex"),
                        name.X="Dose") #unites pour les graphes, pareil pour name.X, donnee en abscisse
plot(saemix.data50,type="b",col="DarkBlue",main="Simulated dataset 2 (N=50)", xlab='dose (mg)', ylab='response (-)')
plot(saemix.data50,individual=TRUE)

#psi contient les parametres (3 parametres ici)
#xidep contient les predicteurs (ici on en a 1)
#Attention aux ordres
modelPD<-function(psi,id,xidep) {
  dose<-xidep[,1]
  E0<-psi[id,1]
  Emax<-psi[id,2]
  ED50<-psi[id,3]
  ypred<-E0+Emax*dose/(ED50+dose)
  return(ypred)
}

#creation du modele type saemix
# Default model, no covariate
saemix.model50<-saemixModel(model=modelPD,
                            description="PD Emax model",
                            psi0=matrix(c(20,200,10),ncol=3,byrow=TRUE, 
                                        dimnames=list(NULL, c("E0","Emax","ED50"))), fixed.estim = c(1,1,1),transform.par=c(1,1,1), 
                            error.model = "proportional", covariance.model = matrix(c(1,0,0,0,0,0,0,0,1),nrow=3, ncol=3, byrow=T))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
saemix.fit50<-saemix(saemix.model50,saemix.data50,saemix.options) #dans saemix : model, data, options eventuellement
sim.pd50.saemix<-saemix.fit50
#recupere que les resultats (sans modele, sans data)
y50<-summary(sim.pd50.saemix)
print(y50$fixed.effects)

print(y50$logLik$X.2xLL[2])

sim.pd50.saemix["results"]["ll.is"]

cat("Parametres individuels (debut du tableau des estimations MAP)\n")
print(head(y50$coefficients$random$map.psi))

# Graphes
plot(sim.pd50.saemix)

plot(sim.pd50.saemix,plot.type="both")
plot(sim.pd50.saemix,plot.type="vpc")
plot(sim.pd50.saemix,plot.type="residuals.scatter")
plot(sim.pd50.saemix,plot.type="residuals.scatter",which.pres="npde")
#regarder les graphes de droites
plot(sim.pd50.saemix,plot.type="residuals.distribution",which.pres="npde")
par(mfrow=c(3,2))
cat("Distribution des parametres individuels estimes\n")
plot(saemix.fit50,plot.type="marginal",indiv.histo=TRUE)







#############################################################
################# PD Emax 200 ID - 4 obs/ID ##################
#############################################################
saemix.data200<-saemixData(name.data=sim.pd200, #header=TRUE,sep=" ",na=NA,
                        name.group=c("ID"),name.predictors=c("dose"),
                        name.response=c("E"),#name.covariates=c("Weight","Sex"),
                        name.X="Dose") #unites pour les graphes, pareil pour name.X, donnee en abscisse
plot(saemix.data200,type="b",col="DarkBlue",main="Simulated dataset 3 (N=200)", xlab='dose (mg)', ylab='response (-)')

#psi contient les parametres (3 parametres ici)
#xidep contient les predicteurs (ici on en a 1)
#Attention aux ordres
modelPD<-function(psi,id,xidep) {
  dose<-xidep[,1]
  E0<-psi[id,1]
  Emax<-psi[id,2]
  ED50<-psi[id,3]
  ypred<-E0+Emax*dose/(ED50+dose)
  return(ypred)
}

#creation du modele type saemix
# Default model, no covariate
saemix.model200<-saemixModel(model=modelPD,
                             description="PD Emax model",
                             psi0=matrix(c(20,200,10),ncol=3,byrow=TRUE, 
                                         dimnames=list(NULL, c("E0","Emax","ED50"))), fixed.estim = c(1,1,1),transform.par=c(1,1,1), 
                             error.model = "proportional", covariance.model = matrix(c(1,0,0,0,0,0,0,0,1),nrow=3, ncol=3, byrow=T))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
saemix.fit200<-saemix(saemix.model200,saemix.data200,saemix.options) #dans saemix : model, data, options eventuellement
sim.pd200.saemix<-saemix.fit200
summary(sim.pd200.saemix) #recupere que les resultats (sans modele, sans data)
y200<-summary(sim.pd200.saemix)
print(y200$fixed.effects)

print(y200$logLik$X.2xLL[2])

sim.pd200.saemix["results"]["ll.is"]
cat("Parametres individuels (debut du tableau des estimations MAP)\n")

print(head(y200$coefficients$random$map.psi))

# Graphes
plot(sim.pd200.saemix)
cat("Tracer des graphes particuliers (voir argument type)\n")
plot(sim.pd200.saemix,plot.type="both")
plot(sim.pd200.saemix,plot.type="vpc")
plot(sim.pd200.saemix,plot.type="residuals.scatter")
plot(sim.pd200.saemix,plot.type="residuals.scatter",which.pres="npde")
#regarder les graphes de droites
plot(sim.pd200.saemix,plot.type="residuals.distribution",which.pres="npde")
par(mfrow=c(3,2))
cat("Distribution des parametres individuels estimes\n")
plot(saemix.fit200,plot.type="marginal",indiv.histo=TRUE)


######################################################################
###################         TEST THEO.SAEMIX        ##################
######################################################################

data(theo.saemix)

#SANS COVARIABLE
saemix.data.theo<-saemixData(name.data=theo.saemix, #header=TRUE,sep=" ",na=NA, 
                        name.group=c("Id"),name.predictors=c("Dose","Time"),
                        name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
                        units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time") #unites pour les graphes, pareil pour name.X, donnee en abscisse
model.theo<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}
saemix.model.theo<-saemixModel(model=model.theo,
                          description="One-compartment model with first-order absorption",
                          psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, 
                                      dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))
saemix.options.theo<-list(seed=632545,save=FALSE,save.graphs=FALSE)

saemix.fit.theo<-saemix(saemix.model.theo,saemix.data.theo,saemix.options.theo)

SIRTheo1 <- saemixSIR(saemix.fit.theo, M=20, m=4)
SIRTheo1@name.param
SIRTheo1@est.mu
SIRTheo1@inflation
SIRTheo1@cov.mat
SIRTheo1@optionll
SIRTheo1@IR
SIRTheo1@resampled.theta
diag(SIRTheo1@cov.mat)
SIRTheo1@sdSIR


#AVEC COVARIABLE

saemix.covmodel.theo<-saemixModel(model=model.theo,
                             description="One-compartment model with first-order absorption",
                             psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
                             transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
                             omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE)) #,error.model="combined")
saemix.cov.theo<-saemix(saemix.covmodel.theo,saemix.data.theo,saemix.options.theo)

SIRTheo1cov <- saemixSIR(saemix.cov.theo, M=20, m=4)
SIRTheo1cov@name.param
SIRTheo1cov@est.mu
SIRTheo1cov@inflation
SIRTheo1cov@cov.mat
SIRTheo1cov@optionll
SIRTheo1cov@IR
SIRTheo1cov@resampled.theta
diag(SIRTheo1cov@cov.mat)
SIRTheo1cov@sdSIR


# AVEC COVARIABLE ET COVARIANCE
saemix.iivmodel.theo<-saemixModel(model=model.theo,
                             description="One-compartment model with first-order absorption",
                             psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
                             transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
                             covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),error.model="combined")
saemix.iiv.theo<-saemix(saemix.iivmodel.theo,saemix.data.theo,saemix.options.theo)

SIRTheo2 <- saemixSIR(saemix.iiv.theo, M=20, m=4)
SIRTheo2@name.param
SIRTheo2@est.mu
SIRTheo2@inflation
SIRTheo2@cov.mat
SIRTheo2@optionll
SIRTheo2@IR
SIRTheo2@resampled.theta
diag(SIRTheo2@cov.mat)
SIRTheo2@sdSIR

