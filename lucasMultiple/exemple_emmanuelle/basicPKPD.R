##################################################################
# Basic PKPD model
##################################################################
# Library to simulate in a multivariate normal distribution
library(mvtnorm)

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

# Simulation de paramètres individuels
simulpar <- function(nsuj, par.mu, par.iiv) {
  psi1<-NULL
  for(i in 1:nsuj) {
    parsuj <- par.mu*exp(rnorm(length(par.mu))*par.iiv)
    psi1 <- rbind(psi1, parsuj)
    #  psi1<-matrix(parsuj,ncol=5)
  }
  return(psi1)
}

# Simulation avec mvtnorm quand matrice non diagonale
simulparcor <- function(nsuj, par.mu, omega) {
  etas <- rmvnorm(nsuj, sigma=omega)
  psi1<-NULL
  for(i in 1:length(par.mu))
    psi1 <- cbind(psi1,par.mu[i]*exp(etas[,i]))
  return(psi1)
}

# Simulation des prédictions
simuly <- function(psi1, xtim, dose) {
  nsuj<- dim(psi1)[1]
  id1 <- rep(1:nsuj,each=length(xtim)*2)
  xidep1<-data.frame(dose=dose, tim=rep(xtim,each=2), ytype=rep(c(1,2),length(xtim)))
  xidep <- do.call(rbind,rep(list(xidep1),nsuj))
  ypred <-modelPKPD(psi1, id1, xidep)
  return(ypred)
}

##################################################################
# Simulation 1

## Modele PK: 1 cpt avec absorption; paramètres ka, V, CL
## Modèle Emax: direct avec E0=0; paramètres Emax et EC50
## pas de covariable dans le modèle

## Paramètres
### paramètres PK: comme théophylline CI
### paramètres PD: Emax=100, EC50=5
#### donnant un effet maximal de 65 avec les paramètres moyens
### IIV 30% sur tous les paramètres
### pas de corrélation entre les paramètres

## Design
### nb de sujets: 10
### nb d'observations: 11/sujet
### mêmes temps de prélèvements pour PK et PD
### temps de 0 à 24h
### même dose, unique
##################################################################
# Choix de paramètres
par.mu <- c(1.,20,2, 100, 5)
par.iiv <- rep(0.3,5)
par.sig <- c(0.1, 3)

# Design - nombre de points
## mêmes points pour PK et PD
xtim <- c(0, 0.25, 0.5, 1, 2, 3, 5, 7, 9, 12, 24)
nsuj <- 10
dose <- 300
ytype <- rep(rep(c(1,2),length(xtim)), nsuj)

## Simulation
psi1<- simulpar(nsuj, par.mu, par.iiv)
ypred <- simuly(psi1, xtim, dose)
nobs<-length(ypred)

## Ajout variabilité résiduelle
gpred <- ypred.nosig <- ypred
gpred[ytype==1] <- par.sig[1]*abs(ypred[ytype==1])
gpred[ytype==2] <- par.sig[2]
ypred <- ypred + gpred*rnorm(nobs)

ypl <- data.frame(id=rep(1:nsuj,each=length(xtim)*2),dose=rep(dose,nsuj*2),tim=rep(rep(xtim,each=2),nsuj),y=ypred, ytype=ytype)

# Graphe des données simulées pour les 10 sujets
par(mfrow=c(1,2))
plot(ypl$tim[ypl$ytype==1],ypl$y[ypl$ytype==1],xlab="Time after dose (hr)", ylab="Theophylline concentration (mg/L)")
plot(ypl$tim[ypl$ytype==2],ypl$y[ypl$ytype==2],xlab="Time after dose (hr)", ylab="Effet (%)")

# Graphe des données prédites pour le premier sujet - sans erreur résiduelle (prédiction pure)
par(mfrow=c(1,1))
plot(ypl$tim[ypl$id==1], ypred.nosig[ypl$id==1],xlab="Time after dose (hr)", ylab="Concentrations (black) and Responses (blue)")
for(i in 1:2) lines(ypl$tim[ypl$ytype==i & ypl$id==1],ypred.nosig[ypl$ytype==i & ypl$id==1], col=ifelse(i==1,"black","blue"))

# Sauver les données en format saemix
write.csv(ypl, "basicPKPD_Sim1.csv", row.names=F, quote=F)

##################################################################
# Simulation 2
# Plus de sujets (devrait stabiliser les estimations)
##################################################################
### nb de sujets: 50 (peut-être nécessaire de monter à 100 ou 200)
##################################################################
nsuj<-50

## Simulation
psi1<- simulpar(nsuj, par.mu, par.iiv)
ypred <- simuly(psi1, xtim, dose)
nobs<-length(ypred)
ytype <- rep(rep(c(1,2),length(xtim)), nsuj)

## Ajout variabilité résiduelle
gpred <- ypred.nosig <- ypred
gpred[ytype==1] <- par.sig[1]*abs(ypred[ytype==1])
gpred[ytype==2] <- par.sig[2]
ypred <- ypred + gpred*rnorm(nobs)

ypl <- data.frame(id=rep(1:nsuj,each=length(xtim)*2),dose=rep(dose,nsuj*2),tim=rep(rep(xtim,each=2),nsuj),y=ypred, ytype=ytype)

# Sauver les données en format saemix
write.csv(ypl, "basicPKPD_Sim2.csv", row.names=F, quote=F)

##################################################################
# Simulation 3
# Modèle d'erreur résiduelle: combiné pour PK, proportionnel pour PD
##################################################################
nsuj<-50
par.sig <- c(0.5, 0.1, 0.2)

## Simulation
psi1<- simulpar(nsuj, par.mu, par.iiv)
ypred <- simuly(psi1, xtim, dose)
nobs<-length(ypred)
ytype <- rep(rep(c(1,2),length(xtim)), nsuj)

## Ajout variabilité résiduelle
gpred <- ypred.nosig <- ypred
gpred[ytype==1] <- sqrt(par.sig[1]**2 + (par.sig[2]*ypred[ytype==1])**2)
gpred[ytype==2] <- par.sig[2]*ypred[ytype==2]
ypred <- ypred + gpred*rnorm(nobs)

ypl <- data.frame(id=rep(1:nsuj,each=length(xtim)*2),dose=rep(dose,nsuj*2),tim=rep(rep(xtim,each=2),nsuj),y=ypred, ytype=ytype)

# Sauver les données en format saemix
ypl <- cbind(ypl, dose=dose)
write.csv(ypl, "basicPKPD_Sim3.csv", row.names=F, quote=F)

##################################################################
# Simulation 4
# Corrélation entre deux paramètres (CL et EC50)
# objectif: tester que ça marche avec une matrice de var-cov non diagonale
##################################################################
nsuj<-50
par.sig <- c(0.1, 3)
ytype <- rep(rep(c(1,2),length(xtim)), nsuj)

omega <- diag(par.iiv**2, ncol=length(par.mu))
omega[2,3]<-omega[3,2]<-0.7*sqrt(omega[2,2])*sqrt(omega[3,3])
psi1<- simulparcor(nsuj, par.mu, omega)
# check (on doit retomber à peu près sur par.iiv)
# apply(psi1,2,sd)/par.mu

ypred <- simuly(psi1, xtim, dose)
nobs<-length(ypred)

## Ajout variabilité résiduelle
gpred <- ypred.nosig <- ypred
gpred[ytype==1] <- par.sig[1]*abs(ypred[ytype==1])
gpred[ytype==2] <- par.sig[2]
ypred <- ypred + gpred*rnorm(nobs)

ypl <- data.frame(id=rep(1:nsuj,each=length(xtim)*2),dose=rep(dose,nsuj*2),tim=rep(rep(xtim,each=2),nsuj),y=ypred, ytype=ytype)

# Sauver les données en format saemix
ypl <- cbind(ypl, dose=dose)
write.csv(ypl, "basicPKPD_Sim4.csv", row.names=F, quote=F)

##################################################################
# Simulation 5
# Pas d'IIV sur ka
# Corrélation entre deux paramètres (CL et EC50)
# objectif: tester que ça marche avec une matrice de var-cov non diagonale et non complète
##################################################################
nsuj<-50
par.mu <- c(1.,20,2, 100, 5)
par.iiv <- c(0,rep(0.3,4))
par.sig <- c(0.1, 3)
ytype <- rep(rep(c(1,2),length(xtim)), nsuj)

## Simulation
omega <- diag(par.iiv[2:5]**2, ncol=(length(par.mu)-1))
omega[2,1]<-omega[1,2]<-0.7*sqrt(omega[2,2])*sqrt(omega[1,1])
psi1<- simulparcor(nsuj, par.mu[2:5], omega)
psi1 <- cbind(rep(par.mu[1],nsuj),psi1)
ypred <- simuly(psi1, xtim, dose)
nobs<-length(ypred)

## Ajout variabilité résiduelle
gpred <- ypred.nosig <- ypred
gpred[ytype==1] <- par.sig[1]*abs(ypred[ytype==1])
gpred[ytype==2] <- par.sig[2]
ypred <- ypred + gpred*rnorm(nobs)

ypl <- data.frame(id=rep(1:nsuj,each=length(xtim)*2),dose=rep(dose,nsuj*2),tim=rep(rep(xtim,each=2),nsuj),y=ypred, ytype=ytype)

# Sauver les données en format saemix
ypl <- cbind(ypl, dose=dose)
write.csv(ypl, "basicPKPD_Sim5.csv", row.names=F, quote=F)

##################################################################
# Simulation 6
# 3 doses différentes
# temps de prélèvements différents pour PK et PD
##################################################################
nsuj<-50
par.mu <- c(1.,20,2, 100, 5)
par.iiv <- rep(0.3,5)
par.sig <- c(0.1, 3)
xtim.PK <- c(0, 0.25, 0.5, 1, 2, 3, 5, 7, 9, 12, 24)
xtim.PD <- c(0, 0.5, 1, 2, 4, 8, 16, 24, 36, 48)
xtim <- c(xtim.PK, xtim.PD)
dose <- c(100,300,1000)

psi1<- simulpar(nsuj, par.mu, par.iiv)
id1 <- rep(1:nsuj,each=length(xtim))
xidep1<-data.frame(dose=rep(dose[2],length(xtim)), tim=xtim, ytype=c(rep(1,length(xtim.PK)), rep(2,length(xtim.PD))))
xidep1 <- xidep1[order(xidep1$tim, xidep1$ytype),]
xidep <- do.call(rbind,rep(list(xidep1),nsuj))
xidep$dose[id1<16] <- dose[1]
xidep$dose[id1>35] <- dose[3]
ypred <-modelPKPD(psi1, id1, xidep)
nobs<-length(ypred)
ytype <- xidep$ytype

## Ajout variabilité résiduelle
gpred <- ypred.nosig <- ypred
gpred[ytype==1] <- par.sig[1]*abs(ypred[ytype==1])
gpred[ytype==2] <- par.sig[2]
ypred <- ypred + gpred*rnorm(nobs)

ypl <- data.frame(id=id1,xidep[,2:3],y=ypred)

# Sauver les données en format saemix
ypl <- cbind(ypl, dose=dose)
write.csv(ypl, "basicPKPD_Sim6.csv", row.names=F, quote=F)

par(mfrow=c(1,1))
plot(ypl$tim[id1 %in% c(1,20,40)], ypred.nosig[id1 %in% c(1,20,40)],xlab="Time after dose (hr)", ylab="Concentrations (black) and Responses (blue)")
for(i in 1:2) lines(ypl$tim[ypl$ytype==i & id1 %in% c(1,20,40)],ypred.nosig[ypl$ytype==i & id1 %in% c(1,20,40)], col=ifelse(i==1,"black","blue"))

##################################################################
# Autres simulations possibles
## Paramètres
### covariables dans le modèle (eg poids sur V et CL)
### varier l'IIV (50 ou 70%)
### IIV différente selon les paramètres (IIV de 30% sur V, CL, Emax et 70% sur ka et EC50)
### varier ratio entre IIV et résiduelle (ie varier sigma pour PK et PD)

## Design
### nb de sujets: varier (N=10 [peu] ou N=50/100 [bcp])
### nb d'observations: varier (n=4), possiblement descendre avec des prélèvements différents selon les sujets (eg: 3 groupes de sujets avec 2 prélèvements PK et 2 prélèvements PD chacun, à des temps répartis intelligemment, possiblement avec des doses différentes)

##################################################################

