##################################################################
# PKPD model with effect compartment
##################################################################

# Library to simulate in a multivariate normal distribution
library(mvtnorm)

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
## Compartiment d'effet: paramètre ke0
## Modèle Emax: direct avec E0=0; paramètres Emax et EC50
## pas de covariable dans le modèle

## Paramètres
### paramètres PK: comme théophylline CI
### ke0 = 0.13 (1/2 vie d'équilibre 5h)
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
npar<-6
par.mu <- c(1.,20,2, log(2)/5, 100, 5)
par.iiv <- rep(0.3,npar)
par.sig <- c(0.1, 3)
dose <- 300

# Design - nombre de points
## mêmes points pour PK et PD
xtim <- c(0, 0.25, 0.5, 1, 2, 3, 5, 7, 9, 12, 24)
ytype <- rep(rep(c(1,2),length(xtim)), nsuj)

psi1<-NULL
nsuj<-50
for(i in 1:nsuj) {
  parsuj <- par.mu*exp(rnorm(npar)*par.iiv)
  psi1 <- rbind(psi1, parsuj)
  #  psi1<-matrix(parsuj,ncol=5)
}
id1 <- rep(1:nsuj,each=length(xtim)*2)
xidep1<-data.frame(dose=300, tim=rep(xtim,each=2), ytype=rep(c(1,2),length(xtim)))
xidep <- do.call(rbind,rep(list(xidep1),nsuj))
ypred <-modelPKPD.cpt(psi1, id1, xidep)

## Ajout variabilité résiduelle
nobs<-length(ypred)
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
write.csv(ypl, "effectCptPKPD_Sim1.csv", row.names=F, quote=F)

##################################################################
# Simulation 2
# Emax 1000 fois plus grand
# objectif: vérifier que l'échelle des données n'a pas d'importance
##################################################################
nsuj<-50
par.mu <- c(1.,20,2, log(2)/5, 100000, 5)
par.sig <- c(0.1, 3000) # en additif il faut aussi changer la résiduelle

psi1<-NULL
nsuj<-10
for(i in 1:nsuj) {
  parsuj <- par.mu*exp(rnorm(npar)*par.iiv)
  psi1 <- rbind(psi1, parsuj)
  #  psi1<-matrix(parsuj,ncol=5)
}
id1 <- rep(1:nsuj,each=length(xtim)*2)
xidep1<-data.frame(dose=300, tim=rep(xtim,each=2), ytype=rep(c(1,2),length(xtim)))
xidep <- do.call(rbind,rep(list(xidep1),nsuj))
ypred <-modelPKPD.cpt(psi1, id1, xidep)

## Ajout variabilité résiduelle
nobs<-length(ypred)
gpred <- ypred.nosig <- ypred
gpred[ytype==1] <- par.sig[1]*abs(ypred[ytype==1])
gpred[ytype==2] <- par.sig[2]
ypred <- ypred + gpred*rnorm(nobs)

ypl <- data.frame(id=rep(1:nsuj,each=length(xtim)*2),dose=rep(dose,nsuj*2),tim=rep(rep(xtim,each=2),nsuj),y=ypred, ytype=ytype)

# Sauver les données en format saemix
write.csv(ypl, "effectCptPKPD_Sim2.csv", row.names=F, quote=F)

##################################################################
# Autres simulations possibles

# Mêmes idées que pour l'autre modèle
## structure d'IIV (certains paramètres sans IIV)
## structure de la matrice de variance-covariance (avec ou sans corrélations)
## erreur résiduelle
## différents temps de prélèvements pour la PK et la PD

## Paramètres
### covariables dans le modèle (eg poids sur V et CL)
### varier l'IIV (50 ou 70%)
### IIV différente selon les paramètres (IIV de 30% sur V, CL, Emax et 70% sur ka et EC50)
### varier ratio entre IIV et résiduelle (ie varier sigma pour PK et PD)

## Design
### nb de sujets: varier (N=10 [peu] ou N=50/100 [bcp])
### nb d'observations: varier (n=4), possiblement descendre avec des prélèvements différents selon les sujets (eg: 3 groupes de sujets avec 2 prélèvements PK et 2 prélèvements PD chacun, à des temps répartis intelligemment, possiblement avec des doses différentes)
##################################################################
