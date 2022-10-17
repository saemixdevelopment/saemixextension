#################################################
saemixDir <- "/home/eco/work/saemix/saemixextension"
workDir <- file.path(saemixDir, "alexandra","ecoJoint")
setwd(workDir)

library(ggplot2)
library(Cairo)
library("viridis")  

# Chargement des fonctions originelles de la librairie
progDir<-file.path(saemixDir, "R")
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"func_plots.R")) # for saemix.plot.setoptions

# model
pkTTE<-function(psi,id,xidep) {
  ka <- psi[id,1] 
  V <- psi[id,2]
  Cl <- psi[id,3]
  
  tim<-xidep[,1] 
  D <- xidep[,2]
  ytype<-xidep$ytype
  
  T<-xidep[ytype==2,1]
  y<-xidep[ytype==2,2] # events (1=event, 0=no event)
  cens <- 30 # common censoring time
#  cens<-which(xidep[ytype==2,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  lambda <- psi[id,4]
  lambda <- lambda[ytype==2]
  beta <- psi[id,5]
  beta <-beta[ytype==2]
  Nj <- length(T)
  
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  hazard <- (beta/lambda)*(T/lambda)^(beta-1) # ln(H')
  H <- (T/lambda)^beta # ln(H)
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))

  ypred <- (D/V)*(ka/(ka-Cl/V))*(exp(-(Cl/V)*tim)-exp(-ka*tim))
  ypred[ytype==2] <- logpdf
  #  ypd <- Emax*ypred/(ypred+EC50)
  #  ypred[ytype==2] <- ypd[ytype==2]
  return(ypred)
}


# Proportional error model
param<-c(2,8,2,10,2)
pkTTE.prop<-saemixModel(model=pkTTE,description="Mock PK+TTE (separate)",modeltype=c("structural","likelihood"),
                            psi0=matrix(param,ncol=5,byrow=TRUE,dimnames=list(NULL, c("ka","V","Cl","Emax","EC50"))),
                            transform.par=c(1,1,1,1,1), covariance.model=diag(c(1,1,1,1,1)),fixed.estim = c(1,1,1,1,1),
                            omega.init = diag(rep(0.5,5)),error.model = "proportional",error.init = c(0,1),
                            name.sigma = c("a","b"))
pkTTE.prop@modeltype
pkTTE.prop@error.model
pkTTE.prop@error.init
pkTTE.prop@name.sigma
