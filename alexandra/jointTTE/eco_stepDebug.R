#################################################
saemixDir <- "/home/eco/work/saemix/saemixextension"
workDir <- file.path(saemixDir, "alexandra","jointTTE")
setwd(workDir)

# library(saemix)
library(ggplot2)
library(Cairo)
library("viridis")  

simPKPD <- FALSE

# Chargement des fonctions originelles de la librairie
progDir<-file.path(saemixDir, "R")
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"main.R"))
source(file.path(progDir,"func_aux.R"))
source(file.path(progDir,"main_initialiseMainAlgo.R"))
source(file.path(progDir,"main_estep.R"))
source(file.path(progDir,"main_mstep.R"))
source(file.path(progDir,"func_FIM.R"))
source(file.path(progDir,"func_plots.R"))
source(file.path(progDir,"func_distcond.R"))
source(file.path(progDir,"func_simulations.R"))
source(file.path(progDir,"compute_LL.R"))
source(file.path(progDir,"func_estimParam.R"))
source(file.path(progDir,"func_npde.R"))
source(file.path(progDir,"backward.R"))
source(file.path(progDir,"forward.R"))
source(file.path(progDir,"stepwise.R"))
source(file.path(progDir,"func_stepwise.R"))
source(file.path(progDir,"func_compare.R"))
source(file.path(progDir,"func_bootstrap.R"))
source(file.path(progDir,"func_exploreData.R"))
source(file.path(progDir,"func_discreteVPC.R"))

#################################################
# Replace main with the beginning of the algorithm
source("copy_estep_print.R")

#################################################
# Reading data
data_pkpd.prop <- read.csv(file.path(workDir, "simAlex_pkpd_proportional.csv"), header=TRUE)

#################################################
# Fitting PK alone
# PK model
modelPK<-function(psi,id,xidep) {
  ka <- psi[id,1] 
  V <- psi[id,2]
  Cl <- psi[id,3]
  tim<-xidep[,1] 
  D <- xidep[,2]
  ypred <- (D/V)*(ka/(ka-Cl/V))*(exp(-(Cl/V)*tim)-exp(-ka*tim))
  return(ypred)
}

data_pkpd <- data_pkpd.prop[data_pkpd.prop$ytype==1,]
saemix.dataPK.prop<-saemixData(name.data=data_pkpd, name.group=c("id"), name.predictors=c("time","dose"), name.response="obs")

pkmodel.prop<-saemixModel(model=modelPK,description="Only PK",modeltype="structural",
                          psi0=matrix(c(2, 8, 0.15),ncol=3,byrow=TRUE,dimnames=list(NULL, c("ka","V","Cl"))),
                          transform.par=c(1,1,1), covariance.model=diag(c(1,1,1)),
                          omega.init = diag(c(0.5,0.5,0.5)),error.model = c("proportional"),error.init = c(0,1))
#                         omega.init = diag(c(1.96,0.02,0.09)),error.model = c("proportional"),error.init = c(0,1))
saemix.options.prop<-list(seed=12345,save=FALSE,save.graphs=FALSE, fim=FALSE, displayProgress=FALSE, nbiter.saemix=c(150,100))

model <- pkmodel.prop
data <- saemix.dataPK.prop
control <- saemix.options.prop

# Without computing delta on eta in the first kernel
source("copy_estep_print.R")
source("libmain_fordebug.R")

postscript(file="convplot_pkseul.eps",horizontal=T)
convplot.infit(allpar,saemix.options$nbiter.saemix[1],niter=(kiter-2))
dev.off()

# Computing delta on eta in the first kernel
if(FALSE) {
  source("copy_estep_delta.R")
  source("libmain_fordebug.R")
  # unchanged
}


######################################################
# Keeping track of AR
source("copy_estep_checkAR.R")
source("libmain_fordebug.R")

# Acceptance ratios
nbeta<-dim(keepAR)[2]-3
colnames(keepAR)<-c("kiter","kernel","uiter",model@name.modpar)
keepAR <- as.data.frame(keepAR)
gtab<-cbind(keepAR[keepAR$kernel==1, 1:4], param="Joint kernel 1")
colnames(gtab)[4]<-"AR"
gtab2<-NULL
for(ipar in 1:nbeta) {
  gtab1<-cbind(keepAR[keepAR$kernel>1 & keepAR$uiter==3, c(1:3,3+ipar)], param=model@name.modpar[ipar])
  xtab <- cbind(keepAR[keepAR$kernel==23 & keepAR$uiter==1, c(1:3,3+ipar)], param=model@name.modpar[ipar])
  colnames(gtab1)[4]<-colnames(xtab)[4]<-"AR"
  gtab<-rbind(gtab, gtab1)
  gtab2<-rbind(gtab2,xtab)
}
gtab<-gtab[!is.na(gtab$AR),]
gpl <- ggplot(gtab, aes(x=kiter, y=AR, group=param, colour=as.factor(param))) + geom_line() + geom_line(data=gtab2, aes(x=kiter, y=AR, group=param, colour=as.factor(param)), linetype="dashed") + ylab("Acceptance ratio") + theme_minimal() + scale_color_viridis(discrete=TRUE) + guides(colour=guide_legend(title="Parameter/Step"))


cairo_ps(file = "acceptanceRates_pkseul250.eps", onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
print(gpl)
dev.off()

#################################################
# PD alone

# PD with individual parameters for PK
modelPD<-function(psi,id,xidep) {
  Emax <- psi[id,1]
  EC50 <- psi[id,2]
  tim<-xidep[,1] 
  D <- xidep[,2]
  ka <- xidep[,3]
  V <- xidep[,4]
  Cl <- xidep[,5]
  ypred <- (D/V)*(ka/(ka-Cl/V))*(exp(-(Cl/V)*tim)-exp(-ka*tim))
  ypred <- Emax*ypred/(ypred+EC50)
  return(ypred)
}

data_pd <- data_pkpd.prop[data_pkpd.prop$ytype==2,]
saemix.dataPD.prop<-saemixData(name.data=data_pd, name.group=c("id"), name.predictors=c("time","dose", "ka","V","CL"), name.response="obs")

pdmodel.prop<-saemixModel(model=modelPD,description="Only PD",modeltype="structural",
                          psi0=matrix(c(100, 5),ncol=2,byrow=TRUE,dimnames=list(NULL, c("Emax","EC50"))),
                          transform.par=c(1,1), covariance.model=diag(c(1,1)),
                          omega.init = diag(c(0.5,0.5)),error.model = c("proportional"),error.init = c(0,1))


model <- pdmodel.prop
data <- saemix.dataPD.prop
control <- saemix.options.prop

# Without computing delta on eta in the first kernel
source("copy_estep_print.R")
source("libmain_fordebug.R")

# Computing delta on eta in the first kernel
if(FALSE) {
  source("copy_estep_delta.R")
  source("libmain_fordebug.R")
}
postscript(file="convplot_pdseul.eps",horizontal=T)
convplot.infit(allpar,saemix.options$nbiter.saemix[1],niter=(kiter-2))
dev.off()

######################################################
# Keeping track of AR
source("copy_estep_checkAR.R")
source("libmain_fordebug.R")

# Acceptance ratios
nbeta<-dim(keepAR)[2]-3
colnames(keepAR)<-c("kiter","kernel","uiter",model@name.modpar)
keepAR <- as.data.frame(keepAR)
gtab<-cbind(keepAR[keepAR$kernel==1, 1:4], param="Joint kernel 1")
colnames(gtab)[4]<-"AR"
gtab2<-NULL
for(ipar in 1:nbeta) {
  gtab1<-cbind(keepAR[keepAR$kernel>1 & keepAR$uiter==3, c(1:3,3+ipar)], param=model@name.modpar[ipar])
  xtab <- cbind(keepAR[keepAR$kernel==23 & keepAR$uiter==1, c(1:3,3+ipar)], param=model@name.modpar[ipar])
  colnames(gtab1)[4]<-colnames(xtab)[4]<-"AR"
  gtab<-rbind(gtab, gtab1)
  gtab2<-rbind(gtab2,xtab)
}
gtab<-gtab[!is.na(gtab$AR),]
gpl <- ggplot(gtab, aes(x=kiter, y=AR, group=param, colour=as.factor(param))) + geom_line() + geom_line(data=gtab2, aes(x=kiter, y=AR, group=param, colour=as.factor(param)), linetype="dashed") + ylab("Acceptance ratio") + theme_minimal() + scale_color_viridis(discrete=TRUE) + guides(colour=guide_legend(title="Parameter/Step"))


cairo_ps(file = "acceptanceRates_pdseul250.eps", onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
print(gpl)
dev.off()

#################################################
# PK/PD with Alexandra's changes

# modified functions: saemix, initialiseMainAlgo, estep, mstep, compute.LLy, compute.Uy => added .alex
# all other functions are taken from saemix

#source("main_2longi.R")

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
  ypred[ytype==2] <- Emax[ytype==2]*ypred[ytype==2]/(ypred[ytype==2]+EC50[ytype==2]) # changed here
  
  return(ypred)
}

# Proportional error model
data_pkpd <- data_pkpd.prop
saemix.dataPKPD.prop<-saemixData(name.data=data_pkpd, name.group=c("id"), name.predictors=c("time","dose"), name.response="obs",name.ytype = "ytype")

pkpdmodel.prop<-saemixModel(model=LONGIjointmodel,description="joint pkpd",modeltype="structural",
                            psi0=matrix(param,ncol=5,byrow=TRUE,dimnames=list(NULL, c("ka","V","Cl","Emax","EC50"))),
                            transform.par=c(1,1,1,1,1), covariance.model=diag(c(1,1,1,1,1)),fixed.estim = c(1,1,1,1,1),
                            #                          omega.init = diag(c(1.96,0.02,0.09,0.02,0.04)),error.model = c("proportional","proportional"),error.init = c(0,1,0,1),
                            omega.init = diag(rep(0.5,5)),error.model = c("proportional","proportional"),error.init = c(0,1,0,1),
                            name.sigma = c("a.1","b.1","a.2","b.2"))

model <- pkpdmodel.prop
data <- saemix.dataPKPD.prop
control <- saemix.options.prop

######################################################
# Keeping track of AR
source("functions_alex_checkAR.R")
source("alexmain_fordebug.R")
postscript(file="convplot_pkpd250.eps",horizontal=T)
convplot.infit(allpar,saemix.options$nbiter.saemix[1],niter=(kiter-2))
dev.off()

# Acceptance ratios
nbeta<-dim(keepAR)[2]-3
colnames(keepAR)<-c("kiter","kernel","uiter",model@name.modpar)
keepAR <- as.data.frame(keepAR)
gtab<-cbind(keepAR[keepAR$kernel==1, 1:4], param="Joint kernel 1")
colnames(gtab)[4]<-"AR"
gtab2<-NULL
for(ipar in 1:nbeta) {
  gtab1<-cbind(keepAR[keepAR$kernel>1 & keepAR$uiter==3, c(1:3,3+ipar)], param=model@name.modpar[ipar])
  xtab <- cbind(keepAR[keepAR$kernel==23 & keepAR$uiter==1, c(1:3,3+ipar)], param=model@name.modpar[ipar])
  colnames(gtab1)[4]<-colnames(xtab)[4]<-"AR"
  gtab<-rbind(gtab, gtab1)
  gtab2<-rbind(gtab2,xtab)
}
gtab<-gtab[!is.na(gtab$AR),]
gpl <- ggplot(gtab, aes(x=kiter, y=AR, group=param, colour=as.factor(param))) + geom_line() + geom_line(data=gtab2, aes(x=kiter, y=AR, group=param, colour=as.factor(param)), linetype="dashed") + ylab("Acceptance ratio") + theme_minimal() + scale_color_viridis(discrete=TRUE) + guides(colour=guide_legend(title="Parameter/Step"))

cairo_ps(file = "acceptanceRates_pkpd250_pb.eps", onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
print(gpl)
dev.off()

######################################################
# Adding test on etas in the first kernel - doesn't change the acceptance rate for the first kernel
model <- pkpdmodel.prop
data <- saemix.dataPKPD.prop
control <- saemix.options.prop

source("functions_alex_checkAR2.R")
source("alexmain_fordebug.R")

nbeta<-dim(keepAR)[2]-3
colnames(keepAR)<-c("kiter","kernel","uiter",model@name.modpar)
keepAR <- as.data.frame(keepAR)
keep1 <- keepAR[keepAR$kernel==1, 1:4]
xlab1<- paste(keep1$kiter, keep1$uiter, sep=".")
keep1<-keep1[!duplicated(xlab1),]
gtab<-cbind(keep1, param="Joint kernel 1")
colnames(gtab)[4]<-"AR"
gtab2<-NULL
for(ipar in 1:nbeta) {
  gtab1<-cbind(keepAR[keepAR$kernel>1 & keepAR$uiter==3, c(1:3,3+ipar)], param=model@name.modpar[ipar])
  xtab <- cbind(keepAR[keepAR$kernel==23 & keepAR$uiter==1, c(1:3,3+ipar)], param=model@name.modpar[ipar])
  colnames(gtab1)[4]<-colnames(xtab)[4]<-"AR"
  gtab<-rbind(gtab, gtab1)
  gtab2<-rbind(gtab2,xtab)
}
gtab<-gtab[!is.na(gtab$AR),]
gpl <- ggplot(gtab, aes(x=kiter, y=AR, group=param, colour=as.factor(param))) + geom_line() + geom_line(data=gtab2, aes(x=kiter, y=AR, group=param, colour=as.factor(param)), linetype="dashed") + ylab("Acceptance ratio") + theme_minimal() + scale_color_viridis(discrete=TRUE) + guides(colour=guide_legend(title="Parameter/Step"))
gpl


######################################################
# Splitting the acceptance for the first kernel - better AR for the first kernel but still poor estimates
model <- pkpdmodel.prop
data <- saemix.dataPKPD.prop
control <- saemix.options.prop

source("alex_estep_outcome.R")
source("functions_alex_checkAR3.R")
source("alexmain_fordebug.R")

######################################################