workDir <- "/home/eco/work/saemix/saemixextension/lucasMultiple"
setwd(workDir)
datDir <- file.path(workDir, "exemple_emmanuelle")

library(tidyverse)
library(Rcpp)

for (i in 1:2){
  fichiers = sort(list.files(path = "R" , pattern = "\\." , full.names = TRUE))
  sapply(fichiers,source)
}

fichiers_multi = list.files(path="R_multi", pattern = "\\." , full.names = TRUE)
  sapply(fichiers_multi,source)

source(file.path(workDir, "R_eco","func_FIM_multi.R"))

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

###########################################################
          # Simulation Numéro 2
############################################################

data_sim2 = read.csv(file.path(datDir,"basicPKPD_Sim2.csv"), header = TRUE)
# Pb avec un modèle exponentiel quand on a des données à t=0 (my bad)
data_sim2 <- data_sim2[data_sim2$tim>0,]

# note: j'ai essayé de ne retirer que les données de PK mais du coup b part à dache 
# vérifier qu'on peut avoir ytype=2 comme lap remière réponse vue parce que bizarre...

saemix.data <- saemixData(name.data=data_sim2,header=TRUE,sep=" ",na=NA,
          name.group=c("id"),name.predictors=c("dose","tim"), name.response=c("y"), name.X="tim",
					units = list(x = "mg/L", y = "hr"),name.ytype = "ytype")

# paramètres mu 
psi_ini =  c(ka = 1, V = 20, CL = 2, Emax = 100, EC50 = 5)
x.grid05 <- gridfuncd(data= data_sim2, psi0=psi_ini, e.obj = 0.05 , class = "ytype", model = modelPKPD)

saemix.model.exact <- saemixModelPBPK(psi0=psi_ini, 
                                x.grid=x.grid05,
                                simulations=unique(data_sim2$id), 
                                output=NULL,
				covariate.model = NULL  , 
				covariance.model=diag(1,5,5), # Matrice diagonale , pas de covariance entre les paramètres
				error.model = c("proportional" , "constant"), # proportionnel pour PK et constant pour PD
				modeltype = c("structural","structural") , 
				model.interpolate = FALSE , # on n'utilise pas la grille
				name.response = c("y1","y2") , 
				model = modelPKPD) # le model structural 

options <- list(seed=54333, map=F, fim=F, ll.is=F,
                nbiter.saemix = c(300,200), nbiter.sa=20,
                nmc.is = 5000, nu.is=5, r.is=1, nbdisplay=25,
                displayProgress=T, print.is=F, save.graphs=F,
		warnings=TRUE, directory = "exactFit2")				
				
time.beg1 <- Sys.time()
fit_sim2.exact =  saemix.multi(saemix.model.exact, saemix.data, options)
time.end1 <- Sys.time()
time.fit2.exact <- time.end1-time.beg1

save(fit_sim2.exact, file=file.path(workDir,"exactFit2","fit_sim2.exact.Robj"))

###########################################################
# FIM
tryfim <- fim.saemix(fit_sim2.exact)
print(tryfim)

# LL IS
tryfim <- llis.saemix(tryfim)
tryfim@results@ll.lin
tryfim@results@ll.is

# LL GQ
tryfim <- llgq.saemix(tryfim)
tryfim@results@ll.gq

# Individual MAP estimates
trymap <- map.saemix(tryfim)

# Conditional distributions
trycond <- conddist.saemix(trymap)

# Predictions
trypred<-saemix.predict(trycond)

# Simulations
object <- trycond
x<-simulate(trycond, nsim=2)

###########################################################
# Graphs
x<-compute.sres(object)
plot(x@data@data$tim,x@results@predictions$npde)
abline(h=0)

# npde
is.y1 <- x@data@data$ytype==1
dat.y1 <- x@data@data[is.y1,]
dat.y2 <- x@data@data[!(is.y1),]
sim.y1 <- x@sim.data@datasim[rep(is.y1,x@sim.data@nsim),-c(2)]
sim.y2 <- x@sim.data@datasim[rep(!(is.y1),x@sim.data@nsim),-c(2)]

npde.y1 <- autonpde(dat.y1, sim.y1, iid="id", ix="tim", iy="y", boolsave=FALSE)
npde.y2 <- autonpde(dat.y2, sim.y2, iid="id", ix="tim", iy="y", boolsave=FALSE)

plot(npde.y1)
plot(npde.y2)

# VPC
plot(npde.y1, plot.type="vpc")
plot(npde.y2, plot.type="vpc")

# Graphes individuels
ypl <- trypred@data@data[,c(2,4,5,9),]
ypl$ipred <- trypred@results@ipred
ggplot(ypl[ypl$ytype==1 & ypl$id<7,],aes(x=tim, y=y, group=id)) + geom_point() + geom_line(aes(x=tim, y=ipred)) + facet_wrap(.~id, nrow=2, ncol=3)
ggplot(ypl[ypl$ytype==2 & ypl$id<7,],aes(x=tim, y=y, group=id)) + geom_point() + geom_line(aes(x=tim, y=ipred)) + facet_wrap(.~id, nrow=2, ncol=3)

###########################################################


saemix.model.exact <- saemixModelPBPK(psi0=psi_ini, 
                                      x.grid=x.grid05,
                                      simulations=unique(data_sim2$id), 
                                      output=NULL,
                                      covariate.model = NULL  , 
                                      covariance.model=diag(1,5,5), # Matrice diagonale , pas de covariance entre les paramètres
                                      error.model = c("proportional" , "constant"), # proportionnel pour PK et constant pour PD
                                      modeltype = c("structural","structural") , 
                                      model.interpolate = FALSE , # on n'utilise pas la grille
                                      name.response = c("y1","y2") , 
                                      model = modelPKPD) # le model structural 

options <- list(seed=54333, map=T, fim=T, ll.is=T,
                nbiter.saemix = c(300,200), nbiter.sa=20,
                nmc.is = 5000, nu.is=5, r.is=1, nbdisplay=25,
                displayProgress=T, print.is=F, save.graphs=F,
                warnings=TRUE, directory = "exactFit2")				

time.beg1 <- Sys.time()
yfit.full <-  saemix.multi(saemix.model.exact, saemix.data, options)
time.end1 <- Sys.time()

