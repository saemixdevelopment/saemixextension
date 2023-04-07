# Jackknife

saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
datDir<-file.path(saemixDir,"data")

source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixRes.R"))
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
source(file.path(progDir,"func_npde.R"))
source(file.path(progDir,"func_estimParam.R"))
source(file.path(progDir,"backward.R"))
source(file.path(progDir,"forward.R"))
source(file.path(progDir,"func_stepwise.R"))
source(file.path(progDir,"stepwise.R"))
source(file.path(progDir,"func_compare.R"))
source(file.path(progDir,"func_bootstrap.R"))
source(file.path(progDir,"func_exploreData.R"))
source(file.path(progDir,"func_discreteVPC.R"))

# Minimal theophylline example
theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T)
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
                        name.group=c("Id"),name.predictors=c("Dose","Time"),
                        name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
                        units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time", verbose = FALSE)

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

# Model with covariate Weight
saemix.model<-saemixModel(model=model1cpt,modeltype="structural",
                          description="One-compartment model with first-order absorption",
                          psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
                          transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE), verbose=FALSE)

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

# Leave-one-out

saemix.object <- saemix.fit
sopt<-saemix.object@options
sopt$displayProgress<-FALSE
sopt$save <- FALSE
sopt$save.graphs <- FALSE

resorig <- summary(saemix.object@results)
origpar <- c(resorig$fixed.effects[,2], resorig$random.effects[,2])

idx.rho<-which(saemix.object@model@covariance.model[lower.tri(saemix.object@model@covariance.model)]==1)
origpar <- c(saemix.object@results@fixed.effects[c(saemix.object@results@indx.fix, saemix.object@results@indx.cov)], diag(saemix.object@results@omega)[saemix.object@results@indx.omega], saemix.object@results@omega[lower.tri(saemix.object@results@omega)][idx.rho], saemix.object@results@respar[saemix.object@results@indx.res])

zesuj<-unique(saemix.object@data@data$index)
xvec <- vector(mode="list", length=length(zesuj))
for(isuj in 1:length(zesuj)) {
  sdata <- subset(saemix.object@data, index!=zesuj[isuj])
  yfit <- saemix(saemix.object@model, sdata, control=sopt)
  xvec[[isuj]] <- summary(yfit@results)
}

# Computing Ii = (\hat(\theta) - \hat(\theta_{-i}))*(N-1)
iaccel <- thetai <- NULL
for(isuj in 1:length(zesuj)) {
  res<-xvec[[isuj]]
  ipar <- c(res$fixed.effects[,2], res$random.effects[,2])
  iaccel <- rbind(iaccel, (origpar-ipar))
  thetai <- rbind(thetai, ipar)
}
iaccel<-iaccel*(saemix.object@data@N-1)
xaccel <- colSums(iaccel**3)/(colSums(iaccel**2)**(3/2))/6

# Bootstrap
boot.case <- saemix.bootstrap(saemix.fit, method="case")

# same ordering as jackknife (will need to change the order to the same as bootstrap +++)
boot.case.order<- boot.case[,c(2:5,9,6:8)]
boot.case.order<- boot.case[,c(2:9)]

alpha <- 0.05 #Desired quantiles
u <- c(alpha/2, 1-alpha/2) 
boot.perc <- apply(boot.case.order, 2, quantile, u)

boot.unbias <- t(t(boot.case.order)-origpar)

par(mfrow=c(3,3))
for(i in 1:8) {
  hist(boot.case.order[,i])
  abline(v=origpar[i])
}

# BCa
# Accelerated Bootstrap CI
zu <- qnorm(u)

boot.bca <- NULL
for(icol in 1:8) {
  z0 <- qnorm(mean(boot.case.order[,icol] <= origpar[icol]))
  uadj <- pnorm(z0 + (z0+zu)/(1-xaccel[icol]*(z0+zu))) 
  boot.bca<-cbind(boot.bca,quantile(boot.case.order[,icol], uadj))
}

# Comparing bootstrap and jackknife estimates for a rough idea
summary(boot.case)
jackdist <- saemix.jackknife(saemix.fit, compute.likelihood = TRUE)
summary(jackdist)

boot.perc <- apply(boot.case.order, 2, quantile, u)

boot.bca
bca.percentile(saemix.fit, boot.case, jackdist)
  
bca.percentile(saemix.fit, boot.case)
