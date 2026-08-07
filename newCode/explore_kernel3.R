# Folders
saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
datDir<-file.path(saemixDir,"data")

# Libraries
library(MASS)
library(rlang)
library(ggplot2)
library(gridExtra)

# Sourcing saemix functions
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixData-methods.R"))
source(file.path(progDir,"SaemixData-methods_covariates.R"))
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

###################################
# Theophylline dataset
theo.saemix<-read.table(file.path(datDir, "theo.saemix.tab"), header=TRUE)
theo.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
                      name.group=c("Id"),name.predictors=c("Dose","Time"),
                      name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
                      units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
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
theo.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",
                        psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE,  dimnames=list(NULL, c("ka","V","CL"))),
                        transform.par=c(1,1,1))
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
theo.fit<-saemix(theo.model,theo.data,saemix.options)

# Alternate E-step - changing 3rd kernel
# still a large increase in omega(ka)
source(file.path(saemixDir,"newCode","main_altinitialiseMainAlgo.R"))
source(file.path(saemixDir,"newCode","main_altEstep.R"))
theo.fit2<-saemix(theo.model,theo.data,saemix.options)

# Not working: inflation of omega (here, omega(ka), associated with a slight increase of ka)
# later simulations show a general problem :-/
coef(theo.fit)$fixed
coef(theo.fit2)$fixed
sqrt(diag(theo.fit@results@omega))
sqrt(diag(theo.fit2@results@omega))

# Even worse when we add more iterations to the 3rd kernel, now ka is very biased
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nbiter.mcmc=c(2,2,10,0))
theo.fit2<-saemix(theo.model,theo.data,saemix.options)
coef(theo.fit2)$fixed

# Just one iteration - increases but less
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nbiter.mcmc=c(2,2,1,0))
theo.fit2<-saemix(theo.model,theo.data,saemix.options)

# No third kernel
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nbiter.mcmc=c(2,2,0,0))
theo.fit3<-saemix(theo.model,theo.data,saemix.options)
coef(theo.fit)$fixed
coef(theo.fit3)$fixed

####### Random walk as kernel 2 but over all etas
# still a large increase in omega(ka) and a small increase in ka
source(file.path(saemixDir,"newCode","main_rwalkEstep.R"))
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
theo.fit.rw<-saemix(theo.model,theo.data,saemix.options)
theo.fit.rw

coef(theo.fit)$fixed
coef(theo.fit.rw)$fixed
sqrt(diag(theo.fit@results@omega))
sqrt(diag(theo.fit.rw@results@omega))

###################################
# Comparing kernels on simulated data
parpop<-c(5,30,500,2)
nampar<-c("E0","Emax","ED50","gamma")
omega<-diag(c(0.09,0.49,0.49))
omega[3,2]<-omega[2,3]<-0.245
respar<-c(0.1)

# Massive increase in Omega(ED50, Emax) and much more variability in estimates of the corresponding fixed effects including many absurd values, as a result mean and median both inflated
tab1<-read.csv(file.path("/home/eco/work/saemix/saemixextension/simulationSuite","cont","hillRichProp/results", "260512_eco_alt4Estep_hillRichProp_defaultTrue.res"), sep=" ", header=T)
summary(tab1)

# Comparison with v35
tab2<-read.csv(file.path("/home/eco/work/saemix/saemixextension/simulationSuite","cont","hillRichProp/results", "260512_eco_v35_hillRichProp_defaultTrue.res"), sep=" ", header=T)
summary(tab2)

# No difference between v35 and v31
tab3<-read.csv(file.path("/home/eco/work/saemix/saemixextension/simulationSuite","cont","hillRichProp/results", "220913_eco_cran31_hillRichProp_defaultTrue.res"), sep=" ", header=T)
for(icol in 1:dim(tab3)[2])
  if(max(abs(tab3[,icol]-tab2[,icol]))>10-6) cat("Check estimates for",namcol(tab3)[icol],"\n")

# Removing 3rd kernel
tab4<-read.csv(file.path("/home/eco/work/saemix/saemixextension/simulationSuite","cont","hillRichProp/results", "260512_eco_no3Estep_hillRichProp_defaultTrue.res"), sep=" ", header=T)
summary(tab4)

# random walk as in 2nd kernel but over all parameters
tab5<-read.csv(file.path("/home/eco/work/saemix/saemixextension/simulationSuite","cont","hillRichProp/results", "260520_eco_rwalkEstep_hillRichProp_defaultTrue.res"), sep=" ", header=T)
summary(tab5)

###################################
