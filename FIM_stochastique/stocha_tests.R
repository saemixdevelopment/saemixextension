rootDir = "C:/Users/lucie.fayette/Documents/Stage_M2/saemixextension"
scriptDir = file.path(rootDir, "FIM_stochastique") 

outPutDir = file.path(scriptDir, "stocha_outputs")
dir.create(outPutDir)

source(file.path(scriptDir,"Source.R"))

dataPath = file.path(rootDir, "data")

# library(rlang) #is_empty
library(matrixcalc)# elimination.matrix

library(dplyr) # %>%

CI.lin_sa <- function(OBJ){
  
  if(OBJ@options$fim&&OBJ@options$fim.sa){
    tab.sa = OBJ@results@conf.int.sa
    colnames(tab.sa) = c("name", "estimate", "se.sa", "cv.sa", "lower.sa", "upper.sa" )
    
    tab.lin_sa = left_join(OBJ@results@conf.int, tab.sa)
  }else if (OBJ@options$fim){
    tab.lin_sa = OBJ@results@conf.int
  }else if (OBJ@options$fim.sa){
    tab.lin_sa = OBJ@results@conf.int.sa
  }else {
    tab.lin_sa = NULL
  }
  
  return(tab.lin_sa)
}

writeCI.lin_sa <- function(OBJ, CIdir, Name){
  
  if(OBJ@options$fim&&OBJ@options$fim.sa){
    tab.sa = OBJ@results@conf.int.sa
    colnames(tab.sa) = c("name", "estimate", "se.sa", "cv.sa", "lower.sa", "upper.sa" )
    
    tab.lin_sa = left_join(OBJ@results@conf.int, tab.sa)
  }else if (OBJ@options$fim){
    tab.lin_sa = OBJ@results@conf.int
  }else if (OBJ@options$fim.sa){
    tab.lin_sa = OBJ@results@conf.int.sa
    colnames(tab.lin_sa) = c("name", "estimate", "se.sa", "cv.sa", "lower.sa", "upper.sa" )
    
  }else {
    tab.lin_sa = NULL
  }
  
  tab.lin_sa = tab.lin_sa %>% 
    mutate_if(is.numeric, round, 2)
  write.table(x = tab.lin_sa, file = paste0(CIdir, "/", Name,  "_CI.txt"), 
              row.names = FALSE, 
              quote = FALSE, fileEncoding = "UTF-8", sep=" & ")
}


##################################################"
###--- Theophylline pharmacokinetics -> too sparse

theo.saemix = read.csv(paste0(dataPath, "/theo.saemix.tab"),
                       header=TRUE,sep=" ")

theo.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
                      name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"), 
                      name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), 
                      name.X="Time")

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


# Model with a covariate effect (effect of Weight on CL)
theo.model<-saemixModel(model=model1cpt,
                        description="One-compartment model with first-order absorption", 
                        psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, 
                                                                                              c("ka","V","CL"))),
                        transform.par=c(1,1,1), 
                        covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE))

theo.option <-list(save=FALSE,save.graphs=FALSE,  fim=TRUE, fim.sa=TRUE)
theo.fit_stoch = saemix.fim_stoch(theo.model,theo.data,theo.option)

CI.lin_sa(theo.fit_stoch)
writeCI.lin_sa(theo.fit_stoch, outPutDir, "Theophylline_CI_lin_sa")


##################################################"
###--- Emax

PD1.saemix = read.csv(paste0(dataPath, "/PD1.saemix.tab"),
                      header=TRUE,sep=" ")

head(PD1.saemix)
unique(PD1.saemix$subject)


saemix.data1<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
                         name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
                         units=list(x="mg",y="-",covariates="-"))

modelemax<-function(psi,id,xidep) {
  # input:
  #   psi : matrix of parameters (3 columns, E0, Emax, EC50)
  #   id : vector of indices 
  #   xidep : dependent variables (same nb of rows as length of id)
  # returns:
  #   a vector of predictions of length equal to length of id
  dose<-xidep[,1]
  e0<-psi[id,1]
  emax<-psi[id,2]
  e50<-psi[id,3]
  f<-e0+emax*dose/(e50+dose)
  return(f)
}

# Compare models with and without covariates with LL by Importance Sampling
model1<-saemixModel(model=modelemax,description="Emax growth model", 
                    psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,
                                                                                   c("E0","Emax","EC50"))), transform.par=c(1,1,1),
                    covariate.model=matrix(c(0,0,0), ncol=3,byrow=TRUE),fixed.estim=c(1,1,1))
model2<-saemixModel(model=modelemax,description="Emax growth model", 
                    psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL, 
                                                                                   c("E0","Emax","EC50"))), transform.par=c(1,1,1),
                    covariate.model=matrix(c(0,0,1), ncol=3,byrow=TRUE),fixed.estim=c(1,1,1))

saemix.options<-list(algorithms=c(0,1,1),nb.chains=3,seed=765754, 
                     nbiter.saemix=c(500,300),save=FALSE,save.graphs=FALSE, displayProgress=FALSE,
                     fim=TRUE, fim.sa=TRUE)
fit1.stoch <- saemix.fim_stoch(model1,saemix.data1,saemix.options)
fit2.stoch <- saemix.fim_stoch(model2,saemix.data1,saemix.options)



# Comparison observed fim (stochastic) and expected fim (linearization)
CI.lin_sa(fit1.stoch)
writeCI.lin_sa(fit1.stoch, outPutDir, "PD1_fit1_CI_lin_sa")

CI.lin_sa(fit2.stoch)
writeCI.lin_sa(fit2.stoch, outPutDir, "PD1_fit2_CI_lin_sa")


##################################################"
###---  Discrete
toenail.saemix = read.csv(paste0(dataPath, "/toenail.saemix.tab"),
                          header=TRUE,sep=" ")

head(toenail.saemix)
length(unique(toenail.saemix$id))


toenail.saemix.data<-saemixData(name.data=toenail.saemix,name.group=c("id"),name.predictors=c("time","y"), name.response="y",
                                name.covariates=c("treatment"),name.X=c("time"))

# Model
binary.model<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  logpdf<-rep(0,length(tim))
  P.obs = (y==0)*(1-pevent)+(y==1)*pevent
  logpdf <- log(P.obs)
  return(logpdf)
}

binary.saemix.model<-saemixModel(model=binary.model,description="Binary model",
                                 modeltype="likelihood",
                                 psi0=matrix(c(-0.5,-.15,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("alpha","beta"))),
                                 transform.par=c(0,0), covariate.model=c(0,1),
                                 covariance.model=matrix(c(1,0,0,1),ncol=2), omega.init=diag(c(0.5,0.3)))

binary.saemix.model.iiv1<-saemixModel(model=binary.model,description="Binary model",
                                      modeltype="likelihood",
                                      psi0=matrix(c(-0.5,-.15,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("alpha","beta"))),
                                      transform.par=c(0,0), covariate.model=c(0,1),
                                      covariance.model=matrix(c(1,0,0,0),ncol=2), omega.init=diag(c(0.5,0.3)))

binary.saemix.model.nocov<-saemixModel(model=binary.model,description="Binary model",
                                       modeltype="likelihood",
                                       psi0=matrix(c(-0.5,-.15,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("alpha","beta"))),
                                       transform.par=c(0,0), covariate.model=c(0,0),
                                       covariance.model=matrix(c(1,0,0,1),ncol=2), omega.init=diag(c(0.5,0.3)))


saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE,
                     nb.chains=10, fim=FALSE, fim.sa=TRUE, map=FALSE)
binary.fit.stoch <- saemix.fim_stoch(binary.saemix.model,toenail.saemix.data,saemix.options)
binary.fit.nocov.stoch <- saemix.fim_stoch(binary.saemix.model.nocov,toenail.saemix.data,saemix.options)
binary.fit.iiv1.stoch <- saemix.fim_stoch(binary.saemix.model.iiv1,toenail.saemix.data,saemix.options)

CI.lin_sa(binary.fit.stoch)
writeCI.lin_sa(binary.fit.stoch, outPutDir, "Binary_fit_CI_lin_sa")

CI.lin_sa(binary.fit.nocov.stoch)
writeCI.lin_sa(binary.fit.nocov.stoch, outPutDir, "Binary_fit_nocov_CI_lin_sa")

CI.lin_sa(binary.fit.iiv1.stoch)
writeCI.lin_sa(binary.fit.iiv1.stoch, outPutDir, "Binary_fit_iiv1_CI_lin_sa")

# model = binary.saemix.model.iiv1
# data= toenail.saemix.data
# control =saemix.options

##################################################"
###---  Oxford boys
oxboys.saemix = read.csv(paste0(dataPath, "/oxboys.saemix.tab"),
                         header=TRUE,sep=" ")
head(oxboys.saemix)
length(unique(oxboys.saemix$Subject))
nrow(oxboys.saemix)/26

oxboy.data<-saemixData(name.data=oxboys.saemix,header=TRUE,
                       name.group=c("Subject"),name.predictors=c("age"),name.response=c("height"),
                       units=list(x="yr",y="cm"))

growth.linear<-function(psi,id,xidep) {
  x<-xidep[,1]
  base<-psi[id,1]
  slope<-psi[id,2]
  f<-base+slope*x
  return(f)
}
oxboy.model<-saemixModel(model=growth.linear,description="Linear model",
                         psi0=matrix(c(140,1),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
                         transform.par=c(1,0),covariance.model=matrix(c(1,1,1,1),ncol=2,byrow=TRUE), 
                         error.model="constant")


saemix.options<-list(algorithms=c(1,1,1),nb.chains=1,seed=201004,save=FALSE,
                     save.graphs=FALSE, displayProgress=FALSE,
                     fim=TRUE, fim.sa=TRUE)
oxboy.fit.fim_stoch = saemix.fim_stoch(oxboy.model,oxboy.data,saemix.options)

CI.lin_sa(oxboy.fit.fim_stoch)
writeCI.lin_sa(oxboy.fit.fim_stoch, outPutDir, "OxBoys_CI_lin_sa")



##################################################"
###---  knee
knee.saemix = read.csv(paste0(dataPath, "/knee.saemix.tab"),
                         header=TRUE,sep=" ")
head(knee.saemix)
length(unique(knee.saemix$id))
nrow(knee.saemix)/127

knee.data<-saemixData(name.data=knee.saemix,name.group=c("id"),
                      name.predictors=c("y", "time"), name.X=c("time"),
                      name.covariates = c("Age","Sex","treatment","Age2"),
                      units=list(x="d",y="", covariates=c("yr","-","-","yr2")), verbose=FALSE)

# Fitting PO model - ok but fit not very good (see diagnostics)
ordinal.model<-function(psi,id,xidep) {
  y<-xidep[,1]
  time<-xidep[,2]
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  alp4<-psi[id,4]
  beta<-psi[id,5]
  
  logit1<-alp1 + beta*time
  logit2<-logit1+alp2
  logit3<-logit2+alp3
  logit4<-logit3+alp4
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  pge4<-exp(logit4)/(1+exp(logit4))
  pobs = (y==1)*pge1+(y==2)*(pge2 - pge1)+(y==3)*(pge3 - pge2)+(y==4)*(pge4 - pge3)+(y==5)*(1 - pge4)
  logpdf <- log(pobs)
  
  return(logpdf)
}
# simulate function
simulateOrdinal<-function(psi,id,xidep) {
  y<-xidep[,1]
  time<-xidep[,2]
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  alp4<-psi[id,4]
  beta<-psi[id,5]
  
  logit1<-alp1 + beta*time
  logit2<-logit1+alp2
  logit3<-logit2+alp3
  logit4<-logit3+alp4
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  pge4<-exp(logit4)/(1+exp(logit4))
  x<-runif(length(time))
  ysim<-1+as.integer(x>pge1)+as.integer(x>pge2)+as.integer(x>pge3)+as.integer(x>pge4)
  return(ysim)
}


knee.model.PO<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                             simulate.function=simulateOrdinal, 
                             psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5, byrow=TRUE, 
                                         dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), transform.par=c(0,1,1,1,1),
                             omega.init=diag(c(100, 1, 1, 1, 1)), covariance.model = diag(c(1,0,0,0,1)), verbose=FALSE)

# Fitting
knee.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, 
                   fim=FALSE, fim.sa=TRUE,
                   nb.chains=10, nbiter.saemix=c(600,100), print=FALSE)

knee.fit_stoch<-saemix.fim_stoch(knee.model.PO,knee.data,knee.options)

CI.lin_sa(knee.fit_stoch)
writeCI.lin_sa(knee.fit_stoch, outPutDir, "Knee_CI_lin_sa")


# model = knee.model.PO
# data = knee.data
# control = knee.options


