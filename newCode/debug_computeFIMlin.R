library(saemix)
library(numDeriv)

saemixDir <- "/home/eco/work/saemix/saemixextension"
source(file.path(saemixDir,"newCode","compute_derivFIMlin.R"))

###################################### Theophylline example (N=12, n=10)
data(theo.saemix)
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
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

# Model with covariates - 3 IIV
saemix.model<-saemixModel(model=model1cpt, description="One-compartment model with first-order absorption",
                          psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, 
                                      dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), 
                          covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
                          covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),
                          omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

saemix.options<-list(seed=39546,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

# Comparing computation of LL on eta and psi scale
if(FALSE) {
  source(file.path(saemixDir,"newCode","compute_derivFIMlin_etascale.R"))
  invfim1<-computeFIMinv.lin(saemix.fit)
  invfim2<-computeFIMinv.linEta(saemix.fit)
  
  cat("SE, new code, psi scale=",sqrt(diag(invfim1)), "\n") # Psi scale
  cat("SE, new code, eta scale=",sqrt(diag(invfim2))*c(saemix.fit@results@fixed.psi, rep(1,dim(fim2)[1]-3)), "\n") # Eta scale
  cat("SE, saemix 3.1=",sqrt(diag(solve(saemix.fit@results@fim))), "\n") # saemix 3.1
}

# Comparing LL and SE
invfim1<-computeFIMinv.lin(saemix.fit)
ll1<-computeLL.lin(saemix.fit)
cat("LL, saemix 3.1=",saemix.fit@results@ll.lin, "\n")
cat("LL, new code=",ll1, "\n")

cat("SE, saemix 3.1=",sqrt(diag(solve(saemix.fit@results@fim))),"\n")
cat("SE, new code=",sqrt(diag(invfim1)),"\n")

# Inversion by bloc versus the whole matrix
cat("bloc inversion", sqrt(diag(solve(solve(invfim1[1:4,1:4])))), sqrt(diag(solve(solve(invfim1[5:10,5:10])))),"\n")

# Comparing resulting confidence intervals
saemix.fit.fim <- computeFIM.lin(saemix.fit)
saemix.fit.fim@results@conf.int
saemix.fit@results@conf.int

#############################
# Model with only 2 IIV
smxmod.noIIVka<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",
                          psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, 
                                      dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), 
                          covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
                          covariance.model=matrix(c(0,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),
                          omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")
fit.noIIVka<-saemix(smxmod.noIIVka,saemix.data,saemix.options)
fit.noIIVka.fim <- computeFIM.lin(fit.noIIVka)

invfim.noIIVka<-computeFIMinv.lin(fit.noIIVka)
cat("SE, saemix 3.1=",sqrt(diag(solve(fit.noIIVka@results@fim))),"\n")
cat("SE, new code=",sqrt(diag(invfim.noIIVka)),"\n")

# Comparing resulting confidence intervals
fit.noIIVka.fim <- computeFIM.lin(fit.noIIVka)
fit.noIIVka.fim@results@conf.int
fit.noIIVka@results@conf.int

#############################
# Single IIV
smxmod.oneIIV<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",
                            psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, 
                                        dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), 
                            covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
                            covariance.model=diag(c(1,0,0)),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")
fit.oneIIV<-saemix(smxmod.oneIIV,saemix.data,saemix.options)
invfim.oneIIV<-computeFIMinv.lin(fit.oneIIV)
cat("SE, saemix 3.1=",sqrt(diag(solve(fit.oneIIV@results@fim))),"\n")
cat("SE, new code=",sqrt(diag(invfim.oneIIV)),"\n")

# Additive error model 
smxmod.addError<-saemixModel(model=model1cpt,psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, 
                                      dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), 
                          covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
                          covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),
                          omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
fit.addError<-saemix(smxmod.addError,saemix.data,saemix.options)

invfim.addError<-computeFIMinv.lin(fit.addError)
cat("SE, saemix 3.1=",sqrt(diag(solve(fit.addError@results@fim))),"\n")
cat("SE, new code=",sqrt(diag(invfim.addError)),"\n")
lladd<-computeLL.lin(fit.addError)
cat("LL saemix 3.1",fit.addError@results@ll.lin,"- new LL",lladd,"\n")

# Inversion by bloc versus the whole matrix
cat("SE, new code=",sqrt(diag(invfim.addError)),"\n")
cat("bloc inversion", sqrt(diag(solve(solve(invfim.addError[1:4,1:4])))), sqrt(diag(solve(solve(invfim.addError[5:9,5:9])))),"\n")

# Proportional error model 
smxmod.propError<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",
                             psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, 
                                         dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), 
                             covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
                             covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),
                             omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="proportional")
fit.propError<-saemix(smxmod.propError,saemix.data,saemix.options)

invfim.propError<-computeFIMinv.lin(fit.propError)
cat("SE, saemix 3.1=",sqrt(diag(solve(fit.propError@results@fim))),"\n")
cat("SE, new code=",sqrt(diag(invfim.propError)),"\n")
llprop<-computeLL.lin(fit.propError)
cat("LL saemix 3.1",fit.propError@results@ll.lin,"- new LL",llprop,"\n")

###################################### Emax, 'N-rich' design, N=100, n=3
# Very similar estimates with the Hessian
data(PD1.saemix)

saemix.pd1<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
                        name.predictors=c("dose"),name.response=c("response"),
                        name.covariates=c("gender"), units=list(x="mg",y="-",covariates=c("-")))

modelemax<-function(psi,id,xidep) {
  dose<-xidep[,1]
  e0<-psi[id,1]
  emax<-psi[id,2]
  e50<-psi[id,3]
  f<-e0+emax*dose/(e50+dose)
  return(f)
}

# Compare models with and without covariates with LL by Importance Sampling
model1<-saemixModel(model=modelemax,description="Emax growth model", 
                    psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,c("E0","Emax","EC50"))), transform.par=c(1,1,1),
                    covariate.model=matrix(c(0,0,0), ncol=3,byrow=TRUE),fixed.estim=c(1,1,1))

model2<-saemixModel(model=modelemax,description="Emax growth model", 
                    psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL, c("E0","Emax","EC50"))), transform.par=c(1,1,1),
                    covariate.model=matrix(c(0,0,1), ncol=3,byrow=TRUE),fixed.estim=c(1,1,1))

emx.options<-list(nb.chains=3,seed=765754, nbiter.saemix=c(500,300),save=FALSE,save.graphs=FALSE, displayProgress=FALSE)

fit1<-saemix(model1,saemix.pd1,emx.options)
fit2<-saemix(model2,saemix.pd1,emx.options)

fit.emax1 <- computeFIMinv.lin(fit1)
fit.emax2 <- computeFIMinv.lin(fit2)
cat("SE, saemix 3.1=",sqrt(diag(solve(fit1@results@fim))),"\n")
cat("SE, new code=",sqrt(diag(fit.emax1)),"\n")
cat("SE, saemix 3.1=",sqrt(diag(solve(fit2@results@fim))),"\n")
cat("SE, new code=",sqrt(diag(fit.emax2)),"\n")

# Comparing resulting confidence intervals
fit1.fim <- computeFIM.lin(fit1)
fit1.fim@results@conf.int
fit1@results@conf.int

