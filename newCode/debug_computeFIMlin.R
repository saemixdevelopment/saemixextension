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
theo.fit<-saemix(saemix.model,saemix.data,saemix.options)

# Comparing resulting confidence intervals
theo.fit.fim <- computeFIM.lin(theo.fit)
theo.fit.fim@results@conf.int
theo.fit@results@conf.int

# Comparing computation of LL on eta and psi scale
if(FALSE) {
  source(file.path(saemixDir,"newCode","compute_derivFIMlin_etascale.R"))
  invfim1<-computeFIMinv.lin(theo.fit)
  invfim2<-computeFIMinv.linEta(theo.fit)
  
  cat("SE, new code, psi scale=",sqrt(diag(invfim1)), "\n") # Psi scale
  cat("SE, new code, eta scale=",sqrt(diag(invfim2))*c(theo.fit@results@fixed.psi, rep(1,dim(fim2)[1]-3)), "\n") # Eta scale
  cat("SE, saemix 3.1=",sqrt(diag(solve(theo.fit@results@fim))), "\n") # saemix 3.1
}

# Comparing LL and SE
invfim1<-computeFIMinv.lin(theo.fit)
ll1<-computeLL.lin(theo.fit)
cat("LL, saemix 3.1=",theo.fit@results@ll.lin, "\n")
cat("LL, new code=",ll1, "\n")

cat("SE, saemix 3.1=",sqrt(diag(solve(theo.fit@results@fim))),"\n")
cat("SE, new code=",sqrt(diag(invfim1)),"\n")

# Inversion by bloc versus the whole matrix
cat("bloc inversion", sqrt(diag(solve(solve(invfim1[1:4,1:4])))), sqrt(diag(solve(solve(invfim1[5:10,5:10])))),"\n")

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

###################################### Oxford boys (linear)
# Almost identical with the Hessian
data(oxboys.saemix)
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

saemix.options<-list(algorithms=c(1,1,1),nb.chains=1,seed=201004,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
oxboy.fit<-saemix(oxboy.model,oxboy.data,saemix.options)

# Comparing resulting confidence intervals
oxboy.fim <- computeFIM.lin(oxboy.fit)
oxboy.fim@results@conf.int
oxboy.fit@results@conf.int

############################################################################ 
# Supposedly equivalent expression using product of gradients... 
## not working, keeps giving a singular matrix that can't be inverted

source(file.path(saemixDir,"newCode","compute_jacobianFIM.R"))

computeJacFIMinv.lin(theo.fit)

############
jacky<-computeJacFIMinv.lin(oxboy.fit)

ifim<-MASS::ginv(hess2)
sqrt(diag(ifim))
pracma::inv(hess2)

# Checking jacobian
lis1 <- initialiseComputationLL(oxboy.fit)
ll0 <- saemix.LLlin(theta=lis1$theta, modeltype=lis1$modeltype, etype=lis1$etype, f0=lis1$f0, DF0=lis1$DF0, z=lis1$z, Mcov=lis1$Mcov, nfix=lis1$nfix, 
                         nomega=lis1$nomega, nind.obs=lis1$nind.obs, ndat.exp=0, idx.fix=lis1$idx.fix,idx.omega=lis1$idx.omega, idx.res=lis1$idx.res, 
                         tr.fix=lis1$tr.fix) 
# Very different results !!! ??? why ?
pracma::jacobian(saemix.LLlin,lis1$theta, modeltype=lis1$modeltype, etype=lis1$etype, f0=lis1$f0, DF0=lis1$DF0, z=lis1$z, Mcov=lis1$Mcov, nfix=lis1$nfix,  nomega=lis1$nomega, nind.obs=lis1$nind.obs, ndat.exp=0, idx.fix=lis1$idx.fix,idx.omega=lis1$idx.omega, idx.res=lis1$idx.res,  tr.fix=lis1$tr.fix) 
numDeriv::jacobian(saemix.LLlin,lis1$theta, modeltype=lis1$modeltype, etype=lis1$etype, f0=lis1$f0, DF0=lis1$DF0, z=lis1$z, Mcov=lis1$Mcov, nfix=lis1$nfix,  nomega=lis1$nomega, nind.obs=lis1$nind.obs, ndat.exp=0, idx.fix=lis1$idx.fix,idx.omega=lis1$idx.omega, idx.res=lis1$idx.res,  tr.fix=lis1$tr.fix, method.args=list(d=1e-6)) # decent
numDeriv::jacobian(saemix.LLlin,lis1$theta, modeltype=lis1$modeltype, etype=lis1$etype, f0=lis1$f0, DF0=lis1$DF0, z=lis1$z, Mcov=lis1$Mcov, nfix=lis1$nfix,  nomega=lis1$nomega, nind.obs=lis1$nind.obs, ndat.exp=0, idx.fix=lis1$idx.fix,idx.omega=lis1$idx.omega, idx.res=lis1$idx.res,  tr.fix=lis1$tr.fix, method.args=list(d=.Machine$double.eps^(1/3))) # horrible, why ? .Machine$double.eps^(1/3)=6.e-6 so within the proper range
jacprac <- pracma::jacobian(saemix.LLlin,lis1$theta, modeltype=lis1$modeltype, etype=lis1$etype, f0=lis1$f0, DF0=lis1$DF0, z=lis1$z, Mcov=lis1$Mcov, nfix=lis1$nfix,  nomega=lis1$nomega, nind.obs=lis1$nind.obs, ndat.exp=0, idx.fix=lis1$idx.fix,idx.omega=lis1$idx.omega, idx.res=lis1$idx.res,  tr.fix=lis1$tr.fix) 

# By hand: not very stable... 1e-6, 1e-7, 1e-8 give similar results
theta0<-lis1$theta
for(deps in c(1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)) {
  dLL<-c()
  for(i in 1:length(theta0)) {
    theta1<-theta0
    theta1[i]<-theta1[i]*(1+deps)
    dLL<-c(dLL,saemix.LLlin(theta=theta1, modeltype=lis1$modeltype, etype=lis1$etype, f0=lis1$f0, DF0=lis1$DF0, z=lis1$z, Mcov=lis1$Mcov, nfix=lis1$nfix, 
                            nomega=lis1$nomega, nind.obs=lis1$nind.obs, ndat.exp=0, idx.fix=lis1$idx.fix,idx.omega=lis1$idx.omega, idx.res=lis1$idx.res, 
                            tr.fix=lis1$tr.fix))
  }
  jac1 <- (dLL-ll0)/(deps*theta0)
  cat("delta=",deps,": grad=",jac1,"\n")
}

# Still singular...
jacfim <- t(jacprac) %*% jacprac
jacfim
oxboy.fit@results@fim # completely different from FIM saemix 3.1
oxboy.fim@results@fim # completely different from FIM from Hessian

solve(jacfim)
sqrt(1/diag(jacfim))

############
jacky<-computeJacFIMinv.lin(theo.fit)

# Computing a couple of the Jacobian term
lis1 <- initialiseComputationLL(theo.fit)
ll0 <- saemix.LLlin(theta=lis1$theta, modeltype=lis1$modeltype, etype=lis1$etype, f0=lis1$f0, DF0=lis1$DF0, z=lis1$z, Mcov=lis1$Mcov, nfix=lis1$nfix, 
                    nomega=lis1$nomega, nind.obs=lis1$nind.obs, ndat.exp=0, idx.fix=lis1$idx.fix,idx.omega=lis1$idx.omega, idx.res=lis1$idx.res, 
                    tr.fix=lis1$tr.fix) 
# Very different results !!! ??? why ?
pracma::jacobian(saemix.LLlin,lis1$theta, modeltype=lis1$modeltype, etype=lis1$etype, f0=lis1$f0, DF0=lis1$DF0, z=lis1$z, Mcov=lis1$Mcov, nfix=lis1$nfix,  nomega=lis1$nomega, nind.obs=lis1$nind.obs, ndat.exp=0, idx.fix=lis1$idx.fix,idx.omega=lis1$idx.omega, idx.res=lis1$idx.res,  tr.fix=lis1$tr.fix) 
numDeriv::jacobian(saemix.LLlin,lis1$theta, modeltype=lis1$modeltype, etype=lis1$etype, f0=lis1$f0, DF0=lis1$DF0, z=lis1$z, Mcov=lis1$Mcov, nfix=lis1$nfix,  nomega=lis1$nomega, nind.obs=lis1$nind.obs, ndat.exp=0, idx.fix=lis1$idx.fix,idx.omega=lis1$idx.omega, idx.res=lis1$idx.res,  tr.fix=lis1$tr.fix, method.args=list(d=1e-6)) 
numDeriv::jacobian(saemix.LLlin,lis1$theta, modeltype=lis1$modeltype, etype=lis1$etype, f0=lis1$f0, DF0=lis1$DF0, z=lis1$z, Mcov=lis1$Mcov, nfix=lis1$nfix,  nomega=lis1$nomega, nind.obs=lis1$nind.obs, ndat.exp=0, idx.fix=lis1$idx.fix,idx.omega=lis1$idx.omega, idx.res=lis1$idx.res,  tr.fix=lis1$tr.fix, method.args=list(d=.Machine$double.eps^(1/3))) 

# By hand: more consistent than with oxboy... 1e-6, 1e-7, 1e-8 give similar results
theta0<-lis1$theta
for(deps in c(1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)) {
  dLL<-c()
  for(i in 1:length(theta0)) {
    theta1<-theta0
    theta1[i]<-theta1[i]*(1+deps)
    dLL<-c(dLL,saemix.LLlin(theta=theta1, modeltype=lis1$modeltype, etype=lis1$etype, f0=lis1$f0, DF0=lis1$DF0, z=lis1$z, Mcov=lis1$Mcov, nfix=lis1$nfix, 
                            nomega=lis1$nomega, nind.obs=lis1$nind.obs, ndat.exp=0, idx.fix=lis1$idx.fix,idx.omega=lis1$idx.omega, idx.res=lis1$idx.res, 
                            tr.fix=lis1$tr.fix))
  }
  jac1 <- (dLL-ll0)/(deps*theta0)
  cat("delta=",deps,": grad=",jac1,"\n")
}

############################################################################ 
# Checking the expression of Ezi

saemixObject <- theo.fit

