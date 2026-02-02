# Loading saemix
library(saemix)

# Libraries needed to compute the FIM by AGQ
library(R6)
library(pracma)
library(compiler)
library(statmod)
library(Matrix)
library(randtoolbox)

# Folders
workDir<-"/home/eco/work/saemix/discreteEval"
setwd(workDir)

saemixDir<-"/home/eco/work/saemix/saemixextension"
dirAGQ <- file.path(saemixDir, "fimAGQ")

scenario <- "binaryIIV"
dir.create(file.path(workDir, scenario))
datDir <- file.path(workDir, scenario, "data")
dir.create(datDir)
resDir <- file.path(workDir, scenario, "results")
dir.create(resDir)

nsim<-200 # Number of simulations

set.seed(711223817)
zeseeds <- trunc(runif(nsim)*1e8)

###################################################### Setup
# Functions
## Model
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

## Function to simulate from model
simulBinary<-function(psi,id,xidep) {
  tim<-xidep[,1]
#  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-1/(1+exp(-logit))
  ysim<-rbinom(length(tim),size=1, prob=pevent)
  return(ysim)
}


# Settings
## Design
months<-c(0, 1, 2, 3, 5.5, 8, 11)
par.fix<-c(-1.71, -0.39, -0.15)
par.om <- c(1, 0.2) # around 50% variability for both
nsuj<-274
ndat<-length(months)
xidep <- data.frame(time=rep(months, nsuj))
xidep$y<-rep(0,dim(xidep)[1])
id1 <- rep(1:nsuj, each=ndat)
covtreat <- rep(0, length(id1))
covtreat[id1>(nsuj/2)]<-1
visit <- rep(1:7, nsuj)

## saemix model for fits
saemix.model.true<-saemixModel(model=binary.model,description="Binary model",modeltype="likelihood",
                          psi0=matrix(c(par.fix[1:2],0,par.fix[3]),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                          transform.par=c(0,0), covariate.model=c(0,1), omega.init=diag(c(par.om**2)))

saemix.model.pop<-saemixModel(model=binary.model,description="Binary model",modeltype="likelihood",
                               psi0=matrix(c(-0.5,-0.19,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                               transform.par=c(0,0), covariate.model=c(0,1), omega.init=diag(c(1,1)))

saemix.model.far<-saemixModel(model=binary.model,description="Binary model",modeltype="likelihood",
                              psi0=matrix(c(0,0,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                              transform.par=c(0,0), covariate.model=c(0,1),omega.init=diag(c(4,0.5)))
## options
saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE)

# Prepare result files
namResFiletrue <- file.path(resDir, paste0("truefit_",scenario,".tab"))
namResFilepop <- file.path(resDir, paste0("popfit_",scenario,".tab"))
namResFilefar <- file.path(resDir, paste0("farfit_",scenario,".tab"))
namcol <- c("irep", "theta1", "theta2", "beta", "omega1", "omega2", "LL.IS")
for(ifile in c(namResFiletrue, namResFilepop, namResFilefar)) {
  write(namcol, file=ifile, ncolumns=length(namcol))
}

###################################################### Loop on simulations

# Loop on simulations
isim<-1
for(isim in 1:nsim) {
  namDatFile <- file.path(datDir, paste0("data_",scenario,"_sim",isim,".tab"))
  
  ## Simulate data
  set.seed(zeseeds[isim])
  psi1<-data.frame(th1=rnorm(nsuj, mean=par.fix[1], sd=par.om[1]), th2=rnorm(nsuj, mean=par.fix[2], sd=par.om[2]))
  psi1[(1+nsuj/2):nsuj,2] <- psi1[(1+nsuj/2):nsuj,2]+par.fix[3]
  
  ysim <- simulBinary(psi1, id1, xidep)
  
  simdat<-cbind(id=id1, xidep, treatment=covtreat, visit=visit)
  simdat$y <- ysim
  
  ## Save data to disk
  write.table(simdat, file=namDatFile, quote=F, row.names=F)
  
  # Fit
  ## Data
  smxdat<-saemixData(name.data=simdat,name.group=c("id"),name.predictors=c("time","y"), name.response="y",
                     name.covariates=c("treatment"),name.X=c("time"))
  
  ## Fitting with 3 different settings
  smxfit.true<-saemix(saemix.model.true,smxdat,saemix.options)
  
  smxfit.pop<-saemix(saemix.model.pop,smxdat,saemix.options)
  
  smxfit.far<-saemix(saemix.model.far,smxdat,saemix.options)
  
  ## Save estimates to disk
  l1<-c(isim, smxfit.true@results@fixed.effects, sqrt(diag(smxfit.true@results@omega)), smxfit.true@results@ll.is)
  l2<-c(isim, smxfit.pop@results@fixed.effects, sqrt(diag(smxfit.pop@results@omega)), smxfit.pop@results@ll.is)
  l3<-c(isim, smxfit.far@results@fixed.effects, sqrt(diag(smxfit.far@results@omega)), smxfit.far@results@ll.is)
  
  write(l1, file=namResFiletrue, ncolumns=length(namcol),append=TRUE)
  write(l2, file=namResFilepop, ncolumns=length(namcol),append=TRUE)
  write(l3, file=namResFilefar, ncolumns=length(namcol),append=TRUE)
}

######################################################  Predicted FIM by MC/AGQ (TODO)
saemix.fit <- smxfit.true

source(file.path(dirAGQ,"default_settings.R"))
source(file.path(dirAGQ,"helper_functions.R"))
source(file.path(dirAGQ,"integration.R"))
source(file.path(dirAGQ,"model.R"))

# Setting up binomial model
# define longitudinal binary model
model <- Model$new(
  parameter_function = function(mu, b) list(base=mu[1]+b[1], slp=mu[2]+b[2], eff=mu[3]),
  log_likelihood_function = function(y, design, base, slp, eff) {
    p <- mapply(function(time, trt) 1/(1+exp(-(base+(slp+eff*trt)*time))), design$time, design$trt)
    log(ifelse(y==1,p,1-p))
  }, 
  simulation_function = function(design, base, slp, eff) {
    p <- mapply(function(time, trt) 1/(1+exp(-(base+(slp+eff*trt)*time))), design$time, design$trt)
    as.numeric(p > runif(nrow(design)))
  },
  inverse_simulation_function = function(design, urand, base, slp, eff) {
    if(is.null(urand)) return(seq_along(design$time))
    p <- mapply(function(time, trt) 1/(1+exp(-(base+(slp+eff*trt)*time))), design$time, design$trt)
    qbinom(urand, 1, prob = p)
  },
  mu = saemix.fit@results@fixed.effects,
  omega = saemix.fit@results@omega)

# define settings (agq with 3 grid points, quasi random monte-carlo and 500 samples)
settings <- defaults.agq(gq.quad_points = 3,  y_integration.method = "qrmc", y_integration.n_samples = 500, seed =1234)


#### Design (built-in data set)
# 2 groups of 137 subjects, with respectively risk=0 and 1
# same 7 times for all subjects c(0, 1, 2, 3, 5.5, 8, 11)
# create design as data frame for 1 subject, risk=0
design <- data.frame(time=c(0, 1, 2, 3, 5.5, 8, 11), trt=0)
# use only half of the samples for each group
settings$y_integration.n_samples <- ceil(settings$y_integration.n_samples/2)
# calculate fim for no risk group
fim_trt0 <- calc_fim(model, design, settings)
# calculate fim for risk group => update design with risk=1
fim_trt1 <- calc_fim(model, transform(design, trt=1), settings)
# add fims together (500 subjects for each group)
fim <-274*(fim_trt0 + fim_trt1)/2
print(fim)
# calculate rse
rse <- calc_rse(model, fim)
print(rse)

est.se<-sqrt(diag(solve(fim)))
df <- data.frame(param=c(model$mu,diag(model$omega)),se=est.se)
df$rse <- abs(df$se/df$param*100)

write.table(df, file=file.path(resDir, "seFIMAGQ_true.res"), quote=F, row.names=TRUE, col.names=TRUE)

