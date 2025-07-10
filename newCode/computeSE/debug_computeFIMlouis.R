library(saemix)
library(numDeriv)

saemixDir <- "/home/eco/work/saemix/saemixextension"
# new linearised FIM by Jacobian and Hessian (comparison)
source(file.path(saemixDir,"newCode","compute_derivFIMlin.R")) 

# Louis' method
source(file.path(saemixDir,"newCode","compute_FIMLouis.R")) 

###################################### Simulated Emax model (N=100, n=3, N-rich)
data(PD1.saemix)

data.pd1<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
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

# Model without covariates
model.pd1<-saemixModel(model=modelemax,description="Emax growth model", 
                    psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,c("E0","Emax","EC50"))), transform.par=c(1,1,1),
                    covariate.model=matrix(c(0,0,0), ncol=3,byrow=TRUE),fixed.estim=c(1,1,1))

# Model with covariates, no IIV on ka and correlation between EC50 and Emax
covmodel1<-diag(c(0,1,1))
covmodel1[2,3]<-covmodel1[3,2]<-1
modelcov.pd1<-saemixModel(model=modelemax,description="Emax growth model", 
                    psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL, c("E0","Emax","EC50"))), transform.par=c(1,1,1),
                    covariance.model = covmodel1,
                    covariate.model=matrix(c(0,0,1), ncol=3,byrow=TRUE),fixed.estim=c(1,1,1))

emx.options<-list(nb.chains=3,seed=765754, nbiter.saemix=c(500,300),save=FALSE,save.graphs=FALSE, displayProgress=FALSE)

fit1<-saemix(model.pd1,data.pd1,emx.options)
fit2<-saemix(modelcov.pd1,data.pd1,emx.options)

# Conditional distribution
nsamp <- 10 
# nsamp<-300
cond.fit1 <- conddist.saemix(fit1, nsamp = nsamp)
cond.fit2 <- conddist.saemix(fit1, nsamp = nsamp)
dim(cond.fit1@results@phi.samp)

# converting to a single array
# phi.try <- matrix(do.call(rbind,list(cond.fit1@results@phi.samp)), ncol=saemixObject@model@nb.parameters) # nope...
cond.fit <- cond.fit2
phi.samp <- NULL
for(isamp in 1:dim(cond.fit@results@phi.samp)[3])
  phi.samp<-rbind(phi.samp, cond.fit@results@phi.samp[,,isamp])

###################################### 
# create args, Dargs, DYF, pres
saemixObject <- fit2
saemix.data<-saemixObject["data"]
saemix.model<-saemixObject["model"]
N<-saemix.data["N"]
nb.parameters<-saemixObject["model"]["nb.parameters"]
# using several Markov chains
chdat<-new(Class="SaemixRepData",data=saemixObject["data"], nb.chains=nsamp)
NM<-chdat["NM"]
IdM<-chdat["dataM"]$IdM
yM<-chdat["dataM"]$yM
XM<-chdat["dataM"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
io<-matrix(data=0,nrow=N,ncol=max(saemixObject["data"]["nind.obs"]))
for(i in 1:N)
  io[i,1:saemixObject["data"]["nind.obs"][i]]<-1
ioM<-do.call(rbind,rep(list(io),nsamp))
ind.ioM <- which(t(ioM)!=0)
DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])
Dargs<-list(transform.par=saemixObject["model"]["transform.par"], structural.model=saemixObject["model"]["model"],IdM=IdM,XM=XM,yM=yM,etype.exp=which(saemixObject["model"]["error.model"] == "exponential"), modeltype=saemixObject["model"]["modeltype"])
pres<-saemixObject["results"]["respar"]
args<-list(ind.ioM=ind.ioM)

ll1 <- compute.LLy(phi.samp,args,Dargs,DYF,pres)
length(ll1)

###################################### 
# Model without covariates
saemixObject<-new(Class="SaemixObject",data.pd1,model.pd1,emx.options)
# Model with covariates
saemixObject<-new(Class="SaemixObject",data.pd1,modelcov.pd1,emx.options)

if(!saemixObject["options"]$warnings) options(warn=-1)
saemix.options<-saemixObject["options"]
saemix.model<-saemixObject["model"]
saemix.data<-saemixObject["data"]
saemix.data@ocov<-saemix.data@ocov[saemix.data@data[,"mdv"]==0,,drop=FALSE]
saemix.data@data<-saemix.data@data[saemix.data@data[,"mdv"]==0,]
saemix.data@ntot.obs<-dim(saemix.data@data)[1]

xinit<- initialiseMainAlgo(saemix.data,saemix.model,saemix.options)

# Elements from the list
saemix.model<-xinit$saemix.model
Uargs<-xinit$Uargs
varList<-xinit$varList
opt<-xinit$opt

saemix.data<-saemixObject["data"]
saemix.model<-saemixObject["model"]
N<-saemix.data["N"]
nb.parameters<-saemixObject["model"]["nb.parameters"]
# using several Markov chains


chdat<-new(Class="SaemixRepData",data=saemixObject["data"], nb.chains=nsamp)
NM<-chdat["NM"]
IdM<-chdat["dataM"]$IdM
yM<-chdat["dataM"]$yM
XM<-chdat["dataM"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
io<-matrix(data=0,nrow=N,ncol=max(saemixObject["data"]["nind.obs"]))
for(i in 1:N)
  io[i,1:saemixObject["data"]["nind.obs"][i]]<-1
ioM<-do.call(rbind,rep(list(io),nsamp))
ind.ioM <- which(t(ioM)!=0)
DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])
Dargs<-list(transform.par=saemixObject["model"]["transform.par"], structural.model=saemixObject["model"]["model"],IdM=IdM,XM=XM,yM=yM,etype.exp=which(saemixObject["model"]["error.model"] == "exponential"), modeltype=saemixObject["model"]["modeltype"])
pres<-saemixObject["results"]["respar"]
args<-list(ind.ioM=ind.ioM)

# From the fitted object
saemixObject <- fit2

# Compute mean.phi (=mu+betaCOV)
COV<-matrix(nrow=saemixObject["data"]["N"],ncol=0)
for(j in 1:saemixObject["model"]["nb.parameters"]) {
  jcov<-which(saemixObject["model"]["betaest.model"][,j]==1)
  aj<-as.matrix(saemixObject["model"]["Mcovariates"][,jcov])
  COV<-cbind(COV,aj)
}
NCOV <- do.call(rbind,rep(list(COV),nsamp))

# Then we can get nsamp*N values of phi_pop (= mu_phi + Sum beta.COV)
mean.phi <- NCOV%*%saemixObject["results"]["MCOV"]
diff.phi <- phi.samp - mean.phi

# Extract theta from results
theta <- c(saemixObject@results@fixed.effects, convert.cov2vec(saemixObject@results@omega, removeZero = TRUE), saemixObject@results@respar[saemixObject@results@indx.res])

# nb of parameters
npar <- saemixObject@model@nb.parameters
nfix <- length(saemixObject@results@betas) # nb of fixed parameters (mu + beta)
nomega <- length(saemixObject@results@indx.omega) # nb of parameters with IIV
nvar <- sum(saemixObject@model@covariance.model[lower.tri(saemixObject@model@covariance.model,diag=TRUE)]) # nb of parameters in Omega (var+covar)
nres<-length(saemixObject@results@indx.res)
ntheta <- nfix + nvar + nres

# Indices
# Elements of covariance matrix (variance refer to diagonal, while off-diagonal terms refer to lower.triangular)
indx.omega <- saemixObject@results@indx.omega
indx.ltri.covar <- which(saemixObject@model@covariance.model[lower.tri(saemixObject@model@covariance.model)]==1) # index of covariances (off-diagonal terms)
indx.utri.covar <- which(saemixObject@model@covariance.model[upper.tri(saemixObject@model@covariance.model)]==1) # index of covariances (off-diagonal terms) # same...if generalisable, would make conversion much simpler
indx.res <- saemixObject@results@indx.res

# additional index mapping covariance terms back to omega
# indx.covar.inomega 
# Keeping track of indices **in theta**
idx.mu <- saemixObject@results@indx.fix # index of mu
idx.beta <- saemixObject@results@indx.cov # index of beta (covariate effects)
idx.omega <- nfix+1:length(indx.omega) # index of variances (diagonal terms) in theta
if(length(1:length(indx.ltri.covar))) idx.covar <- nfix + length(indx.omega) + 1:length(indx.ltri.covar) else idx.covar<-c() # index of covariances (off-diagonal terms) in theta
if(saemixObject@model@modeltype=="structural") idx.res <-nfix+nvar+saemixObject@results@indx.res else indx.res<-c()  # index of residual error parameters
# Elements of mcov to be estimated
idx.mcov <- Uargs$j.covariate

# compute.LLtheta
lldiff2 <- compute.LLtheta(mu=mean.phi[,saemixObject@results@indx.omega,drop=FALSE], 
                           omega=fit2@results@omega[saemixObject@results@indx.omega, saemixObject@results@indx.omega, drop=FALSE], 
                           phi=phi.samp[,saemixObject@results@indx.omega,drop=FALSE]) 

############################### 
# inputs to compute.completeLL(theta, phi.samp, MCOV, NCOV, nfix=nfix, idx.mcov, tr.fix, npar=npar, idx.mu=idx.mu, indx.omega=indx.omega, idx.omega=idx.omega, idx.covar=idx.covar, indx.ltri.covar= indx.ltri.covar)
theta <- c(saemixObject@results@fixed.effects, convert.cov2vec(saemixObject@results@omega, removeZero = TRUE), saemixObject@results@respar[saemixObject@results@indx.res])
MCOV <- saemixObject["results"]["MCOV"]
tr.fix <- saemixObject@model@transform.par
phi<-phi.samp

###############################
# initialise with saemixObject
omega <- saemixObject@results@omega[saemixObject@results@indx.omega, saemixObject@results@indx.omega, drop=FALSE]
mean.phi <- NCOV%*%saemixObject["results"]["MCOV"]
indx.omega <- saemixObject@results@indx.omega

# initialise with theta
if(length(idx.covar)>0) 
  omega <- convert.vec2cov(size=npar, variances=theta[idx.omega], indx.var=indx.omega, covariances = theta[idx.covar ], indx.covar = indx.ltri.covar) else
  omega <- convert.vec2cov(size=indx.omega, variances=theta[idx.omega], indx.omega=indx.omega)

phi <- phi.samp

parfix <- theta[1:nfix]
parfix[idx.mu]<-transpsi(t(parfix[idx.mu]), tr.fix)
MCOV[idx.mcov]<-parfix
mean.phi <- NCOV%*%MCOV

pres<-c(0,0)
pres[indx.res]<-theta[idx.res]

args<-list(ind.ioM=ind.ioM, MCOV=MCOV, NCOV=NCOV, npar=npar, nfix=nfix, idx.mu=idx.mu, idx.omega=idx.omega, idx.res=idx.res, idx.mcov=idx.mcov, idx.covar=idx.covar, indx.omega=indx.omega, indx.res=indx.res, indx.ltri.covar=indx.ltri.covar)

# compute.LLtheta
lltheta <- compute.LLtheta(mu=mean.phi[,indx.omega,drop=FALSE], 
                           omega=omega[indx.omega, indx.omega, drop=FALSE], 
                           phi=phi[,indx.omega,drop=FALSE]) 
# compute.LLy
lly <- compute.LLy(phi,args,Dargs,DYF,pres)

# compute.completeLL
LLcomp <- compute.completeLL(theta, phi, args, Dargs, DYF)
summary(LLcomp-(lly+lltheta)))

# premier terme à 0 :-/ 
## (paramètre sans IIV, passe dans MCOV mais ensuite non pris en compte dans compute.completeLL et n'apparaît pas dans compute.LLy donc logique...)
jacobLL<-try(numDeriv::jacobian(compute.completeLL, theta, phi=phi, args=args, Dargs=Dargs, DYF=DYF, method.args=list(d=1e-6)))

# and we can't get the Hessian because we need a scalar valued function and here we have a vector...
# much too costly to do for each subject so probably the end of the road... ?
# hessLL<-try(numDeriv::hessian(compute.completeLL, theta, phi=phi, args=args, Dargs=Dargs, DYF=DYF, method.args=list(d=1e-6)))


############################### Testing compute.LLtheta
# compute LLtheta, on parameters with IIV
omega <- fit2@results@omega[saemixObject@results@indx.omega, saemixObject@results@indx.omega, drop=FALSE]
invOmega <- solve(omega)
lldiff <- c()
for(i in 1:dim(diff.phi)[1])
  lldiff<-c(lldiff,
            diff.phi[i,saemixObject@results@indx.omega,drop=FALSE] %*% invOmega %*% t(diff.phi[i,saemixObject@results@indx.omega,drop=FALSE]))
lldiff <- (lldiff-log(det(omega)))/2
lldiff2 <- compute.LLtheta(mu=mean.phi[,saemixObject@results@indx.omega,drop=FALSE], 
                           omega=fit2@results@omega[saemixObject@results@indx.omega, saemixObject@results@indx.omega, drop=FALSE], 
                           phi=phi.samp[,saemixObject@results@indx.omega,drop=FALSE]) 
summary(lldiff-lldiff2)
# check if matrix way to do this and if faster ?

################################
# En attente

compute.phipop <- function() {
  pfix<-matrix(data=0,nrow=1,ncol=nb.parameters)
  LCOV<-MCOV<-matrix(data=0,nrow=nb.betas,ncol=nb.parameters)
  j1<-1
  COV<-matrix(nrow=dim(Mcovariates)[1],ncol=0)
  mean.phi<-matrix(data=0,nrow=saemix.newdata["N"],ncol=nb.parameters)
  fixed.ini<-saemixObject["model"]["betaest.model"]*0
  fixed.ini[saemixObject["model"]["betaest.model"]==1]<-pop.par
  for(j in 1:nb.parameters) {
    jcov<-which(saemixObject["model"]["betaest.model"][,j]==1)
    lambdaj<-fixed.ini[jcov,j]
    aj<-as.matrix(Mcovariates[,jcov])
    COV<-cbind(COV,aj)
    nlj<-length(lambdaj)
    j2<-j1+nlj-1
    LCOV[j1:j2,j]<-matrix(data=1,nrow=nlj,ncol=1)
    j1<-j2+1
    if(length(jcov)<=1) mean.phi[,j]<-aj*lambdaj else mean.phi[,j]<-aj%*%lambdaj
    pfix[j]<-length(lambdaj)
  }
}

################################
# Some matrix algebra
mat1<-matrix(1:16, ncol=4)
lower.tri(mat1, diag=TRUE)
mat1[lower.tri(mat1)]

phi1 <- phi.samp[1:10,]
omega <- diag(c(0.1,0.2,0.5))

# Function to convert a symmetric matrix to a vector
convert.cov2vec <- function(xmat, removeZero = FALSE) {
  # Convert a symmetric non-negative matrix into a vector (no check)
  if(removeZero) { # keep only non-zero lines/columns
    x <- which(abs(mydiag(xmat))>.Machine$double.xmin)
    xmat <- xmat[x, x, drop=FALSE]
  }
  # return first the diagonal then the off-diagonal elements
  xoff<-xmat[lower.tri(xmat)]
  if(removeZero) c(mydiag(xmat), xoff[abs(xoff)>.Machine$double.xmin]) else c(mydiag(xmat), xoff)
}
convert.cov2vec(saemixObject@model@covariance.model)
convert.cov2vec(saemixObject@model@covariance.model, removeZero = TRUE)

# Function to create a symmetric matrix with specified variances and covariance terms
convert.vec2cov <- function(size, variances, idx.var, covariances=NULL, idx.covar=NULL) {
  # input: 
  ## size: matrix size
  ## variances: vector of elements on the diagonal
  ## idx.var: indexes of the elements
  ## covariances: vector of off-diagonal elements
  ## idx.covar: indexes of the off-diagonal elements
  xmat<-matrix(data=0, nrow=size, ncol= size)
  diag(xmat)[idx.var] <- variances
  if(!is.null(covariances)) {
    xmat[lower.tri(xmat)][idx.covar]<-covariances
    xmat[upper.tri(xmat)][idx.covar]<-covariances
  }
  return(xmat)
} 
xmat1 <- convert.vec2cov(3,c(1:3), c(1:3))
xmat1
xmat2 <- convert.vec2cov(3,c(1:3), c(1:3), c(4:5),c(1,3))
xmat2

xmat3 <- convert.vec2cov(3,c(0.5,0.2), c(2:3), c(0.1),c(3))
xmat3


# convert.vec2cov <- function(variances, size=NULL, idx.var=NULL, covariances=NULL, idx.covar=NULL) {
#   if(is.null(idx.var) & is.null(size)) { # assume we want a matrix with variances as a diagonal
#     mat<-mydiag(variances)
#   } else {
#     if(!is.null(size)) {
#       xmat<-matrix(data=0, nrow=size, ncol= size)
#       if(!is.null(idx.var)) mydiag(xmat)[idx.var]<-variances
#     }
#     
#   }
# }  
