########################################################### Setup
saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
progDirExt<-file.path(saemixDir,"Rext")
source(file.path(progDir,"aaa_generics.R"))

# Data class
source(file.path(progDirExt,"SaemixOutcome.R"))
source(file.path(progDirExt,"SaemixOutcome-methods.R"))
source(file.path(progDirExt,"SaemixData.R"))
source(file.path(progDirExt,"SaemixCovariateModel.R"))
source(file.path(progDirExt,"SaemixCovariate.R"))

# Model class
source(file.path(progDirExt,"SaemixVarLevel.R"))
source(file.path(progDirExt,"SaemixParameter.R"))
source(file.path(progDirExt,"SaemixParameter-methods.R"))
source(file.path(progDirExt,"SaemixModel.R"))
source(file.path(progDirExt,"SaemixModel-methods.R"))

########################################################### Covariates
# Creating an SaemixData object with PK/PD data by hand - Responses

# Generating covariates assuming wt and treatment varies with occasion
nsuj<-10
nocc<-3 # same nb of occasions for each subject
cov.gender<-rbinom(nsuj,size=1,p=0.5)
age<- round(rnorm(nsuj, mean=50, sd=5))
wt<- rnorm(nsuj,mean=60, sd=4)
wt[cov.gender==1]<-wt[cov.gender==1]*1.2
wt<-round(wt, digits=1)
cyp <- rbinom(nsuj,size=1,p=0.5)
wt2<-rep(wt, each=3)
wt2<-wt2+round(runif(length(wt2),-5,5))
covmat.id <- data.frame(id=1:nsuj, sex=cov.gender, lage=log(age/50), cyp=cyp)
covmat.occ <- data.frame(id=rep(1:nsuj, each=3), occ=rep(c(1:3), nsuj), lwt=log(wt2/60), trt=rep(LETTERS[1:3], nsuj))
# transform categorical covariate
xmat<-covmat.occ
name.var <- "trt"
refcat<-"A"
xfac <- factor(xmat[,name.var])
ncat<-length(unique(xfac))
dumcov<-xmat[,name.var,drop=FALSE]
for(ilev in levels(xfac)) {
  if(as.character(ilev)!=refcat) {
    idum<-as.integer(xfac==ilev)
    name.dummy<-paste0(name.var,ilev)
    dumcov[,name.dummy]<-idum
  }
}
# idx2 <- which(colSums(dumcov[,-c(1)])==0) # test if some categories not represented ? (or before)

# replacing columns in covmat.occ
i1<-which(colnames(covmat.occ)==name.var)
covmat.occ<-covmat.occ[,-c(i1)]
covmat.occ<-cbind(covmat.occ, dumcov[,-c(1),drop=FALSE])

# startsWith(colnames(covmat.occ),"trt")

## TODO: will need to build covmat.id and covmat.occ automatically from the data
### check if several variability levels
### ie create a covariate matrix per id and a covariate matrix per id x occ
### if same value of cov for the different occasions (same size removing duplicates ?) then per id, otherwise per occ and keep only first value

# List of covariate matrices
covmat<-list(id=covmat.id, occ=covmat.occ)

########################################################### Parameter model
# Mean value for population parameters
pkpd.psi0<-c(ka=1, vd=20, cl=5, ic50=3)

# Model parameters
lpar <- list(ka=saemixPar(mu.start=pkpd.psi0[1], omega.start=0.2, covariate=c(sex=binCov(beta=0.3))),
             cl=saemixPar(mu.start=pkpd.psi0[2], omega.start=c(0.3,0.4), covariate=c(sex=binCov(beta=0.2), lwt=contCov(name="wt",transform=log, beta=0.75, beta.fix=1), lage=contCov(name="age", transform=log, beta=-0.5))),
             vd=saemixPar(mu.start=pkpd.psi0[3], omega.start=c(0.3,0.4), covariate=c(lwt=contCov(name="wt",transform=log, beta=1, beta.fix=1), cyp=binCov(beta=0.5)), rho.param=c("cl"), rho=0.5),
             ic50=saemixPar(mu.start=pkpd.psi0[4], omega.start=0.4, covariate=c(lage=contCov(name="age", transform=log, beta=0.5), trt=catCov(beta=c(0.1,0.3)))))
npar<-length(lpar)


transform.par<-rep(1,npar) # to do automagically

# Covariance models - change to add moore than one level of variability
omega1<-diag(rep(0.3,4))
omega1[1:3,1:3]<-vec2mat(c(0.3,0.8*.3,0,0.3,0,0.5))
omega2<-diag(rep(0,4))
omega2[1:3,1:3]<-vec2mat(c(0.4,0.2*.4,0,0.4,0,0))
var.list<-list(saemixVarModel(omega=omega1, variable="id"), saemixVarModel(omega=omega2, variable="occ", verbose=T))

varlevel<-c()
for(i in 1:length(var.list)) varlevel<-c(varlevel, var.list[[i]]@variable)

# nb of parameters 
## need to transform covariates before (when matching model+data)
## need to cross-check covariate type and covariate
nvarlevel<-length(varlevel)
nb.betas<-rep(0, nvarlevel)
notfound<-c()
for(ipar in 1:npar) {
  if(length(lpar[[ipar]]@covariate)>0) {
    for(icov in 1:length(lpar[[ipar]]@covariate)) {
      name.var<-names(lpar[[ipar]]@covariate)[icov]
      inotfound<-0
      if(lpar[[ipar]]@covariate[[icov]]@type=="categorical" && (length(lpar[[ipar]]@covariate[[icov]]@covariate.transform@ncat)==0 || lpar[[ipar]]@covariate[[icov]]@covariate.transform@ncat>2)) {
        for(ivar in 1:nvarlevel) {
          idx2<-which(startsWith(colnames(covmat[[ivar]]),name.var))
          if(length(idx2)>0) {
            nb.betas[ivar]<-nb.betas[ivar]+length(idx2)
            inotfound<-1
          }
        }
      } else {
        for(ivar in 1:nvarlevel) {
          if(name.var %in% colnames(covmat[[ivar]])) {
            nb.betas[ivar]<-nb.betas[ivar]+1
            inotfound<-1
          }
        }
      }
      if(inotfound==0) notfound<-c(notfound, name.var)
    }
  }
}
if(length(notfound)>0) cat("Covariates",notfound,"not found in the dataset\n")
# Covariates not found in the data should be removed from the model when creating covariate.model, with an error message

# List of covariate models associated with each variability level (to create automatically, including covariate transformations and dummy covariates for categories)
covmodel1 <- matrix(data=0, nrow=3, ncol=npar)
covmodel2 <- matrix(data=0, nrow=3, ncol=npar)
colnames(covmodel1)<-colnames(covmodel2)<-names(lpar)
rownames(covmodel1)<-c("sex","lage","cyp")
covmodel1[1,1]<-covmodel1[1,2]<-1
covmodel1[2,2]<-covmodel1[2,4]<-1
covmodel1[3,3]<-1
rownames(covmodel2)<-c("lwt","trtB","trtC")
covmodel2[1,3]<-covmodel2[1,2]<-1
covmodel2[2,4]<-covmodel2[3,4]<-1
covariate.model<-list(id=covmodel1, occ=covmodel2)
# To automate also
betas.start<-betas.fix<-name.covbetas<-c()
for(ipar in 1:npar) {
  if(length(lpar[[ipar]]@covariate)>0) {
    for(icov in 1:length(lpar[[ipar]]@covariate)) {
      betas.start<-c(betas.start,lpar[[ipar]]@covariate[[icov]]@beta)
      betas.fix<-c(betas.fix,lpar[[ipar]]@covariate[[icov]]@beta.fix)
      xnam<-rep(names(lpar[[ipar]]@covariate)[icov],length(lpar[[ipar]]@covariate[[icov]]@beta))
      name.covbetas<-c(name.covbetas, xnam)
    }
  }
}
name.covbetas[name.covbetas=="trt"]<-paste0("trt",c("B","C"))
# creating lists
list1<-list2<-list3<-vector(mode="list", length=nvarlevel)
for(ivar in 1:nvarlevel) {
  id1<-which(!is.na(match(name.covbetas, rownames(covariate.model[[ivar]]))))
  list1[[ivar]]<-name.covbetas[id1]
  list2[[ivar]]<-betas.start[id1]
  list3[[ivar]]<-betas.fix[id1]
}
name.covbetas<-list1
betas.start<-list2
betas.fix<-list3

# Complete mu+covariate model:
rbind(mu=rep(1,npar),covariate.model[[1]], matrix(nrow=0,ncol=4), covariate.model[[2]])
# TODO: if fixed period effects, need to change the 0 matrix...

######################################################################## Individual parameters
# Generating subject level individual parameters
## By hand
mean.fixedphi <- transpsi(matrix(pkpd.psi0,nrow=1), transform.par)
phi.id<-do.call(rbind, rep(list(mean.fixedphi), nsuj))
phi.id[,1]<-phi.id[,1] + lpar$ka@covariate$sex@beta * covmat.id$sex
phi.id[,2]<-phi.id[,2] + lpar$cl@covariate$sex@beta * covmat.id$sex + lpar$cl@covariate$lage@beta * covmat.id$lage
phi.id[,3]<-phi.id[,3] + lpar$vd@covariate$cyp@beta * covmat.id$cyp
phi.id[,4]<-phi.id[,4] + lpar$ic50@covariate$lage@beta * covmat.id$lage

# Generating occasion level parameters
## by hand
phi.occ <- matrix(data=0, nrow = nsuj*nocc, ncol=length(pkpd.psi0))
phi.occ[,2]<-phi.occ[,2] + lpar$cl@covariate$lwt@beta * covmat.occ$lwt
phi.occ[,3]<-phi.occ[,3] + lpar$vd@covariate$lwt@beta * covmat.occ$lwt  
phi.occ[,4]<-phi.occ[,4] + lpar$ic50@covariate$trt@beta[1] * covmat.occ$trtB + lpar$ic50@covariate$trt@beta[1] * covmat.occ$trtC

phi<-list(phi.id, phi.occ)

# 'population' individual parameters (eg mu + sum_k beta_ki for each level)
phimat <- cbind(repmatrix(phi.id, times=3), phi.occ)
colnames(phimat) <- paste(rep(names(lpar),2),rep(c("id","occ"),each=4),sep=".")

# Generate individual phi_i and phi_ik adding etas
chol.omega1 <- chol(omega1)
eta1<-matrix(rnorm(nsuj*4),ncol=4)%*%chol.omega1
phi.id <- phi.id + eta1
phi1 <- repmatrix(phi.id[,1:4], times=3)
chol.omega2 <- chol(omega2[1:2,1:2])
eta2<-matrix(rnorm(dim(phi.occ)[1]*2),ncol=2)%*%chol.omega2
phi.occ[,1:2] <- phi.occ[,1:2] + eta2

phi.ik<-phi1 + phi.occ[,1:4]
psi.ik <- transphi(phi.ik, transform.par)
summary(psi.ik)

## put everything in a structure (SaemixIndividualModel ?)
## see what else is needed for computations - especially optimising fixed effects ? parameters without IIV ?

######################################################################## Observations
# Generating occasion-level observations

# Generating associated PK/PD data
model1cptdirect<-function(psi,id,xidep) { 
  tim<-xidep[,1]
  dose<-xidep[,2]
  ytype<-xidep$ytype
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  ic50<-psi[id, 4]
  k<-CL/V
  ypk<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  ypd<-100*(1-ypk/(ypk+ic50))
  ypk[ytype==2]<-ypd[ytype==2]
  return(ypk)
}

# Outcomes
xout<-list(conc=saemixOutcome(unit="mg/L", model="combined2", start=c(1, 0.2)), 
           effect=saemixOutcome(unit="%", start=1))

# PK/PD design
sampleTimes<-c(0,1,2,4,8)
xidep1<-data.frame(tim=rep(sampleTimes,each=2), dose=100, ytype=rep(c(1,2), length(sampleTimes)))
xidep<-do.call(rbind, rep(list(xidep1), nsuj*nocc))
id.suj <- rep(1:nsuj, each=length(sampleTimes)*2*nocc)
id.occ <- rep(1:(nsuj*nocc), each=length(sampleTimes)*2)
fpred <- model1cptdirect(psi.ik,id.occ,xidep)

ypred<-fpred
for(iout in 1:length(xout)) {
  idx1<-which(xidep$ytype==iout)
  if(xout[[iout]]@type=="continuous") {
    if(xout[[iout]]@error.model=="exponential") fpred[idx1]<-log(cutoff(fpred[idx1]))
    gpred<-xout[[iout]]@error.function(fpred[idx1], xout[[iout]]@error.parameters)
    ypred[idx1] <- ypred[idx1] + gpred
  }
}

par(mfrow=c(1,2))
for(iout in 1:length(xout)) {
  plot(xidep[xidep$ytype==iout,1], ypred[xidep$ytype==iout], pch=20, xlab="Time", ylab=names(xout)[iout])
  for(isuj in unique(id.occ)) 
    lines(xidep[xidep$ytype==iout & id.occ==isuj,1], ypred[xidep$ytype==iout  & id.occ==isuj])
}

######################################################################## Dataset
## covariates already transformed
cov1<-repmatrix(covmat.id[,2:4], times=rep(length(sampleTimes)*2*nocc, nsuj))
colnames(cov1)<-colnames(covmat.id)[2:4]
cov2<-repmatrix(covmat.occ[,3:5], times=rep(length(sampleTimes)*2, nsuj*nocc))
colnames(cov2)<-colnames(covmat.occ)[3:5]
pkpd.dat <- data.frame(id=id.suj, xidep, y=ypred, occ=rep(1:3, each=length(sampleTimes)*2), cov1, cov2)

# matching the lines in pkpd.dat with the lines in phimat
var.match <- cbind(id.suj, id.occ)

########################################################################
