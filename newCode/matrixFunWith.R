# Creating a matrix - values enter in columns
xmat<-matrix(1:20,ncol=5)

# Indices are the same as the matrix values here
summary(xmat[1:18]-1:18)

# Changing an index
xmat1<-xmat
xmat1[14]<-0
xmat1

########## Var-cov matrix
# Diagonal square matrix
xvar<-diag(c(1:5)*2)
colnames(xvar)<-rownames(xvar)<-paste0("par",1:dim(xvar)[1])
xvar1<-xvar

# Symmetric matrix with covariances
xvar<-diag(c(1:5)*2)
xvar[1,2]<-xvar[2,1]<-0.5
xvar[1,3]<-xvar[3,1]<-0.6
xvar[3,2]<-xvar[2,3]<-0.4
colnames(xvar)<-rownames(xvar)<-paste0("par",1:dim(xvar)[1])
xvar2<-xvar

# One IIV missing
xvar3<-xvar2
xvar3[4,4]<-0

# Full symmetric matrix with covariances
xvar<-diag(c(1:5)*2)
xvar[lower.tri(xvar)]<-seq(0.1,1, 0.1) # fill in the lower tri covariance terms
t(xvar)-diag(diag(xvar)) 
xvar<-xvar+ (t(xvar)-diag(diag(xvar))) # to symmetrise the matrix based on its lower triangular form
colnames(xvar)<-rownames(xvar)<-paste0("par",1:dim(xvar)[1])
xvar4<-xvar 

########## Indices
# Associated model
xvar<-xvar4
xvar.mod<-(ifelse(xvar>0,1,0)) # "model" (0/1)

# Variance
xvar.modvar<-(upper.tri(xvar,diag = TRUE)-upper.tri(xvar))*xvar.mod
indx.var<-which(xvar.modvar>0)
param.var <- xvar[indx.var]
names.param.var <- paste0("omega2.",colnames(xvar))

# Covariances
xvar.modcovar<-upper.tri(xvar)*xvar.mod + lower.tri(xvar)*xvar.mod
indx.covar<-which(xvar.modcovar>0)
indx.covar.low<-which(lower.tri(xvar)*xvar.mod>0)
indx.covar.sup <- setdiff(indx.covar, indx.covar.low) # pb: will not be filled in the same order, so useless
param.covar <- xvar[indx.covar.low]
names.param.covar <- c()
for(j in 1:(dim(xvar)[2]-1)) {
  for(i in (j+1):dim(xvar)[1])
    names.param.covar <- c(names.param.covar,paste0("cov.",colnames(xvar)[j],".", rownames(xvar)[i]))
}
paste(names.param.covar,"=",param.covar)
# Variance matrix in index form
indx.varcov <- c(indx.var, indx.covar.low)
param.varcov <- c(param.var, param.covar)
names.param.varcov <- c(names.param.var, names.param.covar)

########## Covariate model matrix
xmubeta.mod <- rbind(rep(1,4), c(1,0,1,0), c(0,1,1,0),c(0,0,1,1), c(0,0,0,1))
rownames(xmubeta.mod) <- c("mu", "sex", "wt", "age","cyp")
colnames(xmubeta.mod) <- c("ka", "V", "CL", "IC50")
xmubeta <- xmubeta.mod

xmu.mod<-xmubeta.mod*0
xmu.mod[1,]<-xmubeta.mod[1,]
indx.mu <- which(xmu.mod>0)
indx.fixpar <- which(xmubeta.mod>0)
indx.beta <- setdiff(indx.fixpar, indx.mu)

xmubeta[indx.mu] <- c(0.5, 20, 5, 100)
xmubeta[indx.beta] <- c(0.3, 1, 0.2, 0.75, -0.5, 0.5, 0.8)
fixedpar <- xmubeta[xmubeta.mod==1]
names.fixedpar <-c()
for(j in 1:dim(xmubeta)[2]) {
  for(i in 1:dim(xmubeta)[1])
    if(xmubeta.mod[i,j]==1) {
      if(i==1) xnam<-paste0("mu.",colnames(xmubeta.mod)[j]) else 
        xnam<-paste0("beta.",colnames(xmubeta.mod)[j],",",rownames(xmubeta.mod)[i])
      names.fixedpar <-c(names.fixedpar, xnam)
    }
}
paste(names.fixedpar,"=",fixedpar)

########## MCOV/LCOV
nb.parameters <- dim(xmubeta.mod)[2]
nb.fixedpar <- length(fixedpar)
LCOV<-matrix(data=0,nrow=nb.fixedpar,ncol=nb.parameters)
colnames(LCOV)<-colnames(xmubeta.mod)
rownames(LCOV)<-names.fixedpar
j1<-1
mean.phi<-matrix(data=0,nrow=nsuj,ncol=nb.parameters) # population parameters for subject i = mu + sum_q beta_q cov_i,q
for(jpar in 1:nb.parameters) {
  npar1 <- sum(xmubeta.mod[,jpar])
  LCOV[j1:(j1+npar1-1),jpar]<-1
  j1<-j1+npar1
}
MCOV<-LCOV
MCOV[LCOV==1]<-fixedpar

########## Design matrix for covariates
nsuj<-10
cov.gender<-rbinom(nsuj,size=1,p=0.5)
age<- round(rnorm(nsuj, mean=50, sd=5))
wt<- rnorm(nsuj,mean=60, sd=4)
wt[cov.gender==1]<-wt[cov.gender==1]*1.2
wt<-round(wt, digits=1)
cyp <- rbinom(nsuj,size=1,p=0.5)
covmat <- data.frame(id=1:nsuj, sex=cov.gender, lwt=log(wt/60), lage=log(age/50), cyp=cyp)

design.cov <- cbind(rep(1,nsuj), covmat$sex, rep(1,nsuj), covmat$lwt, rep(1,nsuj), covmat$sex, covmat$lwt, covmat$lage, rep(1,nsuj), covmat$lage, covmat$cyp)

mean.phi <- design.cov %*% MCOV

# population value for phi


