# Libraries and functions
library(saemix)

# Trace
tr <- function(xmat)
  sum(diag(xmat))

# Linear model
linmod <- function(t,a,b) {
  a+b*t
}
linmod.smx<-function(psi,id,xidep) { 
  tim<-xidep[,1]  
  a<-psi[id,1]
  b<-psi[id,2]
  ypred <- a+b*tim
  return(ypred)
}

##########################################################################################
# Generating data from linear mixed effect model

param<-c(15, 0.3)
omega <- c(0.5,0.1) # 30% variability
sigma <- 0.1

nsuj<-50
xtim <- c(0,1, 2, 5, 10, 24)
ipar <- data.frame(a=rnorm(nsuj, mean=param[1], sd=omega[1]), b=rnorm(nsuj, mean=param[2], sd=omega[2]))
xtab<-NULL
for(isuj in 1:nsuj) {
  ypred<-linmod(xtim, ipar[isuj,1], ipar[isuj, 2])+rnorm(length(xtim),sd=sigma)
  xtab<-rbind(xtab,
              data.frame(id=isuj, x=xtim, y=ypred))
}
plot(xtab$x, xtab$y)
for(isuj in 1:nsuj)
  lines(xtab$x[xtab$id==isuj],xtab$y[xtab$id==isuj])

# Fit by SAEM
xdat <- saemixData(xtab, name.group="id", name.predictors = "x", name.response = "y")
xmodel <- saemixModel(model=linmod.smx,modeltype="structural",  description="Linear model", 
                     psi0=matrix(c(15,1),ncol=2, byrow=TRUE, dimnames=list(NULL, c("a","b"))))
#                                    omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")

yfit <- saemix(xmodel, xdat)
yfit@results@fim
invFIM <- solve(yfit@results@fim)
sqrt(diag(invFIM))

# Individual Fisher Information Matrix for the LME (Wand 2007)
# in our model X=Z=matrix with first column=1 and second colum=xtim
Xmat <- matrix(c(rep(1, length(xtim)),xtim),ncol=2)
Xmat %*% param # should be equal to a+b*t
summary(c(Xmat %*% param)-(param[1]+param[2]*xtim))

# Xmat*beta + Xmat*u = Xmat*EBE

# All subjects have the same Vi (?) since the design is the same and the individual random effects don't intervene
Vi <- Xmat %*% yfit@results@omega %*% t(Xmat) + diag(rep(yfit@results@respar[1],length(xtim)))
invVi<-solve(Vi)

blocA <- nsuj*(t(Xmat) %*% invVi %*% Xmat)
lxmat <- list(Xmat, Xmat, diag(rep(1,length(xtim))))
lxmat[[1]][,2]<-0
lxmat[[2]][,1]<-0

blocB<-matrix(data=0, nrow=3, ncol=3)
for(i in 1:3) {
  for(j in 1:3) {
    mymat <- t(lxmat[[i]]) %*% invVi %*% lxmat[[j]]
    blocB[i,j] <- 0.5*tr(mymat %*% t(mymat))
  }
}
blocB<-blocB*nsuj

# Comparing FIM by saemix and by the exact formula for the fixed effects => very close
yfit@results@fim[1:2,1:2]
blocA
sqrt(diag(solve(blocA)))
sqrt(diag(solve(yfit@results@fim[1:2,1:2])))

# Comparing FIM by saemix and by the exact formula for the variances => very close
yfit@results@fim[3:4,3:4]
blocB[1:2,1:2]
sqrt(diag(solve(blocB)))[1:2]
sqrt(diag(solve(yfit@results@fim)))[3:4]

# The formula computes the SE for sigma2 while we compute it for sigma => use delta method to get var(sigma2) from var(sigma) => somewhat different (factor 5)
sqrt(diag(solve(yfit@results@fim)))[5]
varsig2<-2*yfit@results@respar[1]*diag(solve(yfit@results@fim))[5]
sqrt(varsig2)
sqrt(diag(solve(blocB)))[3]

# but... yields a very high SE (0.00891951) compared to the sigma2=0.007959306 (RSE=112%), which doesn't seem right for such a simple model
seSig2 <- sqrt(diag(solve(blocB))[3])
seSig2/(yfit@results@respar[1]**2) # very high SE ?

##########################################################################################
# Alternate computation using the derivation of Wang and Merkle for merDeriv (JSS 2018) 
## ie the expression 0.5*tr(V^{-1} . dV/dsig_k . V^{-1} . dV/dsig_l))

# derivation of V w/r to omega2 and sigma
Vderiv <- list(Xmat %*% diag(c(1,0)) %*% t(Xmat), Xmat %*% diag(c(0,1)) %*% t(Xmat), diag(2*yfit@results@respar[1], nrow=length(xtim), ncol=length(xtim)))
# derivation of V w/r to omega2 and sigma2
Vderiv <- list(Xmat %*% diag(c(1,0)) %*% t(Xmat), Xmat %*% diag(c(0,1)) %*% t(Xmat), diag(1, nrow=length(xtim), ncol=length(xtim)))

blocBbis<-matrix(data=0, nrow=3, ncol=3)
for(i in 1:3) {
  for(j in 1:3) {
    blocBbis[i,j] <- 0.5*tr(invVi %*% Vderiv[[i]] %*% invVi %*% Vderiv[[j]])
  }
}
blocBbis<-blocBbis*nsuj

# Virtually identical to the previous computation
blocBbis-blocB

# Comparing to the FIM - very close for the variances but a factor 5 for sigma2... (or 10 for sigma when using the derivative w/r sigma instead of sigma2)
sqrt(diag(solve(yfit@results@fim)))[3:5]
varsig2<-2*yfit@results@respar[1]*diag(solve(yfit@results@fim))[5]
sqrt(varsig2)

sqrt(diag(solve(blocBbis)))

##########################################################################################
# Trying with different parameters

param<-c(15, 0.3)
omega <- param*.3 # 30% variability
sigma <- 3

# Now much closer !! probably sigma was too small... (but means saemix too optimistic by far about its ability to identify the SE)
# note: very poor estimate of omega(slope) for some reason (CV=582% with saemix, 230% with the exact version), maybe need a better design with later time points (can compute this with the expected FIM)
