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
omega <- param*.3 # 30% variability
sigma <- 3

nsuj<-50
xtim <- c(0, 2, 5, 10, 24, 50)
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

sqrt(diag(invFIM))[5]/yfit@results@respar[1]

# note: very poor estimate of omega(slope) for some reason (CV=582% with saemix, 230% with the exact version), maybe need a better design with later time points (can compute this with the expected FIM)

##########################################################################################
# Computing FIM version Maud, at the end of the iterations, using the sufficient statistics

muhat <- yfit@results@fixed.effects # mu
omegahat <- sqrt(diag(yfit@results@omega)) # omega (=sqrt(omega2))
sighat <- yfit@results@respar[1] # sigma
ntot <- yfit@data@ntot.obs

# s(y,phi)
s1<-data.frame(id=yfit@data@data$id, s1=((yfit@data@data$y-yfit@results@ypred)**2)/2)
s3<-s2<-data.frame(id=rep(1:nsuj, 2), s2=c(yfit@results@cond.mean.psi))
s3[,2]<-s2[,2]**2


# dPsi/dtheta
dpsi <- c( muhat/(omegahat**2), nsuj/(2*omegahat)-(muhat**2/(2*(omegahat**4))), ntot/sighat)

# dPsi/dtheta versin individuelle
dpsi <- c( muhat/(omegahat**2), 1/(2*omegahat)-(muhat**2/(2*(omegahat**4))), length(xtim)/sighat)

# dPhi/dtheta
# same dphi for each subject... (same design)
dphi1.sig <- rep(2/sighat**3, length(xtim))
dphi2.mu <- matrix(c(1/(omegahat[1]**2),0,0,1/(omegahat[2]**2)), ncol=2)
dphi2.omega2 <- matrix(c(-muhat[1]/(omegahat[1]**4), 0,0,-muhat[2]/(omegahat[2]**4)), ncol=2)
dphi3.omega2 <- matrix(c(1/(2*omegahat[1]**4),0,0,1/(2*omegahat[2]**4)), ncol=2)

isuj<-1
delta.mu <- dphi2.mu %*% s2[s2$id==isuj,2]
delta.omega <- dphi2.omega2 %*% s2[s2$id==isuj,2] + dphi3.omega2 %*% s3[s3$id==isuj,2]
delta.sig <- sum(dphi1.sig * s1[s1$id==isuj,2])
dphi <- c(delta.mu, delta.omega, delta.sig)
dphi %*% t(dphi)
dpsi %*% t(dpsi)
try(solve(dphi %*% t(dphi))) # can't really inverse the score for one subject
try(solve(dpsi %*% t(dpsi)))

delta <- - dpsi + dphi

score <- delta %*% t(delta)
try(solve(score))
1/diag(score)

# Over the N subjects
xcal <- matrix(data=0, ncol=5, nrow=5)
for(isuj in 1:nsuj) {
  delta.mu <- dphi2.mu %*% s2[s2$id==isuj,2]
  delta.omega <- dphi2.omega2 %*% s2[s2$id==isuj,2] + dphi3.omega2 %*% s3[s3$id==isuj,2]
  delta.sig <- sum(dphi1.sig * s1[s1$id==isuj,2])
  dphi <- c(delta.mu, delta.omega, delta.sig)
  delta <- - dpsi + dphi
  
  score <- delta %*% t(delta)
  xcal <- xcal + score
}
xcal<-xcal/nsuj
1/diag(xcal)
maudVar <- solve(xcal)
sqrt(diag(maudVar))
sqrt(diag(invFIM))

# RSE
sqrt(diag(maudVar))[1:2]/(muhat)
sqrt(diag(maudVar))[3:4]/(omegahat**2)
sqrt(diag(maudVar))[5]/sighat

# Some checks
muhat/omegahat
muhat**2/omegahat**3
muhat**2/omegahat**4
1/sighat/nsuj

# Quick check
# Delta for mu should be just (phi_i - mu)/omega2 
summary((s2[1:50,2]-muhat[1])/omegahat[1]**2)
sum(((s2[1:50,2]-muhat[1])/omegahat[1]**2)**2)
sum(((s2[51:100,2]-muhat[2])/omegahat[2]**2)**2)
# Delta for omega should be (phi_i - mu)**2/(omega**4) -1/(2*omega) 
summary((s2[1:50,2]-muhat[1])**2/omegahat[1]**4)
summary((s2[51:100,2]-muhat[2])**2/omegahat[2]**4)
sum(((s2[1:50,2]-muhat[1])/omegahat[1]**2)**2)
sum(((s2[51:100,2]-muhat[2])/omegahat[2]**2)**2)

##########################################################################################
# Deriving w/r omega instead of omega2 - similar results

muhat <- yfit@results@fixed.effects # mu
omegahat <- sqrt(diag(yfit@results@omega)) # omega (=sqrt(omega2))
sighat <- yfit@results@respar[1] # sigma
ntot <- yfit@data@ntot.obs

# s(y,phi)
s1<-data.frame(id=yfit@data@data$id, s1=((yfit@data@data$y-yfit@results@ypred)**2)/2)
s3<-s2<-data.frame(id=rep(1:nsuj, 2), s2=c(yfit@results@cond.mean.psi))
s3[,2]<-s2[,2]**2


# dPsi/dtheta versin individuelle
dpsi <- c( muhat/(omegahat**2), 1/omegahat-(muhat**2)/(omegahat**3), length(xtim)/sighat)

# dPhi/dtheta
# same dphi for each subject... (same design)
dphi1.sig <- rep(2/sighat**3, length(xtim))
dphi2.mu <- matrix(c(1/(omegahat[1]**2),0,0,1/(omegahat[2]**2)), ncol=2)
dphi2.omega <- matrix(c(-2*muhat[1]/(omegahat[1]**3), 0,0,-2*muhat[2]/(omegahat[2]**3)), ncol=2)
dphi3.omega <- matrix(c(1/(omegahat[1]**3),0,0,1/(omegahat[2]**3)), ncol=2)

xcal <- matrix(data=0, ncol=5, nrow=5)
for(isuj in 1:nsuj) {
  delta.mu <- dphi2.mu %*% s2[s2$id==isuj,2]
  delta.omega <- dphi2.omega %*% s2[s2$id==isuj,2] + dphi3.omega %*% s3[s3$id==isuj,2]
  delta.sig <- sum(dphi1.sig * s1[s1$id==isuj,2])
  dphi <- c(delta.mu, delta.omega, delta.sig)
  delta <- - dpsi + dphi
  
  score <- delta %*% t(delta)
  xcal <- xcal + score
}
xcal<-xcal/nsuj
1/diag(xcal)
maudVar <- solve(xcal)
sqrt(diag(maudVar))
sqrt(diag(invFIM))

sqrt(diag(maudVar))[1:2]/(muhat)
sqrt(diag(maudVar))[3:4]/(omegahat)
sqrt(diag(maudVar))[5]/sighat

