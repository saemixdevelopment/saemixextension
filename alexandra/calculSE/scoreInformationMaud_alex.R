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
omega <- c(3,0.2) 
sigma <- 2

nsuj<-250
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

#write.table(xtab, file="C:/Users/AlexandraLAVALLEY/Documents/Code_saemix/essais/xtab.txt",row.names = F)

# note: very poor estimate of omega(slope) for some reason (CV=582% with saemix, 230% with the exact version), maybe need a better design with later time points (can compute this with the expected FIM)

##########################################################################################
# Computing FIM version Maud, at the end of the iterations, using the sufficient statistics

muhat <- yfit@results@fixed.effects # mu
omegahat <- sqrt(diag(yfit@results@omega)) # omega (=sqrt(omega2))
sighat <- yfit@results@respar[1] # sigma
pred = predict(yfit)  # moi ?a marche pas yfit@results@ypred


mu0 = muhat[1]
mu1 = muhat[2]
omega0 = omegahat[1]
omega1 = omegahat[2]
sigma = sighat

bigmat=matrix(data=0,nrow = 5,ncol = 5)

for (i in unique(yfit@data@data$id)){
  yobsi = yfit@data@data$y[yfit@data@data$id==i]
  ni = length(yobsi)
  
  dphi_mu0 = matrix(c(rep(0,ni),1/omega0**2,0,0,0))
  dphi_mu1 = matrix(c(rep(0,ni),0,1/omega1**2,0,0))
  dphi_omega0 = matrix(c(rep(0,ni),-mu0/(omega0**4),0,1/(2*(omega0**4)),0))
  dphi_omega1 = matrix(c(rep(0,ni),0,-mu1/(omega1**4),0,1/(2*(omega1**4))))
  dphi_sigma = matrix(c(rep(2/sigma**3,ni),0,0,0,0))
  
  Si = as.matrix(c(1/2 * ((yobsi-pred[yfit@data@data$id==i])**2), yfit@results@cond.mean.psi[i,], yfit@results@cond.mean.psi[i,]**2))
  
  delta_mu0 = sum(Si*dphi_mu0) - mu0/(omega0**2)
  delta_mu1 = sum(Si*dphi_mu1) - mu1/(omega1**2)
  delta_omega0 = sum(Si*dphi_omega0) -1/(2*omega0**2) + mu0**2/(2*omega0**4)
  delta_omega1 = sum(Si*dphi_omega1) -1/(2*omega1**2) + mu1**2/(2*omega1**4)
  delta_sigma = sum(Si*dphi_sigma) -1/(sigma) *ni
  
  mat = matrix(c(delta_mu0,delta_mu1,delta_omega0,delta_omega1,delta_sigma),ncol=1)
  bigmat = bigmat + mat%*%t(mat)
}

#bigmat = bigmat/length(unique(yfit@data@data$id))

sol = solve(bigmat)
sqrt(diag(sol)) 


