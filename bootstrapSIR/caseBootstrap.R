#########################################################################################
# Centering and normalising eta & eps
center.eps<-function(x) {
  x<-x-mean(x)
}
center.eta<-function(x) {
  for(i in 1:dim(x)[2]) x[,i]<-center.eps(x[,i])
  return(x)
}
normalise.eta<-function(eta.est, omega.est) {
  etacondc<-center.eta(eta.est)
  letacond<-t(chol(cov(etacondc)))
  letaest<-t(chol(omega.est))
  Acorr<-t(letaest %*% solve(letacond))
  etacorr<-etacondc %*% Acorr
  return(etacorr)
}

###### TODO - add checks for singular matrices
normalise.eta.svd<-function(eta.est, omega.est) {
  etacondc<-center.eta(eta.est)
  Semp<-1/(dim(etacondc)[1]-1) * t(etacondc) %*% etacondc
  Sempd<-svd(Semp)
  Sestd <- svd(omega.est)
  demp<-1/sqrt(Sempd$d)
  dest<-sqrt(Sestd$d)
  Acorr <- Sempd$v %*% diag(demp,nrow=length(demp)) %*%  diag(dest,nrow=length(dest)) %*% t(Sestd$v)
  etacorr<-etacondc %*% Acorr
  return(etacorr)
}

#########################################################################################
# Functions from saemix package
transphi<-function(phi,tr) {
  psi<-as.matrix(phi)
  #  if(is.null(dim(psi))) psi<-as.matrix(t(psi),nrow=1)
  i1<-which(tr==1) # log-normal
  psi[,i1]<-exp(psi[,i1])
  i2<-which(tr==2) # probit
  psi[,i2]<-normcdf(psi[,i2])
  i3<-which(tr==3) # logit
  psi[,i3]<-1/(1+exp(-psi[,i3]))
  if(is.null(dim(phi))) psi<-c(psi)
  return(psi)
}

cutoff<-function(x,seuil=.Machine$double.xmin) {x[x<seuil]<-seuil; return(x)}

normcdf<-function(x,mu=0,sigma=1)
  cutoff(pnorm(-(x-mu)/sigma,lower.tail=FALSE),1e-30)

error.typ<-function(f,ab) {
  g<-cutoff(ab[1]+ab[2]*abs(f))
  return(g)
}

#########################################################################################
# Bootstrap saemix

# Generating bootstrap datasets
dataGen.case<- function(saemix.fit) {
  idx.boot<-sample(1:saemix.fit@data@N, size=saemix.fit@data@N, replace=T)
  smx.data<-saemix.fit@data
  data.boot<-NULL
  for(id in idx.boot)
    data.boot<-rbind(data.boot, saemix.fit@data@data[saemix.fit@data@data$index==id,])
  smx.data@ntot.obs<-dim(data.boot)[1]
  smx.data@nind.obs <-saemix.fit@data@nind.obs[idx.boot]
  smx.data@data<-data.boot
  orig.idx<-1:saemix.fit@data@N
  smx.data@data[,"index"]<-smx.data@data[,saemix.fit@data@name.group] <- rep(1:smx.data@N,times=smx.data@nind.obs)
  # As yorig and ocov not used in fit, not included in dataset
  return(smx.data)
}
