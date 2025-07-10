
llgq.saemix.ind<-function(saemixObject) {
# llgq.saemix # Compute LL for the whole dataset and update saemixObject
# llgq.saemix.ind # Compute LL up to a constant for each subject and return a vector
  
  # RES = MLXGQ(RES) Estimate the log-likelihood using Gaussian Quadrature (multidimensional grid)
  nnodes.gq<-saemixObject["options"]$nnodes.gq  # number of nodes on each 1-D grid
  nsd.gq<-saemixObject["options"]$nsd.gq  # the integral is computed on the interval [E(eta|y) +- nsd_gq*SD(eta|y)]
  
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  xind<-saemix.data["data"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
  yobs<-saemix.data["data"][,saemix.data["name.response"]]
  
  i1.omega2<-saemixObject["model"]["indx.omega"]
  Omega<-saemix.res["omega"]
  pres<-saemix.res["respar"]
  cond.var.phi<-saemix.res["cond.var.phi"]
  cond.mean.phi<-saemix.res["cond.mean.phi"]
  nphi1<-length(i1.omega2)
  IOmega.phi1<-solve(Omega[i1.omega2,i1.omega2,drop=FALSE])
  mean.phi1<-saemix.res["mean.phi"][,i1.omega2,drop=FALSE]
  
  io<-matrix(0,nrow=saemix.data["N"],ncol=max(saemix.data["nind.obs"]))
  for(isuj in 1:saemix.data["N"])
    io[isuj,1:saemix.data["nind.obs"][isuj]]<-1
  ind.io <- which(t(io)!=0)
  DYF<-matrix(0,nrow=dim(io)[2],ncol=dim(io)[1])
  
  phi<-saemix.res["mean.phi"]
  y<-gqg.mlx(nphi1,nnodes.gq)
  x<-(y$nodes-0.5)*2
  w<-(y$weights)*(2**nphi1)
  # ECO TODO check dimensions (unclear in matlab)
  nx<-dim(x)[1]
  condsd.eta<-sqrt(cond.var.phi[,i1.omega2,drop=FALSE])
  xmin<-cond.mean.phi[,i1.omega2,drop=FALSE]-nsd.gq*condsd.eta
  xmax<-cond.mean.phi[,i1.omega2,drop=FALSE]+nsd.gq*condsd.eta
  a<-(xmin+xmax)/2
  b<-(xmax-xmin)/2
  log.const<-0
  idx.exp<-which(saemixObject["model"]["error.model"]=="exponential")
  if(length(idx.exp)>0)
    #	if(saemixObject["model"]["error.model"]=="exponential")
    log.const<-(-sum(yobs[saemix.data["data"][,"ytype"] %in% idx.exp]))
  
  Q<-0
  for (j in 1:nx) {
    phi[,i1.omega2] <- a+b*matrix(rep(x[j,],saemix.data["N"]),ncol=nphi1,byrow=TRUE)
    psi<-transphi(phi,saemixObject["model"]["transform.par"])
    if(saemixObject["model"]["modeltype"]=="structural"){
      f<-saemixObject["model"]["model"](psi, saemix.data["data"][,"index"], xind)
      for(i in idx.exp) f[saemix.data["data"][,"ytype"]==i]<-log(cutoff(f[saemix.data["data"][,"ytype"]==i]))
      g<-error(f,pres,saemix.data["data"][,"ytype"])
      DYF[ind.io] <- -0.5*((yobs-f)/g)**2 - log(g)
      ly<-colSums(DYF)
    } else {
      f<-saemixObject["model"]["model"](psi, saemix.data["data"][,"index"], xind)
      DYF[ind.io] <- f
      ly<-colSums(DYF)
    }
    dphi1<-phi[,i1.omega2,drop=FALSE]-saemix.res["mean.phi"][,i1.omega2,drop=FALSE]
    lphi1<-(-0.5)*rowSums((dphi1%*%IOmega.phi1)*dphi1)
    ltot<-ly+lphi1
    ltot[is.na(ltot)]<-(-Inf)
    Q<-Q+w[j]*exp(ltot)
  }
  ll.ind <- log(Q)+rowSums((log(b))) - log(det(Omega[i1.omega2,i1.omega2, drop=FALSE]))/2
  # -(nphi1*log(2*pi)+ saemix.data["nind.obs"]*log(2*pi))/2 + log.const* saemix.data["nind.obs"]/(length(yobs[saemix.data["data"][,"ytype"] %in% idx.exp])) # check formula...
  ll.ind2 <- ll.ind -(nphi1*log(2*pi)+ saemix.data["nind.obs"]*log(2*pi))/2 
  if(length(idx.exp)>0) ll.ind2 <- ll.ind2 + log.const*saemix.data["nind.obs"]/(length(yobs[saemix.data["data"][,"ytype"] %in% idx.exp]))
  
#    S<-saemix.data["N"]*log(det(Omega[i1.omega2,i1.omega2, drop=FALSE]))+ saemix.data["N"]*nphi1*log(2*pi)+ saemix.data["ntot.obs"]*log(2*pi)
#    ll<-(-S/2) + sum(log(Q)+rowSums(log(b)))+ log.const
    return(ll.ind2)
  
}
