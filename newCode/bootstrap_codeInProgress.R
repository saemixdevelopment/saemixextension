# Center conditional distribution and return etas to be sampled from
# population: centering method, either at the level of the population or at the level of the individual
centerDist.NPcond2<-function(saemixObject, nsamp, population=TRUE) {
  #  eps.est<-saemixObject@results@iwres # Estimated epsilon_ij
  phicond<-saemixObject@results@cond.mean.phi
  etacond<-saemixObject@results@cond.mean.eta # Estimated eta_i
  Cimu<-phicond-etacond  
  eta.samp<-data.frame() # Conditional samples from eta distribution
  for(i in 1:nsamp) {
    eta.samp<-rbind(eta.samp,saemixObject@results@phi.samp[,,i]-Cimu)
  }
  # TODO center also residuals individually
  epsc<-center.eps(eps.est)
  epscorr<-center.eps(eps.est)/sd(epsc) # Correcting empirical residuals
  
  # Centering residuals at the level of the population or at the level of the individual
  if(population) {
    eta.sampc<-center.eta(eta.samp)
  } else {
    id1<-rep(1:dim(phicond)[1],nsamp)
    eta.sampc<-eta.samp
    #    epscorr<-epsc<-eps.est
    for(isuj in 1:dim(phicond)[1]) {
      idx<-which(id1==isuj)
      eta.sampc[idx,]<-center.eta(eta.samp[idx,])
      # epsc[idx]<-center.eps(eps.est[idx])
      # epscorr[idx]<-center.eps(eps.est[idx])/sd(epsc[idx]) # Correcting empirical residuals
    } 
  }
  return(list(eta.sampc=eta.sampc,epsc=epscorr))
}

