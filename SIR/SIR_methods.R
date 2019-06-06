#####################################################################
####################### FUNCTIONS FOR SaemixSIR #####################
#####################################################################


###### SAMPLING ######

sampling <- function(SaemixObject, M=5000, inflcov.mat, optionll, warnings){
  if(warnings) cat('Sampling and computation of OFVs...\n')
  if ("mvtnorm" %in% rownames(installed.packages())==F) {
    install.packages("mvtnorm")
  }
  library('mvtnorm', character.only = T)
  indx.fix <- SaemixObject['results']['indx.fix']
  indx.cov <- SaemixObject['results']['indx.cov']
  npar.est <- length(estpar.vector(SaemixObject))
  indx.fixed <- c(indx.fix, indx.cov)
  nfix <- length(indx.fixed)

  se.parfix<-inflcov.mat[1:nfix,1:nfix]
  se.parvar<-inflcov.mat[(nfix+1):npar.est,(nfix+1):npar.est]
  parfix <- t(SaemixObject['results']['fixed.effects'])
  par <- estpar.vector(SaemixObject)
  parvar <- par[-indx.fixed]
  OFVi <- c()
  samples <- matrix(ncol=npar.est, nrow=0)
  i<-0 #keeping track of number of tries of sampling
  while(nrow(samples)!=M){
    i<- i+1
    if (i==M*3) {
      cat('Sampling is taking too long, inflation could be too high, or est.mu is not close to reality')
      break}
    samplefix <- rmvnorm(1,mean=parfix,sigma=se.parfix)
    samplevar <- rmvnorm(1, mean=parvar, sigma=se.parvar)
    sample <- cbind(samplefix, samplevar)
    ofvi <- OFVi(SaemixObject, sample, optionll)
    if (ofvi==0) next #ofvi=0 : error while computing OFVi
    samples <- rbind(samples, sample)
    OFVi <- c(OFVi, ofvi)
    l <- nrow(samples)
    if(warnings){
      if(l==M%/%5) cat(paste(l, '/', M, 'done...\n'))
      if(l==2*M%/%5) cat(paste(l, '/', M, 'done...\n'))
      if(l==3*M%/%5) cat(paste(l, '/', M, 'done...\n'))
      if(l==4*M%/%5) cat(paste(l, '/', M, 'done...\n'))
      if(l==M) cat(paste(l, 'samples and OFVs done.\n'))
    }
    

  }
  return(list(sampled.theta=samples, OFVi=OFVi, samptries=i))
}



###### IMPORTANCE RATIO ######

OFVi <- function(SaemixObject, sampled.theta, optionll, warn=FALSE){
  #vector for the final ofvi of each parameter vector
  ofvi <- matrix(nrow = nrow(sampled.theta), ncol=1)
  ll <- rep(0, nrow(sampled.theta))
  
  if(warn)cat('Starting OFV computation for resamples...\n')
  length <- nrow(sampled.theta)
  
  if (optionll=='importance_sampling'){
    for (i in 1:length){
      sobj <- SaemixObject
      sobj <- replacePopPar.saemixObject(sobj, sampled.theta[i,])
      ll[i] <- try(llis.saemix(sobj)['results']['ll.is'],silent=T) #the result of llis.saemix is a saemixObject
      if(!is.numeric(ll[i])){ll[i]<- 0}
      ll <- as.numeric(ll)
      if(warn){
        if(i==length%/%5) cat(paste(i, '/', length, 'OFVs done...\n'))
        if(i==2*length%/%5) cat(paste(i, '/', length, 'OFVs done...\n'))
        if(i==3*length%/%5) cat(paste(i, '/', length, 'OFVs done...\n'))
        if(i==4*length%/%5) cat(paste(i, '/', length, 'OFVs done...\n'))
        if(i==length) cat(paste(i, 'OFVs done.\n'))
      }
    }
  }
  
  if (optionll=='linearisation'){
    for (i in 1:nrow(sampled.theta)){
      sobj <- SaemixObject
      sobj <- replacePopPar.saemixObject(sobj, sampled.theta[i,])
      ll[i] <- try(fim.saemix(sobj)['results']['ll.lin'], silent=T)
      if(!is.numeric(ll[i])){ll[i]<- 0}
      ll <- as.numeric(ll)
      if(warn){
        if(i==length%/%5) cat(paste(i, '/', length, 'OFVs done...\n'))
        if(i==2*length%/%5) cat(paste(i, '/', length, 'OFVs done...\n'))
        if(i==3*length%/%5) cat(paste(i, '/', length, 'OFVs done...\n'))
        if(i==4*length%/%5) cat(paste(i, '/', length, 'OFVs done...\n'))
        if(i==length) cat(paste(i, 'OFVs done.\n'))
      }
    }
  }
  ll <- as.numeric(ll)
  ofvi <- -2*ll #vector of OFVi of the sampled.theta
  return(ofvi)
}


importance.ratio <- function(SaemixSIR){
  warnings <- SaemixSIR['warnings']
  if(warnings)cat('IR computing...')
  SaemixObject <- SaemixSIR['SaemixObject']
  sampled.theta <- SaemixSIR['sampled.theta']
  optionll <- SaemixSIR['optionll']
  est.mu <- SaemixSIR['est.mu']
  inflcov.mat <- SaemixSIR['inflcov.mat']
  warnings <- SaemixSIR['warnings']
  
  if (optionll=='importance_sampling'){
    ll <- SaemixObject["results"]["ll.is"]
  }
  if (optionll=='linearisation'){
    ll <- SaemixObject["results"]["ll.lin"]
  } 

  OFVml <- -2*ll
  OFVi <- SaemixSIR['OFVi']
  
  num <- exp(-(1/2)*(OFVi-OFVml))
  num[OFVi==0] <- 0
  den <- c()
  invinflcov.mat <- solve(inflcov.mat)
  for (i in 1:nrow(sampled.theta)){
    den <- c(den, exp(-(1/2)*(sampled.theta[i,]-est.mu)%*%invinflcov.mat%*%as.matrix(sampled.theta[i,]-est.mu)))
  }
  IR <- num/den
  if(warnings)cat('done.\n')
  return(IR)
}


###### RESAMPLING (without replacement) ######
resample <- function(sampled.theta, IR, m, warnings){
  if(warnings)cat('Resampling...')
  resamptheta <- matrix(nrow = m, ncol=ncol(sampled.theta))
  for (i in 1:m){
    p.resample <- IR/sum(IR) #without replacement
    p.resample[is.na(p.resample)] <- 0
    indice <- sample.int(nrow(sampled.theta), size=1, replace=TRUE, prob=p.resample)
    resamptheta[i,1:ncol(sampled.theta)] <- sampled.theta[indice,]
    IR[indice] <- 0
  }
  colnames(resamptheta) <- colnames(sampled.theta)
  if(warnings)cat('done.\n')
  return(resamptheta)
}




