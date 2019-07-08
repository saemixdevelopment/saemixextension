#' Sampling Importance Resampling (SIR) algorithm
#' 
#' SIR algorithm performs parameter uncertainty estimation for nonlinear mixed effects models
#' 
#' @name saemixSIR
#' @aliases SIR_methods
#' 
#' @param 
#' 
#' @export saemix

saemixSIR<-function(SaemixObject, M, m, inflation, est.mu, prop.distr, optionll, warnings) {
  if (missing(SaemixObject)) {
    cat('Missing SaemixObject argument')
    return("Creation of SaemixSIR object has failed")
  }
  if(class(SaemixObject)!="SaemixObject") {
    return("Please provide a valid SaemixObject object (see the help page for SaemixObject)\n")
  }
  SaemixSIR<-new(Class="SaemixSIR",SaemixObject=SaemixObject, M=M, m=m, inflation=inflation, optionll=optionll, est.mu=est.mu, prop.distr=prop.distr, warnings=warnings)
  
  
  ############################################
  #              Main Algorithm              #
  ############################################
  
  SaemixObject <- SaemixSIR['SaemixObject']
  M <- SaemixSIR['M']
  m <- SaemixSIR['m']
  est.mu <- SaemixSIR['est.mu']
  
  inflprop.distr <- SaemixSIR['inflprop.distr']
  optionll <- SaemixSIR['optionll']
  name.param <- SaemixSIR['name.param']
  
  warnings <- SaemixSIR['warnings']
  
  ####### SAMPLING #######
  list <- sampling(SaemixObject, M, inflprop.distr, optionll, warnings)
  sampled.theta <- list[[1]]
  OFVi <- list[[2]]
  samptries <- list[[3]]
  colnames(sampled.theta) <- name.param
  SaemixSIR@sampled.theta <- sampled.theta
  SaemixSIR@OFVi <- OFVi
  SaemixSIR@samptries <- samptries
  
  ####### IMPORTANCE RATIO ######
  IR <- importance.ratio(SaemixSIR)
  SaemixSIR@IR <- IR
  
  ####### RESAMPLING #######
  resamp <- resample(sampled.theta, IR, m, warnings)
  resampled.theta <- resamp[[1]]
  SaemixSIR@resampled.theta <- resampled.theta
  SaemixSIR@resamples.order <- resamp[[2]] 

  
  ####### HISTOGRAMS OF RESAMPLES #######
  for (i in 1:ncol(resampled.theta)){
    name <- colnames(resampled.theta)[i]
    title <- paste('Histogram of',name,'from SIR resamples')
    hist(resampled.theta[,i], main=title, xlab=name)
  }
  sd.new <- apply(resampled.theta, 2, sd) #standard deviation for each parameter 
  SaemixSIR@sdSIR<- sd.new
  return(SaemixSIR)
}
  

