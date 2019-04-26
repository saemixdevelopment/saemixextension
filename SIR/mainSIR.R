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

saemixSIR<-function(SaemixObject, M, m, inflation, est.mu, cov.mat, optionll, warnings) {
  if (missing(SaemixObject)) {
    cat('Missing SaemixObject argument')
    return("Creation of SaemixSIR object has failed")
  }
  if(class(SaemixObject)!="SaemixObject") {
    return("Please provide a valid SaemixObject object (see the help page for SaemixObject)\n")
  }
  SaemixSIR<-new(Class="SaemixSIR",SaemixObject=SaemixObject, M=M, m=m, inflation=inflation, optionll=optionll, est.mu=est.mu, cov.mat=cov.mat, warnings=warnings)
  
  
  ############################################
  #              Main Algorithm              #
  ############################################
  
  SaemixObject <- SaemixSIR['SaemixObject']
  M <- SaemixSIR['M']
  m <- SaemixSIR['m']
  est.mu <- SaemixSIR['est.mu']
  
  inflcov.mat <- SaemixSIR['inflcov.mat']
  optionll <- SaemixSIR['optionll']
  name.param <- SaemixSIR['name.param']
  
  warnings <- SaemixSIR['warnings']
  
  ####### SAMPLING #######
  list <- sampling(SaemixObject, M, inflcov.mat, optionll, warnings)
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
  theta <- resample(sampled.theta, IR, m, warnings)
  SaemixSIR@resampled.theta <- theta[['resampled_theta']]
  SaemixSIR@order <- as.vector(theta[['order_resampled']])
  
  ####### HISTOGRAMS OF RESAMPLES #######
  for (i in 1:ncol(theta[['resampled_theta']])){
    name <- colnames(theta[['resampled_theta']])[i]
    title <- paste('Histogram of',name,'from SIR resamples')
    hist(theta[['resampled_theta']][,i], main=title, xlab=name)
  }
  sd.new <- apply(theta[['resampled_theta']], 2, sd) #standard deviation for each parameter 
  SaemixSIR@sdSIR<- sd.new
  return(SaemixSIR)
}
  

