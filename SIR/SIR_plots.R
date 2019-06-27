###########################################################################################
#############################    SIR  Diagnostic Graphs    ################################
###########################################################################################


###########################################################################################
###############################   dOFV distribution plot    ###############################

# dOFV distribution plot corresponds to the plot of dOFV distributions obtained from the
# M proposal samples and from the m SIR resamples against a Chi square distribution with 
# degrees of freedom equal to the number of estimated parameters. It informs about the 
# adequacy of the proposal distribution and M/m. 

# DIAGNOSTIC : sampling
# - If the dOFV distribution obtained from the M samples is at or close to the Chi square 
# distribution, it means that the proposal distribution is close to the true distribution.
# - If it is far above or below, it means that it is quite different from the true 
# distribution

# DIAGNOSTIC : resampling
# - If the dOFV ditribution obtained from the m SIR resamples is at or below the Chi square
# distribution, M/m may have been sufficient and should be further investigated with the 
# temporal plots. 
# - If the dOFV distribution is above the Chi square distribution, M/m was not sufficient


dOFVsamp <- function(SaemixSIR){
  SaemixObject <- SaemixSIR['SaemixObject']
  optionll <- SaemixSIR['optionll']
  dOFVvec <- c(0*length(nrow(sample)))
  if (optionll=='linearisation'){
    ll <- SaemixObject["results"]["ll.lin"]
  }
  if (optionll=='importance_sampling'){
    ll <- SaemixObject["results"]["ll.is"]
  }
  if (optionll=='gaussian_quadrature'){
    ll <- SaemixObject["results"]["ll.gq"]
  }
  OFVi <- SaemixSIR['OFVi']
  #compute the difference between OFVi and OFV of the vector with maximum likelihood 
  dOFVvec <- OFVi-(-2*ll)
  dOFVvec[dOFVvec<0] <- dOFVvec[dOFVvec<0]*-1
  dOFVvec <- dOFVvec[dOFVvec!=-2*ll]
  return(dOFVvec)
}


dOFVresamp <- function(SaemixSIR, resample, warn){
  SaemixObject <- SaemixSIR['SaemixObject']
  optionll <- SaemixSIR['optionll']
  dOFVvec <- c(0*length(nrow(sample)))
  if (optionll=='linearisation'){
    ll <- SaemixObject["results"]["ll.lin"]
  }
  if (optionll=='importance_sampling'){
    ll <- SaemixObject["results"]["ll.is"]
  }
  if (optionll=='gaussian_quadrature'){
    ll <- SaemixObject["results"]["ll.gq"]
  }
  
  nresamples <- c()
  for (i in 1:nrow(resample)){
    for (j in 1:nrow(SaemixSIR['sampled.theta'])){
      samp <- SaemixSIR['sampled.theta'][j,]
      resamp <- resample[i,]
      if (sum(samp==resamp)==length(SaemixSIR['name.param'])){
        nresamples <- c(nresamples, j)
      }
    }  
  }
  OFVi <- SaemixSIR['OFVi'][nresamples]
  #compute the difference between OFVi and OFV of the vector with maximum likelihood 
  dOFVvec <- OFVi-(-2*ll)
  dOFVvec[dOFVvec<0] <- dOFVvec[dOFVvec<0]*-1
  dOFVvec <- dOFVvec[dOFVvec!=-2*ll]
  return(dOFVvec)
}




#Plot function of dOFV distribution plot
dOFV_distribution_plot <- function(SaemixSIR, warn=TRUE){
  #degrees of freedom = number of parameters estimated
  SaemixObject <- SaemixSIR['SaemixObject']
  df <- length(SaemixSIR["name.param"])
  sample <- SaemixSIR["sampled.theta"]
  
  dOFV.proposal <- dOFVsamp(SaemixSIR) #dOFV of sampled.theta
  dOFV.proposal <- sort(dOFV.proposal)

  #create a vector of x from 0 to 1 with equally distributed ticks = distribution quantiles
  xprop <- seq(0,1, 1/(length(dOFV.proposal)-1))
  #estimation of degrees of freedom = best estimator = mean 
  dfprop <- round(mean(dOFV.proposal), digits=1)
  
  resample <- SaemixSIR["resampled.theta"]
  dOFV.SIR <- dOFVresamp(SaemixSIR, resample, warn)
  dOFV.SIR <- sort(dOFV.SIR) #dOFV of resampled.theta

  xsir <- seq(0,1, 1/(length(dOFV.SIR)-1))
  dfsir<- round(mean(dOFV.SIR), digits=1)
  
  ref <- paste(df, '(REF)')
  prop <- paste(dfprop, '(PROPOSAL)')
  sir <- paste(dfsir, '(SIR)')
  
  #plot of chi-square distribution
  par(mfrow=c(1,1))
  p<- seq(0,0.9999,0.0005)
  plot(qchisq(p,df)~p, type='l', xlab='Distribution quantiles', xaxp=c(0,1,4), ylab="dOFV", col="black",
       main='dOFV distribution plot')
  abline(h=seq(0,50,10), v=seq(0,1, 0.25), col="gray", lty=3)
  lines(dOFV.proposal~xprop, col='blue', lty=2)
  lines(dOFV.SIR~xsir, col='red')
  
  legend(0, 25, legend=c('Estimated df', ref,prop, sir),
         col=c("black", "black", "blue", "red"), lty=c(0,1,2,1), cex=0.8,
         box.lty=0)

  
}
###########################################################################################
###########################################################################################

###########################################################################################
###############################     Spatial trends plot     ###############################

# The spatial trends plot indicates how the proposal differs from the true distribution,
# but it does not inform whether SIR was able to compensate for these differences

# The plots show the resampling proportion (the number of resampled parameters divided by 
# the number of parameters available from the M samples), in different regions (bins) of 
# the parameter space (x-axis) (all values between the lowest and highest sampled parameter
# value). The parameter space is divided in 10 bins which all contain the same number of 
# samples. 

# - Horizontal trend (no trend) : the proposal distribution is close to the true uncertainty.
# - Bell-shaped trend : observed proportion higher in the center and lower at the ends, 
# parameters close to the final estimates are resampled more often than those further away
# from them : the proposal distribution is wider than the true distribution.
# - U-shaped trend : observed proportion lower in the center and  higher at the ends, the 
# proposal distribution is narrower than the true distribution.
# - Diagonal trend : observed proportion higher at one end, lower at the other : the proposal
# distribution has a different (a)symmetry than the true distribution.

spatial_trend_proportion <- function(samp.param, resamp.param, M, m){
  samp.param <- sort(samp.param)

  proportion <- matrix(ncol=2, nrow=10)
  #second column = number of spatial bin
  proportion[,1] <- 1:10
  
  for (i in 1:10){
    #separation of sampled and resampled parameters in 10 bins
    #number of resampled parameters in bin x that are in bin x of the sampled parameters
    s <- sum(resamp.param %in% samp.param[(((i-1)*M/10)+1):(i*M/10)])
    #proportion of resamples from the bin x of samples (which has M/10 samples)
    proportion[i,2] <- s/(M/10)
  }
  return(proportion)
}


spatial_trend_plot <- function(SaemixSIR){
  samp <- SaemixSIR['sampled.theta']
  resamp <- SaemixSIR['resampled.theta']
  M <- SaemixSIR['M']
  m <- SaemixSIR['m']
  name.param <- SaemixSIR['name.param']
  npar.est <- length(SaemixSIR['est.mu'])
  par(mfrow=c(3,3)) 
  for (i in 1:npar.est){
    samp.param <- samp[,i]
    resamp.param <- resamp[,i]
    proportion <- spatial_trend_proportion(samp.param, resamp.param, M, m)
    title <- name.param[i]
    
    #stochastic noise around expected proportion Binomial(m/M (expected proportion))
    expprop <- m/M
    se <- sqrt((expprop*(1-expprop))/(M/10))
    
    plot(x=NA, y=NA, xlim=c(1,10),
         xlab='Spatial bins of initial samples', ylab='Proportion resampled', 
         ylim=c(0,0.25), panel.first = grid(), las=1)
    title(title, line = 0.5)
    polygon(x=c(1,10,10,1), y=c(expprop+1.96*se, expprop+1.96*se,expprop-1.96*se,expprop-1.96*se), col='lightgray', border=NA)
    lines(proportion[,2]~proportion[,1], pch=20, type='o')
    abline(h=expprop, lty=2)
  }
  mtext("Spatial trends plot: adequacy of proposal", side = 3, line = -1.5, outer = TRUE)
}
###########################################################################################
###########################################################################################


###########################################################################################
###############################    Temporal trends plot     ###############################

# The temporal trends plot indicates, for each parameter, whether M/m was high enough to
# compensate for the differences between the proposal distribution and the true uncertainty

# The plot corresponds for each parameter to resamples from the bin in the spatial trends
# plot which had the highest proportion resampled. The spatial bin used is the region where 
# SIR is most likely to "run out" of  good samples, if M/m is not sufficient since resampling 
# is without replacement. The x-axis corresponds to 5 time bins. The plot shows for each bin 
# the number of resamples that belongs to the top spatial bin.
# - Horizontal trend (no trend) : M/m was sufficient to compensate for potential differences 
# between the proposal distribution and the true uncertainty
# - Downward diagonal trend (decreases over time) : it indicates a depletion of samples in 
# the top bin, there were not enough good samples in the SIR procedure to fully correct the 
# proposal uncertainty: M/m was not sufficient.


temporal_trend_tab <- function(SaemixSIR, max.resamp, i){
  resamp.param <- SaemixSIR['resampled.theta'][,i]
  m <- SaemixSIR['m']
  trend <- matrix(nrow=5, ncol=2)
  trend[,1] <- 1:5
  ordermaxresamp <- rep(0,length(max.resamp))
  for (i in 1:length(max.resamp)){
    ind <- which(resamp.param==max.resamp[i])
    ordermaxresamp[i] <- ind
  }
  for (i in 1:5){
    s <- sum(ordermaxresamp %in% ((((i-1)*(m/5))+1):(i*(m/5))))
    trend[i,2] <- s
  }
  return(trend)
}


temporal_trend_plot <- function(SaemixSIR){
  M <- SaemixSIR['M']
  m <- SaemixSIR['m']
  samp <- SaemixSIR['sampled.theta']
  resamp <- SaemixSIR['resampled.theta']
  name.param <- SaemixSIR['name.param']
  par(mfrow=c(3,3)) 
  for (i in 1:ncol(resamp)){
    samp.param <- samp[,i]
    samp.param <- sort(samp.param)
    resamp.param <- resamp[,i]
    proportion <- spatial_trend_proportion(samp.param, resamp.param, M, m)
    ind.max <- which.max(proportion[,2])
    max <- proportion[ind.max,2]
    max.resamp <- resamp.param[resamp.param %in% samp.param[(((ind.max-1)*M/10)+1):(ind.max*M/10)]]
    trend <- temporal_trend_tab(SaemixSIR, max.resamp, i)
    
    h <- (max*M/10)*(1/5) # number of resamples from the max spatial bin divided by the number of temporal bins
    se <- sqrt((max*(1-max))/h)
    se <- se*h
    title <- name.param[i]
    subtitle <- paste('Spatial bin', ind.max)
    plot(x=NA, y=NA, sub=subtitle,
         xlab='Temporal bin of resamples', ylab='Parameters resampled count', 
         xlim=c(1,5), ylim=c(0,m/10), panel.first = grid(), las=1)
    
    polygon(x=c(1,5,5,1), y=c(h+1.96*se, h+1.96*se,h-1.96*se,h-1.96*se), col='lightgray', border=NA)
    lines(trend[,2]~trend[,1], pch=19, type='o')
    abline(h=h, lty=2)
    title(title, line = 0.5)
  }
  mtext("Temporal trends plot: exhaustion of samples", side = 3, line = -1.5, outer = TRUE)
}
###########################################################################################
###########################################################################################

