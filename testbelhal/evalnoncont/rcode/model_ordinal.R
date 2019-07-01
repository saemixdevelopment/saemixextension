# Model parameters Intercept and beta
parpop<-c(0,-0.37)
nampar<-c("Intercept","beta.time")
omega<-diag(c(.21,.1))

pvrai<-c(0,-0.37,.21,.1)
pfaux<-c(10,-8.37,.3,.3)

# Cat data models
binary.model<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  logpdf<-rep(0,length(tim))
  P.obs = (y==0)*(1-pevent)+(y==1)*pevent
  logpdf <- log(P.obs)logpdf <- log(P.obs)
  return(logpdf)
}
