###########################################
# Folders
# When in my directory - set to FALSE to run from getwd(), must be paperSaemix3 (all paths relative)
EcoAtHome <- FALSE

# Libraries
library(saemix)
library(ggplot2)
library(gridExtra)
# library(MASS)
# library(rlang)

if(EcoAtHome) {
  # Load updated files
  # pending compilation and CRAN upload
  workDir<-"/home/eco/work/saemix/saemixextension/paperSaemix3"
} else {
  workDir<-getwd() 
}

# Setting up relative folders
setwd(workDir)
bootResDir <- file.path(workDir, "bootstrapRes")
figDir <- file.path(workDir,"figs")

###########################################
# Whether to save the plots
saveFigs<-FALSE
rsize.text <- 2
rsize.ticks <- 1.5
theme_set(theme_bw(base_size = 15)) 

###########################################
# Number of bootstrap samples
runBootstrap <- FALSE # to read the results from disk
nboot <-10
# nboot <- 200

# Avoiding the tidyverse...
regroupTable <- function(x, value=1, col=c(2,3), fun=list(mean=mean())) {
  tab<-NULL
  c1<-col[1]
  c2<-col[2]
  for(i1 in sort(unique(x[,c1]))) {
    for(i2 in sort(unique(x[,c2]))) {
      yvec <- x[x[,c1]==i1 & x[,c2]==i2,value]
      l1 <- c(i1, i2)
      for(ifun in 1:length(fun)) l1<-c(l1,fun[[ifun]](yvec))
      tab<-rbind(tab, l1)
    }
  }
  if(is.null(names(fun))) names(fun)<-1:length(fun)
  colnames(tab)<-c(col,names(fun))
  return(as.data.frame(tab))
}

################################################### Data - RAPI
data(rapi.saemix)
rapi.saemix$gender <- ifelse(rapi.saemix$gender=="Men",1,0) # Female=reference class as in Atkins

saemix.data<-saemixData(name.data=rapi.saemix, name.group=c("id"),
                        name.predictors=c("time","rapi"),name.response=c("rapi"),
                        name.covariates=c("gender"),
                        units=list(x="months",y="",covariates=c("")), verbose=FALSE)

xpl <- plotDiscreteData(saemix.data, outcome="count", which.cov="gender", breaks=c(0:9, 16, 25,80))
print(xpl)

if(saveFigs) {
  plot1 <- xpl + theme(text = element_text(size=rel(rsize.text)), legend.text=element_text(size=rel(rsize.text)), strip.text.x = element_text(size=rel(rsize.ticks)), strip.text.y = element_text(size=rel(rsize.ticks)))
  namfig<-"rapi_rawDataPropTime.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(xpl)
  dev.off()
}
# Save the plot 
plotRapi.rawDataPropTime <- xpl

# Simple histogram
hist(rapi.saemix$rapi, main="", xlab="RAPI score", breaks=30)

# Zooming on small values of scores
hist(rapi.saemix$rapi[rapi.saemix$rapi < 10], main="", xlab="RAPI score", breaks=30)

table(rapi.saemix$gender, as.integer(rapi.saemix$rapi > 2))

##############################################################################
################### Poisson model
##############################################################################
## Poisson with a time effect
# Model
count.poisson<-function(psi,id,xidep) { 
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  logp <- -lambda + y*log(lambda) - log(factorial(y))
  return(logp)
}
# Simulation function
countsimulate.poisson<-function(psi, id, xidep) {
  time<-xidep[,1]
  y<-xidep[,2]
  ymax<-max(y)
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  y<-rpois(length(time), lambda=lambda)
  y[y>ymax]<-ymax+1 # truncate to maximum observed value to avoid simulating aberrant values
  return(y)
}

## Poisson
### Model without covariate
saemix.model.poi<-saemixModel(model=count.poisson,description="Count model Poisson",simulate.function=countsimulate.poisson,
                              modeltype="likelihood",   
                              psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                              transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)), verbose=FALSE)

### Gender effect on intercept
saemix.model.poi.cov1<-saemixModel(model=count.poisson,description="Count model Poisson",simulate.function=countsimulate.poisson, 
                                   modeltype="likelihood",   
                                   psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                   transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)),
                                   covariance.model =matrix(data=1, ncol=2, nrow=2),
                                   covariate.model=matrix(c(1,0), ncol=2, byrow=TRUE), verbose=FALSE)

### Gender effect on intercept and slope
saemix.model.poi.cov2<-saemixModel(model=count.poisson,description="Count model Poisson",simulate.function=countsimulate.poisson, 
                                   modeltype="likelihood",   
                                   psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                   transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)),
                                   covariance.model =matrix(data=1, ncol=2, nrow=2),
                                   covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE), verbose=FALSE)

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, fim=FALSE, print=FALSE)

### Fit with saemix
poisson.fit<-saemix(saemix.model.poi,saemix.data,saemix.options)
poisson.fit.cov1<-saemix(saemix.model.poi.cov1,saemix.data,saemix.options)
poisson.fit.cov2<-saemix(saemix.model.poi.cov2,saemix.data,saemix.options)
#print(poisson.fit@results)
summary(poisson.fit.cov2)

### Results
if(FALSE) {
  cat("Poisson parameter at time 0 in base model: lambda_0=", exp(poisson.fit@results@fixed.effects[1]),"\n")
  cat("Poisson parameter at time 24 in base model: lambda_24=", exp(poisson.fit@results@fixed.effects[1]+24*poisson.fit@results@fixed.effects[2]),"\n")
}

# print(exp(poisson.fit@results@fixed.effects))
# exp(poisson.fit.cov2@results@fixed.effects)

### Simulations
nsim<-100
yfit1<-simulateDiscreteSaemix(poisson.fit.cov2, nsim=nsim)
plot1 <- discreteVPC(yfit1, outcome="count", which.cov="gender", breaks=c(0:9, 16, 25,80))
plot(plot1)

if(saveFigs) {
  namfig<-"rapi_poissonVPC.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot1)
  dev.off()
}
plotRapi.poissonVPC <-plot1

# Proportion of zeroes
#ggplot(yfit1@data@data, aes(x=rapi)) + geom_histogram()
#ggplot(yfit1@data@data, aes(x=rapi)) + geom_density() + geom_density(data=yfit1@sim.data@datasim,aes(x=ysim), colour="red") + xlim(0,10)

hist(yfit1@data@data$rapi, xlim=c(0,50), freq=F, breaks=50, xlab="Observed counts", main="")
lines(density(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim<50]), lwd = 2, col = 'red')

cat("Observed proportion of 0's", length(yfit1@data@data$rapi[yfit1@data@data$rapi==0])/yfit1@data@ntot.obs,"\n")
cat("      Poisson model, p=",length(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim==0])/length(yfit1@sim.data@datasim$ysim),"\n")

# Checking proportion of zeroes
yfit<-yfit1
simdat <-yfit@sim.data@datasim
simdat$time<-rep(yfit@data@data$time,nsim)
simdat$gender<-rep(yfit@data@data$gender,nsim)

if(FALSE) { # debug
  irep<-1
  xtab<-simdat[simdat$irep==irep,]
  regroupTable(xtab, value=3, c(4,5), fun=list(nev=function(x) sum(x==0), n=function(x) length(x)))
  xtab1 <- regroupTable(xtab, value="ysim", c("time","gender"), fun=list(nev=function(x) sum(x==0), n=function(x) length(x)))
}

ytab<-NULL
for(irep in 1:nsim) {
  xtab<-simdat[simdat$irep==irep,]
  xtab1 <- regroupTable(xtab, value="ysim", c("time","gender"), fun=list(nev=function(x) sum(x==0), n=function(x) length(x)))
  xtab1<-cbind(xtab1, freq=xtab1[,3]/xtab1[,4])
  ytab<-rbind(ytab,xtab1[,c(1,2,5)])
}

gtab <- regroupTable(ytab, value="freq", c("time","gender"), fun=list(lower=function(x) quantile(x, c(0.05)), median=function(x) quantile(x,0.5), upper=function(x) quantile(x,0.95)))
gtab$gender <- ifelse(gtab$gender==0,"Men","Women")
gtab$freq<-1
gtab1<-cbind(gtab, model="Poisson")

rapipl <- regroupTable(rapi.saemix,  value="rapi", c("time","gender"),fun=list(nev=function(x) sum(x==0), n=function(x) length(x)))
rapipl$freq<-rapipl$nev/rapipl$n
rapipl$sd <- sqrt((1-rapipl$nev/rapipl$n)*(rapipl$nev/rapipl$n**2))
rapipl$lower <- rapipl$freq-1.96*rapipl$sd
rapipl$upper <- rapipl$freq+1.96*rapipl$sd
rapipl$lower[rapipl$lower<0] <-0 # we should use a better approximation for CI
rapipl$gender <- ifelse(rapipl$gender==0,"Men","Women")

plot2 <- ggplot(rapipl, aes(x=time, y=freq, group=gender)) + geom_line() + geom_point() + 
  geom_line(data=gtab, aes(x=time, y=median,  group=gender), linetype=2, colour='lightblue') + 
  geom_ribbon(data=gtab,aes(ymin=lower, ymax=upper,  group=gender), alpha=0.5, fill='lightblue') +
  ylim(c(0,0.5)) + theme_bw() + theme(legend.position = "none") + facet_wrap(.~gender) +
  xlab("Time") + ylab("Proportion of subjects reporting no drinking episodes")

print(plot2)

# observed proportions at time 0 of non-drinking episodes is 41/347=11.8% in women and 46/471=9.7% in men
# clear underestimation in women especially

if(saveFigs) {
  namfig<-"rapi_poissonZeroesVPC.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot2)
  dev.off()
}
# Save the plot 
plotRapi.poissonZeroesVPC <- plot2

##############################################################################
## Zero-inflated Poisson model
##############################################################################
# Model
count.poissonzip<-function(psi,id,xidep) {
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  p0<-psi[id,3] # Probability of zero's
  lambda<- exp(intercept + slope*time)
  logp <- log(1-p0) -lambda + y*log(lambda) - log(factorial(y)) # Poisson
  logp0 <- log(p0+(1-p0)*exp(-lambda)) # Zeroes
  logp[y==0]<-logp0[y==0]
  return(logp)
}
# Simulation function
countsimulate.poissonzip<-function(psi, id, xidep) {
  time<-xidep[,1]
  y<-xidep[,2]
  ymax<-max(y)
  intercept<-psi[id,1]
  slope<-psi[id,2]
  p0<-psi[id,3] # Probability of zero's
  lambda<- exp(intercept + slope*time)
  prob0<-rbinom(length(time), size=1, prob=p0)
  y<-rpois(length(time), lambda=lambda)
  y[prob0==1]<-0
  y[y>ymax]<-ymax+1 # truncate to maximum observed value to avoid simulating aberrant values
  return(y)
}
## Generalized Poisson model with time effect
# Model
count.genpoisson<-function(psi,id,xidep) {
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  delta<-psi[id,3]
  logp <- log(lambda) + (y-1)*log(lambda+y*delta) - lambda - y*delta - log(factorial(y))
  return(logp)
}
# Simulation function - not given, see RNGforGPD

## Negative binomial model with time effect
# Model
count.NB<-function(psi,id,xidep) {
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  k<-psi[id,3]
  lambda<- exp(intercept + slope*time)
  logp <- log(factorial(y+k-1)) - log(factorial(y)) - log(factorial(k-1)) + y*log(lambda) - y*log(lambda+k) + k*log(k) - k*log(lambda+k)
  return(logp)
}
# Simulation function - not given, see RNGforGPD

saemix.model.zip<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                              simulate.function = countsimulate.poissonzip,
                              psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                              transform.par=c(0,0,3), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.3,0)), verbose=FALSE)

### ZIP Poisson with gender on both intercept
saemix.model.zip.cov1<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                                   simulate.function = countsimulate.poissonzip,
                                   psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                                   transform.par=c(0,0,3), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.3,0)),
                                   covariate.model = matrix(c(1,0,0),ncol=3, byrow=TRUE), verbose=FALSE)
### ZIP Poisson with gender on both intercept and slope
saemix.model.zip.cov2<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                                   simulate.function = countsimulate.poissonzip,
                                   psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                                   transform.par=c(0,0,3), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.3,0)),
                                   covariate.model = matrix(c(1,1,0),ncol=3, byrow=TRUE), verbose=FALSE)
### ZIP Poisson with gender on both intercept and slope, correlation between intercept and slope
covmat<-diag(c(1,1,0))
covmat[1,2]<-covmat[2,1]<-1
saemix.model.zip.cov3<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                                   simulate.function = countsimulate.poissonzip,
                                   psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                                   transform.par=c(0,0,3), covariance.model=covmat, omega.init=diag(c(0.5,0.3,0)),
                                   covariate.model = matrix(c(1,1,0),ncol=3, byrow=TRUE), verbose=FALSE)

zippoisson.fit<-saemix(saemix.model.zip,saemix.data,saemix.options)
zippoisson.fit.cov1<-saemix(saemix.model.zip.cov1,saemix.data,saemix.options)
zippoisson.fit.cov2<-saemix(saemix.model.zip.cov2,saemix.data,saemix.options)
zippoisson.fit.cov3<-saemix(saemix.model.zip.cov3,saemix.data,saemix.options)
print(zippoisson.fit.cov2@results)

exp(zippoisson.fit@results@fixed.effects)
exp(zippoisson.fit.cov1@results@fixed.effects)
exp(zippoisson.fit.cov2@results@fixed.effects)
exp(zippoisson.fit.cov3@results@fixed.effects)

ysim.zip2<-simulateDiscreteSaemix(zippoisson.fit.cov2,  100)
discreteVPC(ysim.zip2, outcome="count", breaks=c(0:9,16,25,80), which.cov="gender")

if(saveFigs) {
  namfig<-"rapi_diagnosZIP.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  discreteVPC(ysim.zip2, outcome="count", breaks=c(0:9,16,25,80), which.cov="gender")
  dev.off()
}

# Save the plot 
plotRapi.diagnosZIP <- discreteVPC(ysim.zip2, outcome="count", breaks=c(0:9,16,25,80), which.cov="gender")

### Proportion of zeroes
cat("Observed proportion of 0's", length(yfit1@data@data$rapi[yfit1@data@data$rapi==0])/yfit1@data@ntot.obs,"\n")
cat("      Poisson model, p=",length(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim==0])/length(yfit1@sim.data@datasim$ysim),"\n")
cat("  ZI-Poisson model, p=",length(ysim.zip2@sim.data@datasim$ysim[ysim.zip2@sim.data@datasim$ysim==0])/length(ysim.zip2@sim.data@datasim$ysim),"\n")

par(mfrow=c(1,3))
hist(yfit1@data@data$rapi, xlim=c(0,50), freq=F, breaks=30, xlab="Observed counts", main="")
hist(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim<50], xlim=c(0,50), freq=F, breaks=20, xlab="Simulated counts", main="Poisson model")
hist(ysim.zip2@sim.data@datasim$ysim[ysim.zip2@sim.data@datasim$ysim<50], xlim=c(0,50), freq=F, breaks=20, xlab="Simulated counts", main="ZIP model")

# Checking proportion of zeroes
yfit<-ysim.zip2
simdat <-yfit@sim.data@datasim
simdat$time<-rep(yfit@data@data$time,nsim)
simdat$gender<-rep(yfit@data@data$gender,nsim)

ytab<-NULL
for(irep in 1:nsim) {
  xtab<-simdat[simdat$irep==irep,]
  xtab1 <- regroupTable(xtab, value="ysim", c("time","gender"), fun=list(nev=function(x) sum(x==0), n=function(x) length(x)))
  xtab1<-cbind(xtab1, freq=xtab1[,3]/xtab1[,4])
  ytab<-rbind(ytab,xtab1[,c(1,2,5)])
}

gtab <- regroupTable(ytab, value="freq", c("time","gender"), fun=list(lower=function(x) quantile(x, c(0.05)), median=function(x) quantile(x,0.5), upper=function(x) quantile(x,0.95)))
gtab$gender <- ifelse(gtab$gender==0,"Men","Women")
gtab$freq<-1
gtab2<-cbind(gtab, model="ZIP")
gtab<-rbind(gtab1, gtab2)

plot2 <- ggplot(rapipl, aes(x=time, y=freq, group=gender)) + geom_line() + 
  geom_point() + 
  geom_line(data=gtab, aes(x=time, y=median,  group=gender), linetype=2, colour='lightblue') + 
  geom_ribbon(data=gtab,aes(ymin=lower, ymax=upper,  group=gender), alpha=0.5, fill='lightblue') +
  ylim(c(0,0.5)) + theme_bw() + theme(legend.position = "none") + facet_wrap(model~gender) +
  xlab("Time") + ylab("Proportion of drinking episodes")

print(plot2)

##############################################################################
## Hurdle - 2 models 
##############################################################################
saemix.data1<-saemixData(name.data=rapi.saemix[rapi.saemix$rapi>0,], name.group=c("id"),
                         name.predictors=c("time","rapi"),name.response=c("rapi"),
                         name.covariates=c("gender"),
                         units=list(x="week",y="",covariates=c("")), verbose=FALSE)

rapi.saemix$y0<-as.integer(rapi.saemix$rapi>0)
saemix.data0<-saemixData(name.data=rapi.saemix, name.group=c("id"),
                         name.predictors=c("time","y0"),name.response=c("y0"),
                         name.covariates=c("gender"),
                         units=list(x="week",y="",covariates=c("")), verbose=FALSE)

# Fit Binomial model to saemix.data0
binary.model<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  pobs = (y==0)*(1-pevent)+(y==1)*pevent
  logpdf <- log(pobs)
  return(logpdf)
}
# Associated simulation function
simulBinary<-function(psi,id,xidep) {
  tim<-xidep[,1]
  #  y<-xidep[,2] # not used for simulation
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  ysim<-rbinom(length(tim),size=1, prob=pevent)
  return(ysim)
}
saemix.hurdle0<-saemixModel(model=binary.model,description="Binary model",
                            modeltype="likelihood",simulate.function=simulBinary,
                            #                            psi0=matrix(c(-1.5,-.1,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                            psi0=matrix(c(0.5,-.1,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                            transform.par=c(0,0), covariate.model=c(1,1),
                            covariance.model =matrix(data=1, ncol=2, nrow=2), omega.init=diag(c(1,0.3)), verbose=FALSE)

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE, displayProgress=FALSE, print=FALSE)

hurdlefit0<-saemix(saemix.hurdle0,saemix.data0,saemix.options)

tab1<-table(rapi.saemix$rapi[rapi.saemix$time==0] == 0)  
tab.women <- table(rapi.saemix$rapi[rapi.saemix$time==0 & rapi.saemix$gender==0] == 0)  

if(FALSE) {
  cat("Observed proportion of 0's at time 0:",tab1[2]/sum(tab1),"\n") # 10.6%
  # Not comparable because of the IIV on theta1 => see below for a comparison using data simulated under the model
  cat("Expected proportion of 0's at time 0:",1-1/(1+exp(-hurdlefit0@results@fixed.effects[1])),"\n")
  cat("Observed proportion of 0's at time 0 in women:",tab.women[2]/sum(tab.women),"\n") # 9.7%
}

# Taking into account IIV on theta1 to compute expected proportion of zeroes
xcal<-rnorm(1000, mean=hurdlefit0@results@fixed.effects[1], sd=sqrt(hurdlefit0@results@omega[1,1]))
cat("Expected proportion of 0's at time 0 (in women):",mean(1-1/(1+exp(-xcal))),"\n")
cat("Observed proportion of 0's at time 0 in women:",tab.women[2]/sum(tab.women),"\n") # 9.7%

# Fit truncated Poisson model to saemix.data1

## Poisson with a time effect
# Model
truncated.poisson<-function(psi,id,xidep) { 
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  logp <- y*log(lambda) - log(factorial(y)) - log(exp(lambda)-1)
  return(logp)
}
# simulating from a truncated Poisson using the inverse function
# using method outlined by Keith Goldfield in post https://www.r-bloggers.com/2020/08/generating-data-from-a-truncated-distribution/
# Doesn't work, doesn't simulate any 1's
# rtruncPois <- function(n, range, lambda) {
#   F.a <- ppois(min(range), lambda=lambda)
#   F.b <- ppois(max(range), lambda=lambda)
#   u <- runif(n, min = F.a, max = F.b)
#   qpois(u, lambda=lambda)
# }

# Peter Dalgaard on R-help list: https://stat.ethz.ch/pipermail/r-help/2005-May/070680.html
rtpois <- function(N, lambda) {
  qpois(runif(N, dpois(0, lambda), 1), lambda)
}

# Simulation function
countsimulate.truncatedpoisson<-function(psi, id, xidep) {
  # left: truncate to 0
  # right: truncate to maximum observed value to avoid simulating aberrant values
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  #  y<-rtruncPois(length(time), y, lambda=lambda) 
  y<-rtpois(length(time), lambda=lambda)
  return(y)
}

saemix.hurdle1.cov2<-saemixModel(model=truncated.poisson,description="Count model Poisson",modeltype="likelihood",   
                                 simulate.function = countsimulate.truncatedpoisson,
                                 psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                 transform.par=c(0,0), omega.init=diag(c(0.5, 0.01)),
                                 covariance.model =matrix(data=1, ncol=2, nrow=2),
                                 covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE), verbose=FALSE)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, print=FALSE)

hurdlefit1<-saemix(saemix.hurdle1.cov2,saemix.data1,saemix.options)

summary(hurdlefit0)
summary(hurdlefit1)

# Parameters for the two components in a table form
# Compare to column B in Table 2 of Atkins
yfit0<-hurdlefit0
yfit1<-hurdlefit1

rr.tab<-data.frame(param=c("intercept", "beta.Male.inter", "slope", "beta.Male.slope", "omega.inter","omega.slope"), 
                   poissonNoZero=c(yfit1@results@fixed.effects, c(sqrt(diag(yfit1@results@omega)))),
                   logistic=c(yfit0@results@fixed.effects, c(sqrt(diag(yfit0@results@omega)))))
print(rr.tab)


# Simulate binary data
# proportion of 0's in the data
rapi.tab <- table(rapi.saemix$rapi == 0)

nsim<-100
ysim.hurdle0 <- simulateDiscreteSaemix(hurdlefit0, nsim=nsim)
cat("Observed proportion of 0's overall:",rapi.tab[2]/sum(rapi.tab),"\n")
cat("Simulated proportion of 0's overall:",sum(ysim.hurdle0@sim.data@datasim$ysim)/length(ysim.hurdle0@sim.data@datasim$ysim),"\n")

ysim.hurdle1 <- simulateDiscreteSaemix(hurdlefit1, nsim=nsim)

# Graph of proportion of 0's with time
yfit<-ysim.hurdle0
simdat <-yfit@sim.data@datasim
simdat$time<-rep(yfit@data@data$time,nsim)
simdat$gender<-rep(yfit@data@data$gender,nsim)

ytab<-NULL
for(irep in 1:nsim) {
  xtab<-simdat[simdat$irep==irep,]   # Here we want the proportion of zeroes, not 1's !
  xtab1 <- regroupTable(xtab, value="ysim", c("time","gender"), fun=list(nev=function(x) sum(1-x), n=function(x) length(x)))
  xtab1<-cbind(xtab1, freq=xtab1[,3]/xtab1[,4]) 
  ytab<-rbind(ytab,xtab1[,c(1,2,5)])
}

gtab3 <- regroupTable(ytab, value="freq", c("time","gender"), fun=list(lower=function(x) quantile(x, c(0.05)), median=function(x) quantile(x,0.5), upper=function(x) quantile(x,0.95)))
gtab3$gender <- ifelse(gtab3$gender==0,"Men","Women")
gtab3$freq<-1
gtab3<-cbind(gtab3, model="Hurdle")
gtab3$model <- factor(gtab3$model,levels=c("Poisson", "ZIP", "Hurdle"))

gtab<-rbind(gtab1, gtab2, gtab3)

# Compare proportion of zeroes across the 3 models
plot.prop0 <- ggplot(rapipl, aes(x=time, y=freq, group=gender)) + geom_line() + 
  geom_point() + 
  geom_line(data=gtab, aes(x=time, y=median,  group=gender), linetype=2, colour='lightblue') + 
  geom_ribbon(data=gtab,aes(ymin=lower, ymax=upper,  group=gender), alpha=0.5, fill='lightblue') +
  ylim(c(0,0.5)) + theme_bw() + theme(legend.position = "none") + facet_wrap(model~gender, ncol=2) +
  xlab("Time (months)") + ylab("Proportion of subjects without drinking episodes")

print(plot.prop0)
plot2 <- plot.prop0

if(saveFigs) {
  namfig<-"rapi_comparePropZeroes.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  print(plot2)
  dev.off()
}
# Save the plot 
plotRapi.comparePropZeroes <- plot2

nsim<-100
ysim.hurdle0 <- simulateDiscreteSaemix(hurdlefit0, nsim=nsim)
plot1 <- discreteVPC(ysim.hurdle0, outcome="count", which.cov="gender")
plot(plot1)

# Error with the first category with rtruncPois: values of 1 never simulated
# Done: different definition for the truncated Poisson distribution (rtpois)
ysim.hurdle1 <- simulateDiscreteSaemix(hurdlefit1, nsim=nsim)
plot1 <- discreteVPC(ysim.hurdle1, outcome="count", breaks=c(1:9, 16, 25, 80), which.cov="gender")

plot(plot1)

if(saveFigs) {
  namfig<-"rapi_hurdleVPC_truncatedPoiss.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot1)
  dev.off()
}

# Save the plot 
plotRapi.hurdleVPC.truncatedPoiss <- plot1


##############################################################################
#### Comparing the parameter estimates of the Poisson portion of the three models
##############################################################################
yfit1 <- poisson.fit.cov2
l1 <- c(yfit1@results@fixed.effects,NA, sqrt(diag(yfit1@results@omega)))
l1<-c(l1, yfit1@results@omega[1,2]/l1[6]/l1[7]) # rho

yfit2 <- zippoisson.fit.cov3
l2 <- c(yfit2@results@fixed.effects, sqrt(diag(yfit2@results@omega)[1:2]))
l2<-c(l2, yfit2@results@omega[1,2]/l2[6]/l2[7])

yfit3<-hurdlefit1
l3 <- c(yfit3@results@fixed.effects, NA,sqrt(diag(yfit3@results@omega)[1:2]))
l3<-c(l3, yfit3@results@omega[1,2]/l3[6]/l3[7])

rr.tab<-data.frame(Poisson=l1, ZIPoisson=l2, Hurdle=l3)

rownames(rr.tab)<-c("$\\alpha_0$ (-)","$\\beta_{male,\\alpha_0}$ (-)","$\\alpha_1$ (mo$^{(-1)}$)","$\\beta_{male,\\alpha_1}$ (-)","$p_0$ (-)", "$\\omega_{\\alpha_0}$ (-)","$\\omega_{\\alpha_1}$ (-)", "$\\rho(\\alpha_0,\\alpha_1)$ (-)")
for(icol in 1:3)
  rr.tab[,icol]<-format(rr.tab[,icol],digits=1, nsmall=0, trim=T, scientific=FALSE)
rr.tab[is.na(rr.tab)]<-"-"

##############################################################################
#### Standard errors of estimation
##############################################################################

# Time-patterns
time.pattern <- tapply(rapi.saemix$time, rapi.saemix$id, function(x) paste(x,collapse="-"))
table(time.pattern)

# Time-patterns by gender
print(table(time.pattern, rapi.saemix$gender[!duplicated(rapi.saemix$id)]))

###################################################
# Bootstrap approaches
###################################################
# Poisson with covariates
if(runBootstrap) {
  case.count <- saemix.bootstrap(zippoisson.fit.cov2, method="case", nboot=nboot) 
  cond.count <- saemix.bootstrap(zippoisson.fit.cov2, method="conditional", nboot=nboot) 
} else {
  case.count <- read.table(file.path(bootResDir,"bootstrapCase_rapiCov_Poisson.res"), header=T)
  cond.count <- read.table(file.path(bootResDir,"bootstrapCond_rapiCov_Poisson.res"), header=T)
  nboot<-dim(case.count)[1]
}
case.count <- case.count[!is.na(case.count[,2]),]
cond.count <- cond.count[!is.na(cond.count[,2]),]

corr<-poisson.fit.cov2@results@omega[1,2]/sqrt(poisson.fit.cov2@results@omega[1,1]*poisson.fit.cov2@results@omega[2,2])
par.estim<-c(poisson.fit.cov2@results@fixed.effects,diag(poisson.fit.cov2@results@omega), corr)
df2<-data.frame(parameter=colnames(case.count)[-c(1)], saemix=par.estim)
df2[dim(df2)[1],1]<-"corr.intSlope"
# compute correlation instead of covariance
ncol<-dim(case.count)[2]
case.count[ncol]<-case.count[ncol]/sqrt(case.count[ncol-1]*case.count[ncol-2])
cond.count[ncol]<-cond.count[ncol]/sqrt(cond.count[ncol-1]*cond.count[ncol-2])
for(i in 1:2) {
  if(i==1) {
    resboot1<-case.count
    namboot<-"case"
  } else {
    resboot1<-cond.count
    namboot <-"cNP"
  }
  mean.bootDist<-apply(resboot1, 2, mean)[-c(1)]
  sd.bootDist<-apply(resboot1, 2, sd)[-c(1)]
  quant.bootDist<-apply(resboot1[-c(1)], 2, quantile, c(0.025, 0.975))
  l1<-paste0(format(mean.bootDist, digits=2)," (",format(sd.bootDist,digits=2, trim=T),")")
  l2<-paste0("[",format(quant.bootDist[1,], digits=2),", ",format(quant.bootDist[2,],digits=2, trim=T),"]")
  df2<-cbind(df2, l1, l2)
  i1<-3+2*(i-1)
  colnames(df2)[i1:(i1+1)]<-paste0(namboot,".",c("estimate","CI"))
}
print(df2)
df.Poisson <- df2

# ZIPoisson with covariates
if(runBootstrap) {
  case.count <- saemix.bootstrap(zippoisson.fit.cov2, method="case", nboot=nboot) 
  cond.count <- saemix.bootstrap(zippoisson.fit.cov2, method="conditional", nboot=nboot) 
} else {
  # case.count <- read.table(file.path(saemixDir,"bootstrap","results","rapi_caseBootstrap.res"), header=T)
  # cond.count <- read.table(file.path(saemixDir,"bootstrap","results","rapi_condBootstrap.res"), header=T)
  case.count <- read.table(file.path(bootResDir,"bootstrapCase_rapiCov_ZIP_corr.res"), header=T)
  cond.count <- read.table(file.path(bootResDir,"bootstrapCond_rapiCov_ZIP_corr.res"), header=T)
  nboot<-dim(case.count)[1]
}
case.count <- case.count[!is.na(case.count[,2]),]
cond.count <- cond.count[!is.na(cond.count[,2]),]

corr<-zippoisson.fit.cov3@results@omega[1,2]/sqrt(zippoisson.fit.cov3@results@omega[1,1]*zippoisson.fit.cov3@results@omega[2,2])
par.estim<-c(zippoisson.fit.cov3@results@fixed.effects,diag(zippoisson.fit.cov3@results@omega)[zippoisson.fit.cov3@results@indx.omega], corr)
df2<-data.frame(parameter=colnames(case.count)[-c(1)], saemix=par.estim)
df2[dim(df2)[1],1]<-"corr.intSlope"
case.count[ncol]<-case.count[ncol]/sqrt(case.count[ncol-1]*case.count[ncol-2])
cond.count[ncol]<-cond.count[ncol]/sqrt(cond.count[ncol-1]*cond.count[ncol-2])
for(i in 1:2) {
  if(i==1) {
    resboot1<-case.count
    namboot<-"case"
  } else {
    resboot1<-cond.count
    namboot <-"cNP"
  }
  mean.bootDist<-apply(resboot1, 2, mean)[-c(1)]
  sd.bootDist<-apply(resboot1, 2, sd)[-c(1)]
  quant.bootDist<-apply(resboot1[-c(1)], 2, quantile, c(0.025, 0.975))
  l1<-paste0(format(mean.bootDist, digits=2)," (",format(sd.bootDist,digits=2, trim=T),")")
  l2<-paste0("[",format(quant.bootDist[1,], digits=2),", ",format(quant.bootDist[2,],digits=2, trim=T),"]")
  df2<-cbind(df2, l1, l2)
  i1<-3+2*(i-1)
  colnames(df2)[i1:(i1+1)]<-paste0(namboot,".",c("estimate","CI"))
}
print(df2)

df.Zip <- df2

# Hurdle model
## binary logistic component
if(runBootstrap) {
  case.count <- saemix.bootstrap(zippoisson.fit.cov2, method="case", nboot=nboot) 
  cond.count <- saemix.bootstrap(zippoisson.fit.cov2, method="conditional", nboot=nboot) 
} else {
  # case.count <- read.table(file.path(saemixDir,"bootstrap","results","rapi_caseBootstrap.res"), header=T)
  # cond.count <- read.table(file.path(saemixDir,"bootstrap","results","rapi_condBootstrap.res"), header=T)
  case.count <- read.table(file.path(bootResDir,"bootstrapCase_rapiCov_Hurdle0_corr.res"), header=T)
  cond.count <- read.table(file.path(bootResDir,"bootstrapCond_rapiCov_Hurdle0_corr.res"), header=T)
  nboot<-dim(case.count)[1]
}
case.count <- case.count[!is.na(case.count[,2]),]
cond.count <- cond.count[!is.na(cond.count[,2]),]


corr<-hurdlefit0@results@omega[1,2]/sqrt(hurdlefit0@results@omega[1,1]*hurdlefit0@results@omega[2,2])
par.estim<-c(hurdlefit0@results@fixed.effects,diag(hurdlefit0@results@omega)[hurdlefit0@results@indx.omega], corr)
df2<-data.frame(parameter=colnames(case.count)[-c(1)], saemix=par.estim)
df2[dim(df2)[1],1]<-"corr.intSlope"
case.count[ncol]<-case.count[ncol]/sqrt(case.count[ncol-1]*case.count[ncol-2])
cond.count[ncol]<-cond.count[ncol]/sqrt(cond.count[ncol-1]*cond.count[ncol-2])
for(i in 1:2) {
  if(i==1) {
    resboot1<-case.count
    namboot<-"case"
  } else {
    resboot1<-cond.count
    namboot <-"cNP"
  }
  mean.bootDist<-apply(resboot1, 2, mean)[-c(1)]
  sd.bootDist<-apply(resboot1, 2, sd)[-c(1)]
  quant.bootDist<-apply(resboot1[-c(1)], 2, quantile, c(0.025, 0.975))
  l1<-paste0(format(mean.bootDist, digits=2)," (",format(sd.bootDist,digits=2, trim=T),")")
  l2<-paste0("[",format(quant.bootDist[1,], digits=2),", ",format(quant.bootDist[2,],digits=2, trim=T),"]")
  df2<-cbind(df2, l1, l2)
  i1<-3+2*(i-1)
  colnames(df2)[i1:(i1+1)]<-paste0(namboot,".",c("estimate","CI"))
}
print(df2)

df.Hurdle0 <- df2

## truncated Poisson component

if(runBootstrap) {
  case.count <- saemix.bootstrap(zippoisson.fit.cov2, method="case", nboot=nboot) 
  cond.count <- saemix.bootstrap(zippoisson.fit.cov2, method="conditional", nboot=nboot) 
} else {
  # case.count <- read.table(file.path(saemixDir,"bootstrap","results","rapi_caseBootstrap.res"), header=T)
  # cond.count <- read.table(file.path(saemixDir,"bootstrap","results","rapi_condBootstrap.res"), header=T)
  case.count <- read.table(file.path(bootResDir,"bootstrapCase_rapiCov_Hurdle1_corr.res"), header=T)
  cond.count <- read.table(file.path(bootResDir,"bootstrapCond_rapiCov_Hurdle1_corr.res"), header=T)
  nboot<-dim(case.count)[1]
}
case.count <- case.count[!is.na(case.count[,2]),]
cond.count <- cond.count[!is.na(cond.count[,2]),]


corr<-hurdlefit1@results@omega[1,2]/sqrt(hurdlefit1@results@omega[1,1]*hurdlefit1@results@omega[2,2])
par.estim<-c(hurdlefit1@results@fixed.effects,diag(hurdlefit1@results@omega)[hurdlefit1@results@indx.omega], corr)
df2<-data.frame(parameter=colnames(case.count)[-c(1)], saemix=par.estim)
df2[dim(df2)[1],1]<-"corr.intSlope"
case.count[ncol]<-case.count[ncol]/sqrt(case.count[ncol-1]*case.count[ncol-2])
cond.count[ncol]<-cond.count[ncol]/sqrt(cond.count[ncol-1]*cond.count[ncol-2])
for(i in 1:2) {
  if(i==1) {
    resboot1<-case.count
    namboot<-"case"
  } else {
    resboot1<-cond.count
    namboot <-"cNP"
  }
  mean.bootDist<-apply(resboot1, 2, mean)[-c(1)]
  sd.bootDist<-apply(resboot1, 2, sd)[-c(1)]
  quant.bootDist<-apply(resboot1[-c(1)], 2, quantile, c(0.025, 0.975))
  l1<-paste0(format(mean.bootDist, digits=2)," (",format(sd.bootDist,digits=2, trim=T),")")
  l2<-paste0("[",format(quant.bootDist[1,], digits=2),", ",format(quant.bootDist[2,],digits=2, trim=T),"]")
  df2<-cbind(df2, l1, l2)
  i1<-3+2*(i-1)
  colnames(df2)[i1:(i1+1)]<-paste0(namboot,".",c("estimate","CI"))
}
print(df2)

df.Hurdle1 <- df2
