# Folders
workDir<-getwd() 

# @Eco
workDir<-"/home/eco/work/saemix/saemixextension/paperSaemix3"
saemixDir <- "/home/eco/work/saemix/saemixextension"
setwd(workDir)

# Libraries
library(saemix)

# Libraries needed to compute the FIM by AGQ
library(R6)
library(pracma)
library(compiler)
library(statmod)
library(Matrix)
library(randtoolbox)

# FIM by MC/AGQ (code S. Ueckert)

dirAGQ<-file.path(saemixDir,"fimAGQ")

# Bootstrap code
source(file.path(saemixDir, "bootstrap", "saemix_bootstrap.R"))

# Code to compute the exact FIM by MC/AGQ

# library(ggplot2)
# library(MASS)
# library(rlang)
# library(gridExtra)
library(tidyverse)

# Whether to save the plots
saveFigs<-FALSE
figDir <- getwd()

# Number of bootstrap samples
runBootstrap <- FALSE # to read the results from disk
nboot <-10

################################################### Data
data(rapi.saemix)

saemix.data<-saemixData(name.data=rapi.saemix, name.group=c("id"),
                        name.predictors=c("time","rapi"),name.response=c("rapi"),
                        name.covariates=c("gender"),
                        units=list(x="months",y="",covariates=c("")))

# Simple histogram
hist(rapi.saemix$rapi, main="", xlab="RAPI score", breaks=30)

# Zooming on small values of scores
hist(rapi.saemix$rapi[rapi.saemix$rapi < 10], main="", xlab="RAPI score", breaks=30)

table(rapi.saemix$gender, as.integer(rapi.saemix$rapi > 2))

################################################### Poisson model

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
                              transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)))

### Gender effect on intercept and slope
saemix.model.poi.cov2<-saemixModel(model=count.poisson,description="Count model Poisson",simulate.function=countsimulate.poisson, 
                                   modeltype="likelihood",   
                                   psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                   transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)),
                                   covariance.model =matrix(data=1, ncol=2, nrow=2),
                                   covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, fim=FALSE)

### Fit with saemix
poisson.fit<-saemix(saemix.model.poi,saemix.data,saemix.options)
poisson.fit.cov2<-saemix(saemix.model.poi.cov2,saemix.data,saemix.options)

### Results
if(FALSE) {
  cat("Poisson parameter at time 0 in base model: lambda_0=", exp(poisson.fit@results@fixed.effects[1]),"\n")
  cat("Poisson parameter at time 24 in base model: lambda_24=", exp(poisson.fit@results@fixed.effects[1]+24*poisson.fit@results@fixed.effects[2]),"\n")
}

### Simulations
nsim<-100
yfit1<-simulateDiscreteSaemix(poisson.fit.cov2, nsim=nsim)

hist(yfit1@data@data$rapi, xlim=c(0,50), freq=F, breaks=30, xlab="Observed counts", main="")
lines(density(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim<50]), lwd = 2, col = 'red')

cat("Observed proportion of 0's", length(yfit1@data@data$rapi[yfit1@data@data$rapi==0])/yfit1@data@ntot.obs,"\n")
cat("      Poisson model, p=",length(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim==0])/length(yfit1@sim.data@datasim$ysim),"\n")

# Checking proportion of zeroes
yfit<-yfit1
simdat <-yfit@sim.data@datasim
simdat$time<-rep(yfit@data@data$time,nsim)
simdat$gender<-rep(yfit@data@data$gender,nsim)

ytab<-NULL
for(irep in 1:nsim) {
  xtab<-simdat[simdat$irep==irep,]
  suppressMessages(
    xtab1 <- xtab %>%
      group_by(time, gender) %>%
      summarise(nev = sum(ysim==0), n=n()) %>%
      mutate(freq = nev/n)
  )
  ytab<-rbind(ytab,xtab1[,c("time","gender","freq")])
}
gtab <- ytab %>%
  group_by(time, gender) %>%
  summarise(lower=quantile(freq, c(0.05)), median=quantile(freq, c(0.5)), upper=quantile(freq, c(0.95))) %>%
  mutate(gender=ifelse(gender==0,"Men","Women"))
gtab$freq<-1
gtab1<-cbind(gtab, model="Poisson")

rapipl <- rapi.saemix %>%
  group_by(time, gender) %>%
  summarise(nev = sum(rapi==0), n=n()) %>%
  mutate(freq = nev/n, sd=sqrt((1-nev/n)/nev)) %>%
  mutate(lower=freq-1.96*sd, upper=freq+1.96*sd) 
rapipl$lower[rapipl$lower<0] <-0 # we should use a better approximation for CI

plot2 <- ggplot(rapipl, aes(x=time, y=freq, group=gender)) + geom_line() + 
  geom_point() + 
  geom_line(data=gtab, aes(x=time, y=median,  group=gender), linetype=2, colour='lightblue') + 
  geom_ribbon(data=gtab,aes(ymin=lower, ymax=upper,  group=gender), alpha=0.5, fill='lightblue') +
  ylim(c(0,0.5)) + theme_bw() + theme(legend.position = "none") + facet_wrap(.~gender) +
  xlab("Time") + ylab("Proportion of drinking episodes")

print(plot2)

################################################### ZIP Poisson model
## Zero-inflated Poisson model
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
# Simulation function - TBD, see RNGforGPD ?

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
# Simulation function - TBD ? 

## ZIP base model
saemix.model.zip<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                              simulate.function = countsimulate.poissonzip,
                              psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                              transform.par=c(0,0,3), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.3,0)))

### ZIP Poisson with gender on both intercept
saemix.model.zip.cov1<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                                   simulate.function = countsimulate.poissonzip,
                                   psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                                   transform.par=c(0,0,3), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.3,0)),
                                   covariate.model = matrix(c(1,0,0),ncol=3, byrow=TRUE))
### ZIP Poisson with gender on both intercept and slope
saemix.model.zip.cov2<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                                   simulate.function = countsimulate.poissonzip,
                                   psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                                   transform.par=c(0,0,3), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.3,0)),
                                   covariate.model = matrix(c(1,1,0),ncol=3, byrow=TRUE))

zippoisson.fit<-saemix(saemix.model.zip,saemix.data,saemix.options)
zippoisson.fit.cov1<-saemix(saemix.model.zip.cov1,saemix.data,saemix.options)
zippoisson.fit.cov2<-saemix(saemix.model.zip.cov2,saemix.data,saemix.options)

exp(zippoisson.fit@results@fixed.effects)
exp(zippoisson.fit.cov1@results@fixed.effects)
exp(zippoisson.fit.cov2@results@fixed.effects)

### Simulations
ysim.zip2<-simulateDiscreteSaemix(zippoisson.fit.cov2,  100)

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
  suppressMessages(
    xtab1 <- xtab %>%
      group_by(time, gender) %>%
      summarise(nev = sum(ysim==0), n=n()) %>%
      mutate(freq = nev/n)
  )
  ytab<-rbind(ytab,xtab1[,c("time","gender","freq")])
}
gtab <- ytab %>%
  group_by(time, gender) %>%
  summarise(lower=quantile(freq, c(0.05)), median=quantile(freq, c(0.5)), upper=quantile(freq, c(0.95))) %>%
  mutate(gender=ifelse(gender==0,"Men","Women"))
gtab$freq<-1
gtab2<-cbind(gtab, model="ZIP")
gtab<-rbind(gtab1, gtab2)

rapipl <- rapi.saemix %>%
  group_by(time, gender) %>%
  summarise(nev = sum(rapi==0), n=n()) %>%
  mutate(freq = nev/n, sd=sqrt((1-nev/n)/nev)) %>%
  mutate(lower=freq-1.96*sd, upper=freq+1.96*sd) 
rapipl$lower[rapipl$lower<0] <-0 # we should use a better approximation for CI

plot2 <- ggplot(rapipl, aes(x=time, y=freq, group=gender)) + geom_line() + 
  geom_point() + 
  geom_line(data=gtab, aes(x=time, y=median,  group=gender), linetype=2, colour='lightblue') + 
  geom_ribbon(data=gtab,aes(ymin=lower, ymax=upper,  group=gender), alpha=0.5, fill='lightblue') +
  ylim(c(0,0.5)) + theme_bw() + theme(legend.position = "none") + facet_wrap(model~gender) +
  xlab("Time") + ylab("Proportion of drinking episodes")

print(plot2)

# Grouping data by time and score
yfit <- ysim.zip2
ydat <- yfit@data
ysim <- yfit@sim.data@datasim$ysim
nsim<-length(ysim)/dim(ydat@data)[1]
obsmat<-data.frame(id=ydat@data[,ydat@name.group], x=ydat@data[,ydat@name.X], y=ydat@data[,ydat@name.response], covariate.group=ydat@data[,ydat@name.covariates])

# Regrouping times - not needed here as everyone has the same times

# Regrouping scores - observed data
mybreaks <- c(0:9, 16, 25, 80)
x <- cut(obsmat$y, breaks=mybreaks, include.lowest = TRUE)
obsmat$score.group <- x

# With tidyverse
counting.scores <- obsmat %>%
  group_by(x, covariate.group) %>%
  count(score.group)
number.samples <- obsmat %>%
  group_by(x, covariate.group) %>%
  summarise(n=n())
freq.scores <- number.samples %>%
  left_join(counting.scores, 
            by = c("x","covariate.group")) %>%
  mutate(freq=n.y/n.x) %>%
  mutate(covariate.group=ifelse(covariate.group==0,"Men","Women"))

# ggplot(data=freq.scores, aes(x=x, y=freq, group=covariate.group, colour=as.factor(covariate.group))) + geom_line() + facet_wrap(.~score.group, ncol=4)

# Regrouping scores - simulated data
ysim.tab <- data.frame(irep=rep(1:nsim, each=dim(ydat@data)[1]), x=rep(obsmat$x, nsim), covariate.group=rep(obsmat$covariate.group, nsim), score.group=cut(ysim, breaks=mybreaks, include.lowest = TRUE))
sim.scores <- ysim.tab %>%
  group_by(irep, x, covariate.group) %>%
  count(score.group)
simfreq.scores <- number.samples %>%
  left_join(sim.scores, 
            by = c("x","covariate.group")) %>%
  mutate(freq=n.y/n.x) %>%
  mutate(covariate.group=ifelse(covariate.group==0,"Men","Women"))

simfreq.bands <- simfreq.scores %>%
  group_by(x, covariate.group, score.group) %>%
  summarise(lower=quantile(freq, c(0.05)), median=quantile(freq, c(0.5)), upper=quantile(freq, c(0.95)), freq=mean(freq)) 

plot.counts <- ggplot(data=freq.scores, aes(x=x, y=freq, group=covariate.group, colour=as.factor(covariate.group))) + geom_line() + 
  geom_line(data=simfreq.bands, aes(x=x, y=median, group=covariate.group, colour=as.factor(covariate.group)), linetype="dashed") + 
  geom_ribbon(data=simfreq.bands, aes(ymin=lower, ymax=upper,  group=covariate.group, fill=as.factor(covariate.group)), alpha=0.2) + 
  xlab("Time") + ylab("Proportion of counts") + guides(fill=guide_legend(title='Gender'), colour=guide_legend(title='Gender')) +
  facet_wrap(.~score.group, ncol=4)

print(plot.counts)


################################################### Hurdle model
## Hurdle - 2 models 
saemix.data1<-saemixData(name.data=rapi.saemix[rapi.saemix$rapi>0,], name.group=c("id"),
                         name.predictors=c("time","rapi"),name.response=c("rapi"),
                         name.covariates=c("gender"),
                         units=list(x="week",y="",covariates=c("")))

rapi.saemix$y0<-as.integer(rapi.saemix$rapi==0)
saemix.data0<-saemixData(name.data=rapi.saemix, name.group=c("id"),
                         name.predictors=c("time","y0"),name.response=c("y0"),
                         name.covariates=c("gender"),
                         units=list(x="week",y="",covariates=c("")))

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
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  ysim<-rbinom(length(tim),size=1, prob=pevent)
  return(ysim)
}
saemix.hurdle0<-saemixModel(model=binary.model,description="Binary model",
                            modeltype="likelihood",simulate.function=simulBinary,
                            psi0=matrix(c(-1.5,-.1,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                            transform.par=c(0,0), covariate.model=c(1,1),
                            covariance.model=matrix(c(1,0,0,1),ncol=2), omega.init=diag(c(1,0.3)))

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE, displayProgress=FALSE)

hurdlefit0<-saemix(saemix.hurdle0,saemix.data0,saemix.options)
cat("Expected proportion of 0's at time 0:",1/(1+exp(-hurdlefit0@results@fixed.effects[1])),"\n")
table(rapi.saemix$rapi[rapi.saemix$time==0] == 0) # 10.6%

# Fit Poisson model to saemix.data1
saemix.hurdle1.cov2<-saemixModel(model=count.poisson,description="Count model Poisson",modeltype="likelihood",   
                                 simulate.function = countsimulate.poisson,
                                 psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                 transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)),
                                 covariance.model =matrix(data=1, ncol=2, nrow=2),
                                 covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE))
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)

hurdlefit1<-saemix(saemix.hurdle1.cov2,saemix.data1,saemix.options)

summary(hurdlefit0)
summary(hurdlefit1)

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
  xtab<-simdat[simdat$irep==irep,]
  suppressMessages(
    xtab1 <- xtab %>%
      group_by(time, gender) %>%
      summarise(nev = sum(ysim), n=n()) %>%
      mutate(freq = nev/n)
  )
  ytab<-rbind(ytab,xtab1[,c("time","gender","freq")])
}
gtab <- ytab %>%
  group_by(time, gender) %>%
  summarise(lower=quantile(freq, c(0.05)), median=quantile(freq, c(0.5)), upper=quantile(freq, c(0.95))) %>%
  mutate(gender=ifelse(gender==0,"Men","Women"))
gtab$freq<-1
gtab3<-cbind(gtab, model="Hurdle")
gtab<-rbind(gtab1, gtab2, gtab3)
gtab <- gtab %>%
  mutate(model=factor(model, levels=c("Poisson", "ZIP", "Hurdle")))

rapipl <- rapi.saemix %>%
  group_by(time, gender) %>%
  summarise(nev = sum(y0), n=n()) %>%
  mutate(freq = nev/n, sd=sqrt((1-nev/n)/nev)) %>%
  mutate(lower=freq-1.96*sd, upper=freq+1.96*sd) 
rapipl$lower[rapipl$lower<0] <-0 # we should use a better approximation for CI

# Table form - compare to column B in Table 2
yfit0<-hurdlefit0
yfit1<-hurdlefit1

rr.tab<-data.frame(param=c("intercept", "beta.Male.inter", "slope", "beta.Male.slope", "omega.inter","omega.slope"), 
                   poissonNoZero=c(yfit1@results@fixed.effects, c(sqrt(diag(yfit1@results@omega)))),
                   logistic=c(yfit0@results@fixed.effects, c(sqrt(diag(yfit0@results@omega)))))

print(rr.tab)

# Comparing P(Y=0) in the different models
plot.prop0 <- ggplot(rapipl, aes(x=time, y=freq, group=gender)) + geom_line() + 
  geom_point() + 
  geom_line(data=gtab, aes(x=time, y=median,  group=gender), linetype=2, colour='lightblue') + 
  geom_ribbon(data=gtab,aes(ymin=lower, ymax=upper,  group=gender), alpha=0.5, fill='lightblue') +
  ylim(c(0,0.5)) + theme_bw() + theme(legend.position = "none") + facet_wrap(model~gender, ncol=2) +
  xlab("Time") + ylab("Proportion of subjects without drinking episodes")

print(plot.prop0)

nsim<-100
ysim.hurdle1 <- simulateDiscreteSaemix(hurdlefit1, nsim=nsim)

# Grouping data by time and score
yfit <- ysim.hurdle1
ydat <- yfit@data
ysim <- yfit@sim.data@datasim$ysim
nsim<-length(ysim)/dim(ydat@data)[1]
obsmat<-data.frame(id=ydat@data[,ydat@name.group], x=ydat@data[,ydat@name.X], y=ydat@data[,ydat@name.response], covariate.group=ydat@data[,ydat@name.covariates])

# Regrouping times - not needed here as everyone has the same times

# Regrouping scores - observed data
mybreaks <- c(0:9, 16, 25, 80)
x <- cut(obsmat$y, breaks=mybreaks, include.lowest = TRUE)
obsmat$score.group <- x

# With tidyverse
counting.scores <- obsmat %>%
  group_by(x, covariate.group) %>%
  count(score.group)
number.samples <- obsmat %>%
  group_by(x, covariate.group) %>%
  summarise(n=n())
freq.scores <- number.samples %>%
  left_join(counting.scores, 
            by = c("x","covariate.group")) %>%
  mutate(freq=n.y/n.x) %>%
  mutate(covariate.group=ifelse(covariate.group==0,"Men","Women"))

# ggplot(data=freq.scores, aes(x=x, y=freq, group=covariate.group, colour=as.factor(covariate.group))) + geom_line() + facet_wrap(.~score.group, ncol=4)

# Regrouping scores - simulated data
ysim.tab <- data.frame(irep=rep(1:nsim, each=dim(ydat@data)[1]), x=rep(obsmat$x, nsim), covariate.group=rep(obsmat$covariate.group, nsim), score.group=cut(ysim, breaks=mybreaks, include.lowest = TRUE))
sim.scores <- ysim.tab %>%
  group_by(irep, x, covariate.group) %>%
  count(score.group)
simfreq.scores <- number.samples %>%
  left_join(sim.scores, 
            by = c("x","covariate.group")) %>%
  mutate(freq=n.y/n.x) %>%
  mutate(covariate.group=ifelse(covariate.group==0,"Men","Women"))

simfreq.bands <- simfreq.scores %>%
  group_by(x, covariate.group, score.group) %>%
  summarise(lower=quantile(freq, c(0.05)), median=quantile(freq, c(0.5)), upper=quantile(freq, c(0.95)), freq=mean(freq)) 

plot.counts <- ggplot(data=freq.scores, aes(x=x, y=freq, group=covariate.group, colour=as.factor(covariate.group))) + geom_line() + 
  geom_line(data=simfreq.bands, aes(x=x, y=median, group=covariate.group, colour=as.factor(covariate.group)), linetype="dashed") + 
  geom_ribbon(data=simfreq.bands, aes(ymin=lower, ymax=upper,  group=covariate.group, fill=as.factor(covariate.group)), alpha=0.2) + 
  xlab("Time") + ylab("Proportion of counts") + guides(fill=guide_legend(title='Gender'), colour=guide_legend(title='Gender')) +
  facet_wrap(.~score.group, ncol=4)

print(plot.counts)

################################################### Estimation errors (ZIP Poisson model)
# Bootstrap approaches

if(runBootstrap) {
  case.count <- saemix.bootstrap(zippoisson.fit.cov2, method="case", nboot=nboot) 
  cond.count <- saemix.bootstrap(zippoisson.fit.cov2, method="conditional", nboot=nboot) 
} else {
  case.count <- read.table(file.path(saemixDir,"bootstrap","results","rapi_caseBootstrap.res"), header=T)
  cond.count <- read.table(file.path(saemixDir,"bootstrap","results","rapi_condBootstrap.res"), header=T)
  nboot<-dim(case.count)[1]
}
case.count <- case.count[!is.na(case.count[,2]),]
cond.count <- cond.count[!is.na(cond.count[,2]),]

par.estim<-c(zippoisson.fit.cov2@results@fixed.effects,diag(zippoisson.fit.cov2@results@omega)[zippoisson.fit.cov2@results@indx.omega])
df2<-data.frame(parameter=colnames(case.count)[-c(1)], saemix=par.estim)
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
