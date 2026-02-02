#####################################################################################
# Libraries and folders

workDir<-"/home/eco/work/saemix/saemixextension/paperSaemix3"
saemixDir <- "/home/eco/work/saemix/saemixextension"
#figDir<-file.path(saemixDir,"documentation","figs")
setwd(workDir)

# Libraries
library(saemix)

#library(tidyverse)

#####################################################################################
################ Data
data(rapi.saemix)
rapi.saemix$gender <- ifelse(rapi.saemix$gender=="Men",1,0) # Female=reference class as in Atkins

# Hurdel model setup
saemix.data1<-saemixData(name.data=rapi.saemix[rapi.saemix$rapi>0,], name.group=c("id"),
                         name.predictors=c("time","rapi"),name.response=c("rapi"),
                         name.covariates=c("gender"),
                         units=list(x="week",y="",covariates=c("")), verbose=FALSE)

rapi.saemix$y0<-as.integer(rapi.saemix$rapi>0)
saemix.data0<-saemixData(name.data=rapi.saemix, name.group=c("id"),
                         name.predictors=c("time","y0"),name.response=c("y0"),
                         name.covariates=c("gender"),
                         units=list(x="week",y="",covariates=c("")), verbose=FALSE)

################ Binomial model for P(count=0)
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

# Taking into account IIV on theta1 to compute expected proportion of zeroes
xcal<-rnorm(1000, mean=hurdlefit0@results@fixed.effects[1], sd=sqrt(hurdlefit0@results@omega[1,1]))
cat("Expected proportion of 0's at time 0 (in women):",mean(1-1/(1+exp(-xcal))),"\n")
cat("Observed proportion of 0's at time 0 in women:",tab.women[2]/sum(tab.women),"\n") # 9.7%

################ Hurdle model for count>0
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
rtruncPois <- function(n, range, lambda) {
  F.a <- ppois(min(range), lambda=lambda)
  F.b <- ppois(max(range), lambda=lambda)
  u <- runif(n, min = F.a, max = F.b)
  qpois(u, lambda=lambda)
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
  y<-rtruncPois(length(time), y, lambda=lambda) 
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
nsim<-100
ysim.hurdle1 <- simulateDiscreteSaemix(hurdlefit1, nsim=nsim)

ysim.hurdle0 <- simulateDiscreteSaemix(hurdlefit0, nsim=nsim)

#####################################################################################
################ Debug VPC for hurdle model

# Error with the first category, maybe check if starts at zero in discreteVPC.aux (TODO)
## symptom: category 1 not well predicted
plot1 <- discreteVPC(ysim.hurdle1, outcome="count", breaks=c(1:9, 16, 25, 80), which.cov="gender")

# Simulated data: no 1's simulated in the model
ysim <- ysim.hurdle1@sim.data@datasim$ysim
table(ysim)

# discreteVPC
object <- ysim.hurdle1
outcome <- "categorical"
if(object@sim.data@nsim==0) cat("Missing simulated data\n")
verbose <- TRUE

# discreteVPCcount
max.cat<-10
# breaks=NULL, catlabel=NULL
x<-discreteVPC.aux(object, max.cat=max.cat, verbose=verbose)

# So the problem is in the simulations, not the VPC

nsim<-100
seed <- 123456
predictions <- TRUE

# using rtruncPois
y1<-rtruncPois(10000,c(1,74),2)
y0<-rtruncPois(10000,c(0,74),2)
y2<-rtruncPois(10000, c(1,75),3)-1

# Checking rtruncPois
x0 <- rpois(10000,2)
x1 <- rtpois(10000,2)

gtab <- rbind(data.frame(x=y1, range="1-74"), data.frame(x=y0, range="0-74"), data.frame(x=y2, range="(1-75)-1"), data.frame(x=x0, range="rpois full"), data.frame(x=x0[x0>0], range="rpois trunc"), data.frame(x=x1, range="rtpois"))
ggplot(gtab,aes(x=x, fill=range)) + geom_histogram(position="dodge", bins=20) + lims(x=c(0,9))

# Other definitions
rtpois <- function(N, lambda)
  qpois(runif(N, dpois(0, lambda), 1), lambda)
