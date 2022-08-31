
# @Eco
workDir<-"/home/eco/work/saemix/saemixextension/bootstrap"
saemixDir <- "/home/eco/work/saemix/saemixextension"
setwd(workDir)

library(tidyverse)

#  Code saemix
progDir <- file.path(saemixDir,"R")
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDir,"SaemixData.R"))
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixModel.R"))
source(file.path(progDir,"SaemixObject.R"))
source(file.path(progDir,"main.R"))
source(file.path(progDir,"func_aux.R"))
source(file.path(progDir,"main_initialiseMainAlgo.R"))
source(file.path(progDir,"main_estep.R"))
source(file.path(progDir,"main_mstep.R"))
source(file.path(progDir,"func_FIM.R"))
source(file.path(progDir,"func_npde.R"))
source(file.path(progDir,"func_plots.R"))
source(file.path(progDir,"func_distcond.R"))
source(file.path(progDir,"func_simulations.R"))
source(file.path(progDir,"compute_LL.R"))
source(file.path(progDir,"func_estimParam.R"))
source(file.path(progDir,"backward.R"))
source(file.path(progDir,"forward.R"))
source(file.path(progDir,"stepwise.R"))
source(file.path(progDir,"func_stepwise.R"))
source(file.path(progDir,"func_compare.R"))

# Bootstrap code
source(file.path(saemixDir, "bootstrap", "saemix_bootstrap.R"))

# Number of bootstrap samples
nboot <- 500

# Data
datDir <- file.path(saemixDir, "data")
toenail.saemix <- read.table(file.path(datDir, "toenail.saemix.tab"), header = T)

saemix.data<-saemixData(name.data=toenail.saemix,name.group=c("id"),name.predictors=c("time","y"), name.response="y",
                        name.covariates=c("treatment"))

# Model function
binary.model<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  logpdf<-rep(0,length(tim))
  P.obs = (y==0)*(1-pevent)+(y==1)*pevent
  logpdf <- log(P.obs)
  return(logpdf)
}

# Simulation function
simulBinary<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-1/(1+exp(-logit))
  ysim<-rbinom(length(tim),size=1, prob=pevent)
  return(ysim)
}


saemix.model<-saemixModel(model=binary.model, description="Binary model", simulate.function = simulBinary,
                          modeltype="likelihood",
                          psi0=matrix(c(-0.5,-.15,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                          transform.par=c(0,0), covariate.model=c(0,1),covariance.model=matrix(c(1,0,0,0),ncol=2), omega.init=diag(c(0.5,0.3)))
saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE)
binary.fit<-saemix(saemix.model,saemix.data,saemix.options)

if(FALSE)
  case.bin <- saemix.bootstrap(binary.fit, method="case", nboot=nboot) 

table(binary.fit@data@data$y)

data.boot.case <- dataGen.case(binary.fit)

table(data.boot.case@data$y)
table(binary.fit@data@data$y)
fit.boot.case<-saemix(binary.fit@model, data.boot.case, list(map=F, save.graphs=F, fim=F, ll.is=F))

nsamp<-100
binary.fit<-conddist.saemix(binary.fit, nsamp=nsamp) # estimate conditional distributions and sample residuals
eta.sampc<-centerDist.NPcond(binary.fit, nsamp=nsamp)
data.boot.cond <- dataGen.NP(binary.fit, nsamp=nsamp,eta.sampc=eta.sampc, conditional=TRUE)
# checking that the good results for cNP aren't because the original y.1 column was used => no, identical results when removing that column
data.boot.cond@data<-data.boot.cond@data[,colnames(data.boot.cond@data)!="y.1"]
fit.boot.cond<-saemix(binary.fit@model, data.boot.cond, list(map=F, save.graphs=F, fim=F, ll.is=F))

table(binary.fit@data@data$y)
table(data.boot.case@data$y)
table(data.boot.cond@data$y)

table(binary.fit@data@data$y.1)
table(data.boot.case@data$y.1)
table(data.boot.cond@data$y.1)

fit.boot.case@results@fixed.effects
fit.boot.cond@results@fixed.effects

# No obvious discrepancy between datasets but estimates very different
# Plotting 

table(cut(data.boot.case@data$time, breaks=c(-1,0.25,1.25,2.25, 3.25, 7,10,20)))
data.boot.case@data$visit<-cut(data.boot.case@data$time, breaks=c(-1,0.25,1.25,2.25, 3.25, 7,10,20))
data.boot.cond@data$visit<-cut(data.boot.cond@data$time, breaks=c(-1,0.25,1.25,2.25, 3.25, 7,10,20))
orig.data<-binary.fit@data
orig.data@data$visit<-cut(orig.data@data$time, breaks=c(-1,0.25,1.25,2.25, 3.25, 7,10,20))

# Explore data - again very similar and both look like the original
toe1 <- orig.data@data %>%
  group_by(visit, treatment) %>%
  summarise(nev = sum(y), n=n()) %>%
  mutate(freq = nev/n, sd=sqrt((1-nev/n)/nev)) %>%
  mutate(lower=freq-1.96*sd, upper=freq+1.96*sd)
toe1$lower[toe1$lower<0] <-0 # we should use a better approximation for CI
toe1$treatment <- factor(toe1$treatment, labels=c("A","B"))
ypl <- cbind(toe1, sample="Original")

toe1 <- data.boot.case@data %>%
  group_by(visit, treatment) %>%
  summarise(nev = sum(y), n=n()) %>%
  mutate(freq = nev/n, sd=sqrt((1-nev/n)/nev)) %>%
  mutate(lower=freq-1.96*sd, upper=freq+1.96*sd)
toe1$lower[toe1$lower<0] <-0 # we should use a better approximation for CI
toe1$treatment <- factor(toe1$treatment, labels=c("A","B"))
ypl <- rbind(ypl, cbind(toe1, sample="Case"))

toe1 <- data.boot.cond@data %>%
  group_by(visit, treatment) %>%
  summarise(nev = sum(y), n=n()) %>%
  mutate(freq = nev/n, sd=sqrt((1-nev/n)/nev)) %>%
  mutate(lower=freq-1.96*sd, upper=freq+1.96*sd)
toe1$lower[toe1$lower<0] <-0 # we should use a better approximation for CI
toe1$treatment <- factor(toe1$treatment, labels=c("A","B"))
ypl <- rbind(ypl, cbind(toe1, sample="cNP"))

plot1<-ggplot(ypl, aes(x=visit, y=freq, group=treatment)) + geom_line(aes(colour=treatment)) + 
  geom_point(aes(colour=treatment)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=treatment), alpha=0.2) +
  ylim(c(0,1)) + theme_bw() + theme(legend.position = "top") +
  xlab("Visit number") + ylab("Observed frequency of infection") + facet_wrap(.~sample)

print(plot1)
