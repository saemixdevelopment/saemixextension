########################################################################
### Create datasets forn the new discrete outcome models
### Prepare figures for the LaTeX documentation
########################################################################
# Folders
saemixDir<-"/home/eco/work/saemix/saemixextension"
datDir<-file.path(saemixDir, "data")
docDir<-file.path(saemixDir, "documentation")
figDir<-file.path(docDir,"figs")
progDir<-file.path(saemixDir,"R")

createDataSaemix<-FALSE
testMode<-FALSE
########################################################################
# Load library
library(ggplot2)
library(gridExtra)
library(tidyverse)

if(!testMode) {
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
} else
  library(saemix)

######################################################################## BINARY
# Binary data - Toenail 
# library(HSAUR3)
# data(toenail)
# toe<-transform(toenail,y=as.integer(toenail$outcome=="moderate or severe"))
# saemix.data<-saemixData(name.data=toe,name.group=c("patientID"),name.predictors=c("time","y"), name.response="y",
#                        name.covariates=c("treatment"),name.X=c("time"))

if(createDataSaemix) {
  library(prLogistic)
  data(Toenail)
  toenail.saemix<-Toenail
  colnames(toenail.saemix)<-c("id","y","treatment","time","visit")
  toenail.saemix<-toenail.saemix[,c(1,4,2,3,5)]
  write.table(toenail.saemix,file.path(datDir,"toenail.saemix.tab"), quote=F, row.names=F)
}
if(testMode)
  data(toenail.saemix) else
  toenail.saemix<-read.table(file.path(datDir,"toenail.saemix.tab"),header=TRUE)
toe<-toenail.saemix

# Data
saemix.data<-saemixData(name.data=toenail.saemix,name.group=c("id"),name.predictors=c("time","y"), name.response="y",
                        name.covariates=c("treatment"),name.X=c("time"))

# Explore data
toe1 <- toenail.saemix %>%
  group_by(visit, treatment) %>%
  summarise(nev = sum(y), n=n()) %>%
  mutate(freq = nev/n, sd=sqrt((1-nev/n)/nev)) %>%
  mutate(lower=freq-1.96*sd, upper=freq+1.96*sd)
toe1$lower[toe1$lower<0] <-0 # we should use a better approximation for CI
toe1$treatment <- factor(toe1$treatment, labels=c("A","B"))


plot1<-ggplot(toe1, aes(x=visit, y=freq, group=treatment)) + geom_line(aes(colour=treatment)) + 
  geom_point(aes(colour=treatment)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=treatment), alpha=0.2) +
  ylim(c(0,1)) + theme_bw() + theme(legend.position = "top") +
  xlab("Visit number") + ylab("Observed frequency of infection")

# theme(legend.position = "none") + scale_x_discrete(name = "Treatment", labels = c("A", "B")) +
#labs(colour = "Treatment") +

namfig<-"toenail_rawdata.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(plot1)
dev.off()

plot2 <- ggplot(toe1, aes(x=visit, y=nev, group=treatment, fill=treatment)) + 
  geom_bar(stat="identity",position=position_dodge()) +
  theme_bw() + theme(legend.position = "top") +
  xlab("Visit number") + ylab("Number of infected subjects")

namfig<-"toenail_barplotData.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(plot2)
dev.off()


if(FALSE) {
  sum(toe$y[toe$visit==1])/length(toe$y[toe$visit==1])
  sum(toe$y[toe$visit==1 & toe$treatment==1])/length(toe$y[toe$visit==1 & toe$treatment==1])
  sum(toe$y[toe$visit==1 & toe$treatment==0])/length(toe$y[toe$visit==1 & toe$treatment==0])
  sum(toe$y[toe$visit==7])/length(toe$y[toe$visit==7])
}

# Fit with saemix
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

saemix.model<-saemixModel(model=binary.model,description="Binary model",
                          modeltype="likelihood",
                          psi0=matrix(c(0,-.5,0,0.5),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                          transform.par=c(0,0), covariate.model=c(0,1),covariance.model=matrix(c(1,0,0,1),ncol=2))
saemix.model.paper<-saemixModel(model=binary.model,description="Binary model",
                          modeltype="likelihood",
                          psi0=matrix(c(-5,-.1,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                          transform.par=c(0,0), covariate.model=c(0,1),covariance.model=matrix(c(1,0,0,0),ncol=2))
saemix.model.paper2<-saemixModel(model=binary.model,description="Binary model",
                                modeltype="likelihood",
                                psi0=matrix(c(-5,-0.2,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                                transform.par=c(0,0), covariate.model=c(0,1),covariance.model=matrix(c(1,0,0,1),ncol=2))

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE)

binary.fit<-saemix(saemix.model,saemix.data,saemix.options)
binary.fit.paper<-saemix(saemix.model.paper,saemix.data,saemix.options)
binary.fit.paper2<-saemix(saemix.model.paper2,saemix.data,saemix.options)

binary.fit.paper@results@fixed.effects
sqrt(binary.fit.paper@results@omega[1,1])
# Estimates with Monolix
# $\theta_1$ & 1.76 & 0.33 & 18.6 \\
# $\theta_2$ & 0.36 & 0.039 & 11 \\
# $\beta$ & 0.19 & 0.064 & 33.3\\
# $\omega_1$ & 4.05 & 0.37 & 9.04\\

# Simulations
## Simulation file
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

yfit<-binary.fit.paper
nsim<-100
yfit <- simulateDiscreteSaemix(yfit, simulBinary, nsim=nsim)

yfit<-binary.fit.paper2
nsim<-100
yfit <- simulateDiscreteSaemix(yfit, simulBinary, nsim=nsim)

#yfit1<-simulate.SaemixObject(yfit, nsim=nsim, predictions=FALSE, uncertainty=uncertainty)

## Population-level prediction interval for the proportion, per visit across treatments
head(yfit@sim.data@datasim)
simdat <-yfit@sim.data@datasim
simdat$visit<-rep(toenail.saemix$visit,nsim)
simdat$treatment<-rep(toenail.saemix$treatment,nsim)

ytab<-NULL
for(irep in 1:nsim) {
  xtab<-simdat[simdat$irep==irep,]
  xtab1 <- xtab %>%
    group_by(visit) %>%
    summarise(nev = sum(ysim), n=n()) %>%
    mutate(freq = nev/n)
  ytab<-rbind(ytab,xtab1[,c(1,4)])
}
gtab<-data.frame(visit=sort(unique(ytab$visit)),matrix(unlist(tapply(ytab$freq, ytab$visit, quantile, c(0.05,0.5, 0.95))), byrow=T, ncol=3))
colnames(gtab)[2:4]<-c("lower","median","upper")
gtab$treatment<-1
gtab$freq<-1

plot1 <- ggplot(toe1, aes(x=visit, y=freq, group=treatment)) + geom_line(aes(colour=treatment)) + 
  geom_point(aes(colour=treatment)) + 
  geom_line(data=gtab, aes(x=visit, y=median), linetype=2, colour='lightblue') + 
  geom_ribbon(data=gtab,aes(ymin=lower, ymax=upper), alpha=0.5, fill='lightblue') +
  ylim(c(0,0.5)) + theme_bw() + theme(legend.position = "top") +
  xlab("Visit number") + ylab("Observed frequency of infection")

plot1

namfig<-"toenail_globalpropVPC.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(plot1)
dev.off()


ytab<-NULL
for(irep in 1:nsim) {
  xtab<-simdat[simdat$irep==irep,]
  xtab1 <- xtab %>%
    group_by(visit, treatment) %>%
    summarise(nev = sum(ysim), n=n()) %>%
    mutate(freq = nev/n)
  ytab<-rbind(ytab,xtab1[,c("visit","freq","treatment")])
}
gtab <- ytab %>%
  group_by(visit, treatment) %>%
  summarise(lower=quantile(freq, c(0.05)), median=quantile(freq, c(0.5)), upper=quantile(freq, c(0.95))) %>%
  mutate(treatment=ifelse(treatment==1,"B","A"))
gtab$freq<-1

plot2 <- ggplot(toe1, aes(x=visit, y=freq, group=treatment)) + geom_line(aes(colour=treatment)) + 
  geom_point(aes(colour=treatment)) + 
  geom_line(data=gtab, aes(x=visit, y=median), linetype=2, colour='lightblue') + 
  geom_ribbon(data=gtab,aes(ymin=lower, ymax=upper), alpha=0.5, fill='lightblue') +
  ylim(c(0,0.5)) + theme_bw() + theme(legend.position = "none") + facet_wrap(.~treatment) +
  xlab("Visit number") + ylab("Frequency of infection")


namfig<-"toenail_globalVPC.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(plot2)
dev.off()


##################### 
## Binomial model in saemix for the toenail data - complete data 
## Fit only to the subjects with complete observations over the 7 visits :

# Some graphs
barplot(table(toenail.saemix$y,toenail.saemix$visit),beside=T)
barplot(table(toenail.saemix$y,toenail.saemix$visit),beside=F)

# Graph of raw data
nobs<-tapply(toenail$visit,toenail$patientID,length)
isuj.nomiss<-names(nobs)[nobs==7]
toe2<-toe[toenail$patientID %in% isuj.nomiss,]

saemix.data2<-saemixData(name.data=toe2,name.group=c("patientID"),name.predictors=c("time","y"), name.response="y",
                        name.covariates=c("treatment"),name.X=c("time"))

binary.fit2<-saemix(saemix.model,saemix.data2,saemix.options)

######################################################################## ORDINAL
# Categorical data - Knee
library(catdata)
data(knee)

# Transformation to long format and dichotomisation of response
# knee <- reshape(knee, direction="long", varying=list(5:8), v.names="R",timevar="Time")
knee2<-cbind(knee[,c(1:5)], time=0)
xtim<-c(0,3,7,10)
for(icol in 6:8) {
  xprov<-cbind(knee[,c(1:4,icol)], time=xtim[icol-4])
  colnames(xprov)<-colnames(knee2)
  knee2<-rbind(knee2, xprov)
}
knee2<-knee2[order(knee2$N, knee2$time),c(1,6,5,2:4)]
colnames(knee2)[3]<-"y"
knee2$RD<-as.integer(knee2$y>2)
knee2$Age <- knee2$Age - 30
knee2$treatment<-knee2$Th-1
knee2$Age2 <- knee2$Age^2

knee.saemix<-knee2[,-c(4,7)]
colnames(knee.saemix)[1]<-"id"
barplot(table(knee.saemix$y,knee.saemix$time),beside=T)
ggplot(knee.saemix, aes(x=y))

if(FALSE) 
  write.table(knee.saemix, file.path(datDir, "knee.saemix.tab"), row.names=F, quote=F)

saemix.data<-saemixData(name.data=knee.saemix,name.group=c("id"),
                        name.predictors=c("y", "time"), name.X=c("time"),
                        name.covariates = c("Age","Sex","treatment","Age2"),
                        units=list(x="d",y="", covariates=c("yr","-","-","yr2")))

# Barplot
gtab <- knee.saemix %>%
  group_by(time, y) %>%
  summarise(n=length(y)) %>%
  mutate(y=as.factor(y))

kneebp <- ggplot(data = gtab, aes(x = time, y=n, group=y, fill=y)) + 
  geom_bar(stat="identity", position = "dodge") + theme_bw() + 
  scale_fill_brewer(palette = "Blues") + theme(legend.position = "top") +
  labs(fill = "Score") + xlab("Time (d)") + ylab("Counts")

namfig<-"knee_barplotData.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(kneebp)
dev.off()

ordinal.model<-function(psi,id,xidep) {
  y<-xidep[,1]
  time<-xidep[,2]
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  alp4<-psi[id,4]
  beta<-psi[id,5]
  
  logit1<-alp1 + beta*time
  logit2<-logit1+alp2
  logit3<-logit2+alp3
  logit4<-logit3+alp4
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  pge4<-exp(logit4)/(1+exp(logit4))
  pobs = (y==1)*pge1+(y==2)*(pge2 - pge1)+(y==3)*(pge3 - pge2)+(y==4)*(pge4 - pge3)+(y==5)*(1 - pge4)
  logpdf <- log(pobs)
  
  return(logpdf)
}

covmodel4<-covmodel3<-covmodel2<-covmodel<-matrix(data=0,ncol=5,nrow=4)
#covmodel[1:2,1]<-1
covmodel[,1]<-1
covmodel[,5]<-1
covmodel2[1,1]<-covmodel2[4,5]<-1
covmodel3[3,5]<-covmodel3[1,1]<-1
covmodel4[3,5]<-covmodel4[4,1]<-1

saemix.model<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                          psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                          transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)))

saemix.model.cov<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                              psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                              transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)),
                              covariate.model = covmodel)
saemix.model.cov2<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                               psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                               transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)),
                               covariate.model = covmodel2)
saemix.model.cov3<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                               psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                               transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)),
                               covariate.model = covmodel3)
saemix.model.cov4<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                               psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                               transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)),
                               covariate.model = covmodel4)

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, nb.chains=10)
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, nb.chains=10, fim=FALSE)

ord.fit<-saemix(saemix.model,saemix.data,saemix.options)
ord.fit.cov<-saemix(saemix.model.cov,saemix.data,saemix.options)
ord.fit.cov2<-saemix(saemix.model.cov2,saemix.data,saemix.options)
ord.fit.cov3<-saemix(saemix.model.cov3,saemix.data,saemix.options)
ord.fit.cov4<-saemix(saemix.model.cov4,saemix.data,saemix.options)
BIC(ord.fit)
BIC(ord.fit.cov)
BIC(ord.fit.cov2)
BIC(ord.fit.cov3)
BIC(ord.fit.cov4)

# why does this return NULL (rerun with same covariates)
compare.saemix(ord.fit.cov, ord.fit.cov2, ord.fit.cov3, ord.fit.cov4)

### Simulations for VPC
simulateOrdinal<-function(psi,id,xidep) {
  y<-xidep[,1]
  time<-xidep[,2]
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  alp4<-psi[id,4]
  beta<-psi[id,5]
  
  logit1<-alp1 + beta*time
  logit2<-logit1+alp2
  logit3<-logit2+alp3
  logit4<-logit3+alp4
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  pge4<-exp(logit4)/(1+exp(logit4))
  x<-runif(length(time))
  ysim<-1+as.integer(x>pge1)+as.integer(x>pge2)+as.integer(x>pge3)+as.integer(x>pge4)
  return(ysim)
}

nsim<-100
yfit<-ord.fit
yfit<-simulateDiscreteSaemix(yfit, simulateOrdinal, nsim=nsim)

simdat <-yfit@sim.data@datasim
simdat$time<-rep(yfit@data@data$time,nsim)
simdat$treatment<-rep(yfit@data@data$treatment,nsim)


ytab<-NULL
for(irep in 1:nsim) {
  xtab<-simdat[simdat$irep==irep,]
  xtab1 <- xtab %>%
    group_by(time, treatment, ysim) %>%
    summarise(n=length(ysim))
  ytab<-rbind(ytab,xtab1[,c("time","ysim","n","treatment")])
}
gtab <- ytab %>%
  group_by(time, treatment, ysim) %>%
  summarise(lower=quantile(n, c(0.05)), n=quantile(n, c(0.5)), upper=quantile(n, c(0.95))) %>%
  mutate(y=as.factor(ysim))

knee2 <- knee.saemix %>%
  group_by(time, treatment, y) %>%
  summarise(n=length(y)) %>%
  mutate(y=as.factor(y))


kneevpc <- ggplot(data = knee2, aes(x = time, y=n, fill=y, group=treatment)) + 
  geom_ribbon(data=gtab, aes(x=time, ymin=lower, ymax=upper), alpha=0.9, colour="lightblue") +
  geom_col(position = "dodge", width=0.5, colour="lightblue") + theme_bw() + 
  scale_fill_brewer(palette = "Blues") + theme(legend.position = "top") +
  labs(fill = "Score") + xlab("Time (d)") + ylab("Counts") + facet_wrap(treatment~y, nrow=2)

# VPC for median score in each group
knee3 <- knee.saemix %>%
  group_by(time, treatment) %>%
  summarise(mean=mean(y))

ytab<-NULL
for(irep in 1:nsim) {
  xtab<-simdat[simdat$irep==irep,]
  xtab1 <- xtab %>%
    group_by(time, treatment) %>%
    summarise(mean=mean(ysim))
  ytab<-rbind(ytab,xtab1[,c("time","treatment","mean")])
}
gtab <- ytab %>%
  group_by(time, treatment) %>%
  summarise(lower=quantile(mean, c(0.05)), mean=median(mean), upper=quantile(mean, c(0.95)))

kneeMedvpc <- ggplot(data = knee3, aes(x = time, y=mean, group=treatment)) + 
  geom_ribbon(data=gtab, aes(x=time, ymin=lower, ymax=upper), alpha=0.5, fill="lightblue") +
  geom_point(colour='blue') + theme_bw() + 
  scale_fill_brewer(palette = "Blues") + theme(legend.position = "top") +
  labs(fill = "Score") + xlab("Time (d)") + ylab("Median value of score over time") + facet_wrap(.~treatment)

namfig<-"knee_medianScoreVPC.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(kneeMedvpc)
dev.off()

#########################
#  ysim[x==0]<-1 
#  ysim[x==1]<-5
# Debug
nsim<-100
object<-ord.fit.cov3
simulate.function<-simulateOrdinal
object<-simulate.SaemixObject(object, nsim=nsim, predictions=FALSE, uncertainty=uncertainty)
simpar<-object@sim.data@sim.psi

# Simulate observations using these parameters and the simulate.function to simulate from the same model
xidep<-object@data@data[,object@data@name.predictors]
id1<-object@data@data[,"index"]
nsuj<-object@data@N
datasim<-object@sim.data@datasim
datasim$ysim<-NA
for(irep in 1:nsim) {
  psi1<-simpar[(1+(irep-1)*nsuj):(irep*nsuj),-c(1)]
  ysim<-simulate.function(psi1, id1, xidep)
  #    if(sum(is.na(ysim))>0) cat(irep,"\n")
  datasim$ysim[datasim$irep==irep]<-ysim
}
object@sim.data@datasim<-datasim

#######################################
# Analysis of dichotomised response :-/
if(FALSE) { # install glmmML, same factors found significant
  library(glmmML)
  kneeGHQ <- glmmML(RD ~ as.factor(Th) + as.factor(Sex) + Age + Age2, data=knee2,
                    family=binomial(), method="ghq", n.points=20, cluster=N)
  summary(kneeGHQ)
}

kneePQL <- glmmPQL(RD ~ as.factor(Th) + as.factor(Sex) + Age + Age2, data=knee2,
                   random = ~ 1|N, family=binomial())
summary(kneePQL)

kneePQL2 <- glmmPQL(RD ~ as.factor(Th) + Age2, data=knee2,
                    random = ~ 1|N, family=binomial())
summary(kneePQL2)

######################################################################## COUNT
# - Vraies données: epilepsy and RAPI
# - contacté David Atkins (tutorial in 2013 on analysing count data with GLMM and GEE): dataset on gender differences in drinking patterns that would be great to use as an example in saemix $\Rightarrow$ accepted ! lovely :-)
# - Salamanders data from the glmmTMB package
# - fit successful when using only the data for one species
# - but error when using more than one species with a recurrent error message (solve.default...) **TODO** investigate
# - note: error in the previous version of Poisson model (factorial(y) instead of log(factorial(y))) ?
#   

# Count data - Epilepsy
library(MASS)
data(epil)

saemix.data<-saemixData(name.data=epil, name.group=c("subject"),
                        name.predictors=c("period","y"),name.response=c("y"),
                        name.covariates=c("trt","base", "age"),
                        units=list(x="day (scaled)",y="",covariates=c("","","yr")))

countmodel.poisson<-function(psi,id,xidep) { 
  y<-xidep[,2]
  lambda<-psi[id,1]
  logp <- -lambda + y*log(lambda) - log(factorial(y))
  return(logp)
}

saemix.model<-saemixModel(model=countmodel.poisson,description="count model Poisson",modeltype="likelihood",   
                          psi0=matrix(c(0.5),ncol=1,byrow=TRUE,dimnames=list(NULL, c("lambda"))), 
                          transform.par=c(1)) #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE))

saemix.model.cov<-saemixModel(model=countmodel.poisson,description="count model Poisson",modeltype="likelihood",   
                              psi0=matrix(c(0.5),ncol=1,byrow=TRUE,dimnames=list(NULL, c("lambda"))), 
                              covariate.model = matrix(c(1,rep(1,2)), ncol=1, byrow=T),
                              transform.par=c(1)) #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE))

saemix.model.cov2<-saemixModel(model=countmodel.poisson,description="count model Poisson",modeltype="likelihood",   
                              psi0=matrix(c(0.5),ncol=1,byrow=TRUE,dimnames=list(NULL, c("lambda"))), 
                              covariate.model = matrix(c(0,1,0), ncol=1, byrow=T),
                              transform.par=c(1)) #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
poisson.fit<-saemix(saemix.model,saemix.data,saemix.options)
poisson.fit.cov<-saemix(saemix.model.cov,saemix.data,saemix.options)
poisson.fit.cov2<-saemix(saemix.model.cov2,saemix.data,saemix.options)
plot(poisson.fit, plot.type="convergence")

## ZIP Poisson model
countmodel.zip<-function(psi,id,xidep) {
  y<-xidep[,2]
  lambda<-psi[id,1]
  p0<-psi[id,2]
  logp <- log(1-p0) -lambda + y*log(lambda) - log(factorial(y))
  logp0 <- log(p0+(1-p0)*exp(-lambda))
  logp[y==0]<-logp0[y==0]
  return(logp)
}

saemix.model.zip<-saemixModel(model=countmodel.zip,description="count model ZIP",modeltype="likelihood",   
                              psi0=matrix(c(0.5,0.2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","p0"))), 
                              transform.par=c(1,3), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                              covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE))
saemix.model.zipcov<-saemixModel(model=countmodel.zip,description="count model ZIP",modeltype="likelihood",   
                                 psi0=matrix(c(0.5,0.2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","p0"))), 
                                 transform.par=c(1,3), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                                 covariate.model = matrix(c(1,rep(1,5)), ncol=2, byrow=T),
                                 covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE))
zippoisson.fit<-saemix(saemix.model.zip,saemix.data,saemix.options)
zippoisson.fitcov<-saemix(saemix.model.zipcov,saemix.data,saemix.options)
plot(zippoisson.fit, plot.type="convergence")

BIC(poisson.fit)
BIC(zippoisson.fit)

##Generalized Poisson model
countmodel.genpoisson<-function(psi,id,xidep) {
  y<-xidep[,1]
  delta<-psi[id,1]
  lambda<-psi[id,2]
  logp <- -lambda
  pos.ind <- which(y>0)
  logp[pos.ind] <- log(lambda) + (y-1)*log(lambda+y*delta) - (lambda+y*delta) - log(factorial(y))
  return(logp)
}

saemix.model.gp<-saemixModel(model=countmodel.zip,description="Generalised Poisson model",modeltype="likelihood",   
                              psi0=matrix(c(0.5,0.2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("delta","lambda"))), 
                              transform.par=c(1,1), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                              covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE))
genpoisson.fit<-saemix(saemix.model.gp,saemix.data,saemix.options)

########################################################################
# Count data - RAPI
atkinsDir<-"/home/eco/work/saemix/examples/coutRegTutorialAtkins"
rapi.df <- read.csv(file.path(atkinsDir,"RAPI.Final.csv"), header = TRUE)

### Set categorical variables as factors
rapi.df <- within(rapi.df, {
  gender <- factor(gender, 0:1, c("Women","Men"))
})

### Basic file summaries
summary(rapi.df)

### Number of participants
length(unique(rapi.df$id)) # 818

### How much data per person?
table(table(rapi.df$id)) # 561 have all 5 assessments

rapi.saemix<-rapi.df[,c(1,4,2,3)]

if(FALSE) 
  write.table(rapi.saemix, file.path(datDir, "rapi.saemix.tab"), quote=FALSE, row.names=FALSE)

# Data
saemix.data<-saemixData(name.data=rapi.saemix, name.group=c("id"),
                        name.predictors=c("time","rapi"),name.response=c("rapi"),
                        name.covariates=c("gender"),
                        units=list(x="months",y="",covariates=c("")))
hist(rapi.saemix$rapi, main="", xlab="RAPI score", breaks=30)

# Models
# Poisson with a time effect
count.poisson<-function(psi,id,xidep) { 
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  logp <- -lambda + y*log(lambda) - log(factorial(y))
  return(logp)
}

## ZIP Poisson model with time effect
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

## Generalized Poisson model with time effect
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


# C^n_p = n! /(p! (n-p)!)

## Negative binomial model with time effect
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

# Fits
## Poisson
### Model without covariate
saemix.model.poi<-saemixModel(model=count.poisson,description="Count model Poisson",modeltype="likelihood",   
                              psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                              transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)))

### Gender effect on intercept and slope
saemix.model.poi.cov2<-saemixModel(model=count.poisson,description="Count model Poisson",modeltype="likelihood",   
                                   psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                   transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)),
                                   covariance.model =matrix(data=1, ncol=2, nrow=2),
                                   covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE))
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)

### Fit with saemix
poisson.fit<-saemix(saemix.model.poi,saemix.data,saemix.options)
poisson.fit.cov2<-saemix(saemix.model.poi.cov2,saemix.data,saemix.options)
exp(poisson.fit@results@fixed.effects)
exp(poisson.fit.cov2@results@fixed.effects)

### Simulations
saemix.simulatePoisson<-function(psi, id, xidep) {
  time<-xidep[,1]
  y<-xidep[,2]
  intercept<-psi[id,1]
  slope<-psi[id,2]
  lambda<- exp(intercept + slope*time)
  y<-rpois(length(time), lambda=lambda)
  return(y)
}

yfit1<-simulateDiscreteSaemix(poisson.fit.cov2, 10, saemix.simulatePoisson)

hist(yfit1@data@data$rapi, xlim=c(0,50), freq=F, breaks=30, xlab="Observed counts", main="")
lines(density(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim<50]), lwd = 2, col = 'red')

cat("Observed proportion of 0's", length(yfit1@data@data$rapi[yfit1@data@data$rapi==0])/yfit1@data@ntot.obs,"\n")
cat("      Poisson model, p=",length(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim==0])/length(yfit1@sim.data@datasim$ysim),"\n")

## ZIP
### base model
saemix.model.zip<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                              psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                              transform.par=c(0,0,3), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.3,0)))

### ZIP Poisson with gender on both intercept
saemix.model.zip.cov1<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                                   psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                                   transform.par=c(0,0,3), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.3,0)),
                                   covariate.model = matrix(c(1,0,0),ncol=3, byrow=TRUE))
### ZIP Poisson with gender on both intercept and slope
saemix.model.zip.cov2<-saemixModel(model=count.poissonzip,description="count model ZIP",modeltype="likelihood",   
                                   psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","p0"))), 
                                   transform.par=c(0,0,3), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.3,0)),
                                   covariate.model = matrix(c(1,1,0),ncol=3, byrow=TRUE))

zippoisson.fit<-saemix(saemix.model.zip,saemix.data,saemix.options)
zippoisson.fit.cov1<-saemix(saemix.model.zip.cov1,saemix.data,saemix.options)
zippoisson.fit.cov2<-saemix(saemix.model.zip.cov2,saemix.data,saemix.options)

### Simulations

yfit2<-simulateDiscreteSaemix(zippoisson.fit.cov2, 10, saemix.simulatePoissonZIP)
if(FALSE) {
  par(mfrow=c(1,2))
  hist(yfit@sim.data@datasim$ysim[yfit@sim.data@datasim$ysim<50], xlim=c(0,50), freq=F, breaks=20, xlab="Simulated counts", main="By hand")
  hist(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim<50], xlim=c(0,50), freq=F, breaks=20, xlab="Simulated counts", main="Function")
}

par(mfrow=c(1,3))
hist(yfit1@data@data$rapi, xlim=c(0,50), freq=F, breaks=30, xlab="Observed counts", main="")
hist(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim<50], xlim=c(0,50), freq=F, breaks=20, xlab="Simulated counts", main="Poisson model")
hist(yfit2@sim.data@datasim$ysim[yfit2@sim.data@datasim$ysim<50], xlim=c(0,50), freq=F, breaks=20, xlab="Simulated counts", main="ZIP model")

# As densities - very weird, must be doing something wrong with the plot :-/
par(mfrow=c(1,1))
hist(yfit1@data@data$rapi, xlim=c(0,50), freq=F, breaks=30, xlab="Observed counts", main="")
lines(density(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim<50]), lwd = 2, col = 'red')
lines(density(yfit2@sim.data@datasim$ysim[yfit2@sim.data@datasim$ysim<50]), lwd = 2, col = 'blue')

cat("      ZIP model,     p=",length(yfit2@sim.data@datasim$ysim[yfit2@sim.data@datasim$ysim==0])/length(yfit2@sim.data@datasim$ysim),"\n")


## Generalised Poisson - not sure the results make much sense
saemix.model.genpoi<-saemixModel(model=count.genpoisson,description="Generalised Poisson model",modeltype="likelihood",   
                                 psi0=matrix(c(1.5, 0.01, 0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","delta"))), 
                                 transform.par=c(0,0,1), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.5,0.5)))

genpoisson.fit<-saemix(saemix.model.genpoi,saemix.data,saemix.options)

## Negative binomial
saemix.model.NB<-saemixModel(model=count.NB,description="Negative binomial model",modeltype="likelihood",   
                             psi0=matrix(c(1.5, 0.01, 2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("intercept", "slope","k"))), 
                             transform.par=c(0,0,1), covariance.model=diag(c(1,1,0)), omega.init=diag(c(0.5,0.5,0.5)))

negbin.fit<-saemix(saemix.model.NB,saemix.data,saemix.options)

# Comparing parameters across fits - guessing intercept and slope shouldn't change much
poisson.fit@results@fixed.effects
zippoisson.fit@results@fixed.effects
zippoisson.fit.cov2@results@fixed.effects
genpoisson.fit@results@fixed.effects
negbin.fit@results@fixed.effects

## Hurdle - 2 models 
saemix.data1<-saemixData(name.data=rapi.saemix[rapi.saemix$rapi>0,], name.group=c("id"),
                         name.predictors=c("time","rapi"),name.response=c("rapi"),
                         name.covariates=c("gender"),
                         units=list(x="week",y="",covariates=c("")))

rapi.saemix$y0<-as.integer(rapi.saemix$rapi>0)
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
  logpdf<-rep(0,length(tim))
  P.obs = (y==0)*(1-pevent)+(y==1)*pevent
  logpdf <- log(P.obs)
  return(logpdf)
}

saemix.hurdle0<-saemixModel(model=binary.model,description="Binary model",
                            modeltype="likelihood",
                            psi0=matrix(c(-1.5,-.1,0,0),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                            transform.par=c(0,0), covariate.model=c(1,1),
                            covariance.model=matrix(c(1,0,0,1),ncol=2), omega.init=diag(c(1,0.3)))

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE, displayProgress=FALSE)

hurdlefit0<-saemix(saemix.hurdle0,saemix.data0,saemix.options)

# proportion of 0's in the data
rapi.tab <- table(rapi.saemix$rapi == 0)
rapi.tab/sum(rapi.tab)

# Fit Poisson model to saemix.data1
saemix.hurdle1.cov2<-saemixModel(model=count.poisson,description="Count model Poisson",modeltype="likelihood",   
                                 psi0=matrix(c(log(5),0.01),ncol=2,byrow=TRUE,dimnames=list(NULL, c("intercept","slope"))), 
                                 transform.par=c(0,0), omega.init=diag(c(0.5, 0.5)),
                                 covariance.model =matrix(data=1, ncol=2, nrow=2),
                                 covariate.model=matrix(c(1,1), ncol=2, byrow=TRUE))
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)

hurdlefit1<-saemix(saemix.hurdle1.cov2,saemix.data1,saemix.options)

summary(hurdlefit0)
summary(hurdlefit1)

# Table form - compare to column B in Table 2
yfit0<-hurdlefit0
yfit1<-hurdlefit1

rr.tab<-data.frame(param=c("intercept", "beta.Male.inter", "slope", "beta.Male.slope", "omega.inter","omega.slope"), 
                   poissonNoZero=c(yfit1@results@fixed.effects, c(sqrt(diag(yfit1@results@omega)))),
                   logistic=c(yfit0@results@fixed.effects, c(sqrt(diag(yfit0@results@omega)))))

print(rr.tab)

######################################################################## COUNT other
####################################### Other
# Count data - Encephalitis
if(FALSE) {
  data(encephalitis)
  encephalitis$BAV<-encephalitis$country
  encephalitis$BAV[encephalitis$BAV==2]<-0
  
  # log-linear Poisson model, with dependence on country and year
  enc1 <- glm(count ~ year+I(year^2)+BAV+year*BAV, family = poisson, data=encephalitis)
  summary(enc1)
  
  enc2 <- glm(count ~ year+BAV, family = poisson, data=encephalitis)
  summary(enc2)
}

# homi <- read.table("http://www.stat.ufl.edu/~aa/glm/data/Homicides.dat",  header = TRUE)

## Salamander data from glmmTMB
library(glmmTMB)
data(Salamanders)
summary(Salamanders)

barplot(Salamanders$DOY, Salamanders$count)

table(Salamanders$spp, Salamanders$count)

ecl<-Salamanders[Salamanders$spp=="EC-L",]
twospecies<-Salamanders[Salamanders$spp %in% c("EC-L","DF"),]
twospecies$spp <- as.character(twospecies$spp)

saemix.data<-saemixData(name.data=ecl, name.group=c("site"),
                        name.predictors=c("DOY","count"),name.response=c("count"),
                        name.covariates=c("mined","cover", "DOP", "Wtemp"),
                        units=list(x="day (scaled)",y="",covariates=c("","","","")))

countmodel.poisson<-function(psi,id,xidep) { 
  y<-xidep[,2]
  lambda<-psi[id,1]
  logp <- -lambda + y*log(lambda) - log(factorial(y))
  return(logp)
}

saemix.model<-saemixModel(model=countmodel.poisson,description="count model Poisson",modeltype="likelihood",   
                          psi0=matrix(c(0.5),ncol=1,byrow=TRUE,dimnames=list(NULL, c("lambda"))), 
                          transform.par=c(1)) #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE))

saemix.model.cov<-saemixModel(model=countmodel.poisson,description="count model Poisson",modeltype="likelihood",   
                          psi0=matrix(c(0.5),ncol=1,byrow=TRUE,dimnames=list(NULL, c("lambda"))), 
                          covariate.model = matrix(c(1,rep(0,3)), ncol=1, byrow=T),
                          transform.par=c(1)) #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
poisson.fit<-saemix(saemix.model,saemix.data,saemix.options)
poisson.fit.cov<-saemix(saemix.model.cov,saemix.data,saemix.options)
plot(poisson.fit, plot.type="convergence")

## ZIP Poisson model
countmodel.zip<-function(psi,id,xidep) {
  y<-xidep[,2]
  lambda<-psi[id,1]
  p0<-psi[id,2]
  logp <- log(1-p0) -lambda + y*log(lambda) - log(factorial(y))
  logp0 <- log(p0+(1-p0)*exp(-lambda))
  logp[y==0]<-logp0[y==0]
  return(logp)
}

saemix.model.zip<-saemixModel(model=countmodel.zip,description="count model ZIP",modeltype="likelihood",   
                          psi0=matrix(c(0.5,0.2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","p0"))), 
                          transform.par=c(1,3), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                          covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE))
saemix.model.zipcov<-saemixModel(model=countmodel.zip,description="count model ZIP",modeltype="likelihood",   
                              psi0=matrix(c(0.5,0.2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","p0"))), 
                              transform.par=c(1,3), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                              covariate.model = matrix(c(1,rep(0,7)), ncol=2, byrow=T),
                              covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE))
zippoisson.fit<-saemix(saemix.model.zip,saemix.data,saemix.options)
zippoisson.fitcov<-saemix(saemix.model.zipcov,saemix.data,saemix.options)
plot(zippoisson.fit, plot.type="convergence")

BIC(poisson.fit)
BIC(zippoisson.fit)

if(FALSE) {
  salamander.saemix<-ecl
  salamander.saemix$mined<-as.integer(salamander.saemix$mined=="yes")
  write.table(salamander.saemix, file.path(saemixDir, "testbelhal", "salamander.saemix.tab"), quote=F, row.names=F)
  
  salamander.saemix<-Salamanders
  species<-as.character(unique(Salamanders$spp))
  for(isp in species) {
    xcov<-as.integer(salamander.saemix$spp==isp)
    salamander.saemix<-cbind(salamander.saemix, sp=xcov)
    colnames(salamander.saemix)[dim(salamander.saemix)[2]]<-isp
  }
  salamander.saemix$mined<-as.integer(salamander.saemix$mined=="yes")
  salamander.saemix$spp<-NULL
  salamander.saemix<-salamander.saemix[,c("site","DOY","count","mined","cover","DOP","Wtemp","sample",colnames(salamander.saemix)[9:15])]
  write.table(salamander.saemix, file.path(saemixDir, "testbelhal", "sal.saemix.tab"), quote=F, row.names=F)
  
  saemix.data<-saemixData(name.data=salamander.saemix, name.group=c("site"),
                          name.predictors=c("DOY","count"),name.response=c("count"),
                          name.covariates=c("mined","PR","DM","EC-A","EC-L","DES-L","DF"), # reference class will be GP
                          units=list(x="day (scaled)",y=""))
  saemix.model.zipcov<-saemixModel(model=countmodel.zip,description="count model ZIP",modeltype="likelihood",   
                                   psi0=matrix(c(0.5,0.2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","p0"))), 
                                   transform.par=c(1,3), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                                   covariate.model = matrix(c(rep(1,2*length(saemix.data@name.covariates))), ncol=2, byrow=T),
                                   covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE))
  zippoisson.fitcov<-try(saemix(saemix.model.zipcov,saemix.data,saemix.options)) # fails with Lapack problem => really need to debug this...
  
}

# Same, with the 2 most frequent species, using spp as a covariate - fails

saemix.data<-saemixData(name.data=twospecies, name.group=c("site"),
                        name.predictors=c("DOY","count"),name.response=c("count"),
                        name.covariates=c("spp","mined","cover", "DOP", "Wtemp"),
                        units=list(x="day (scaled)",y="",covariates=c("","","","","")))

countmodel.poisson<-function(psi,id,xidep) { 
  y<-xidep[,2]
  lambda<-psi[id,1]
  logp <- -lambda + y*log(lambda) - log(factorial(y))
  return(logp)
}

saemix.model<-saemixModel(model=countmodel.poisson,description="count model Poisson",modeltype="likelihood",   
                          psi0=matrix(c(1),ncol=1,byrow=TRUE,dimnames=list(NULL, c("lambda"))), 
                          covariate.model = matrix(c(1,rep(0,4)),ncol=1, byrow=T),
                          transform.par=c(1)) #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, fim=FALSE)
poisson.fit<-try(saemix(saemix.model,saemix.data,saemix.options)) # fails, pb with solving Lapack

saemix.model.zip<-saemixModel(model=countmodel.zip,description="count model ZIP",modeltype="likelihood",   
                              psi0=matrix(c(0.5,0.2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","p0"))), 
                              transform.par=c(1,3), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                              covariate.model = matrix(c(1,rep(0,9)),ncol=2, byrow=T),
                              covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE))
zippoisson.fit<-try(saemix(saemix.model.zip,saemix.data,saemix.options))

## Fit with glmmTMB
zipm1 = glmmTMB(count~spp + mined + (1|site), Salamanders, family="poisson")
summary(zipm1)
zipm2 = glmmTMB(count~spp * mined + (1|site), Salamanders, family="poisson")
summary(zipm2)

zipm3 = glmmTMB(count~spp * mined + (1|site), zi=~spp * mined, Salamanders, family="poisson")
summary(zipm3)

### 

## nest data from glmmTMB

library(glmmTMB)
data(Owls)
Owls$time<-Owls$ArrivalTime-min(Owls$ArrivalTime)
Owls<-Owls[order(Owls$Nest, Owls$ArrivalTime),]
Owls[!duplicated(Owls$Nest),]
colnames(Owls)[5]<-"NbNego"

owls.saemix<-Owls

?# saemix
saemix.data<-saemixData(name.data=owls.saemix,header=TRUE,name.group=c("Nest"),
                        name.predictors=c("time","NbNego"),name.response=c("NbNego"),
                        name.covariates=c("FoodTreatment", "SexParent", "NegPerChick", "BroodSize"),
                        units=list(x="min",y="",covariates=c("","","","")))

#Basic Poisson model
countmodel.poisson<-function(psi,id,xidep) { 
  y<-xidep[,1]
  lambda<-psi[id,1]
  dummy<-psi[id,2]
  logp <- -lambda + y*log(lambda) - factorial(y)
  return(logp)
}

saemix.model<-saemixModel(model=countmodel.poisson,description="count model",modeltype="likelihood",   
                          psi0=matrix(c(0.5,1),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","dummy"))), 
                          transform.par=c(1,1), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                          covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE),
                          fixed.estim=c(1,0))

##Generalized Poisson model
countmodel.genpoisson<-function(psi,id,xidep) {
  y<-xidep[,1]
  delta<-psi[id,1]
  lambda<-psi[id,2]
  logp <- -lambda
  pos.ind <- which(y>0)
  logp[pos.ind] <- log(lambda) + (y-1)*log(lambda+y*delta) - (lambda+y*delta) - factorial(y)
  return(logp)
}

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
poisson.fit<-saemix(saemix.model,saemix.data,saemix.options)
plot(poisson.fit, plot.type="convergence")

# Grouseticks Data on red grouse ticks from Elston et al. 2001 (Number of ticks on the heads of red grouse chicks sampled in the field (grouseticks) and an aggregated version (grouseticks_agg); see original source for more details)

library(lme4)
data("grouseticks")
form <- TICKS~YEAR+HEIGHT+(1|BROOD)+(1|INDEX)+(1|LOCATION)
full_mod1  <- glmer(form, family="poisson",data=grouseticks)



####### Simulated count data

# Settings  
param <- c(39.1, 0.0388, 0.1 )
omega<-c(0.5, 0.5) # SD=50%
paramSimul<-c(param, omega)
parnam<-c("alpha","beta","risk","omega.alpha","omega.beta")

nsuj<-40
xtim<-c(0.0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)

partab<-as.data.frame(matrix(data=0,nrow=nsuj,ncol=2,dimnames=list(NULL,parnam[1:2])))
for(i in 1:2) partab[,i]<-rnorm(nsuj,mean=log(param[i]),sd=omega[i])
partab[(1+nsuj/2):nsuj,2]<-partab[(1+nsuj/2):nsuj,2]+param[3]
for(i in 1:2) partab[,i]<-exp(partab[,i])

psim<-data.frame()
for(itim in xtim) {
  lambda<-partab[,1]*exp(-partab[,2]*itim)
  psim<-rbind(psim,lambda)
}
datsim<-data.frame(id=rep(1:nsuj,each=length(xtim)),time=rep(xtim,nsuj),lambda=unlist(psim))
rownames(datsim)<-NULL
ysim<-rpois(dim(datsim)[1], lambda=datsim$lambda)
#  summary(datsim)
datsim$y<-ysim
datsim$risk<-ifelse(datsim$id>(nsuj/2),1,0)


saemix.data<-saemixData(name.data=datsim,name.group=c("id"),name.predictors=c("time","y"), name.covariates=c("risk"),name.X=c("time"))

# Model
countData.model<-function(psi,id,xidep) {
  tim <- xidep[,1]
  y <- xidep[,2] 
  alpha <- psi[id,1]
  beta <- psi[id,2]
  lambda <- alpha*exp(-beta*tim)
  
  logpdf <- rep(0,length(tim))
  logpdf <- -lambda + y*( (log(alpha) - beta*tim )) - log(factorial(y))
  return(logpdf) 
}
saemix.model.true<-saemixModel(model=countData.model,description="Count data model", modeltype="likelihood",
                               psi0=matrix(c(param[1:2],0,param[3]),ncol=2,byrow=TRUE,dimnames=list(NULL,parnam[1:2])),
                               covariate.model=matrix(c(0,1),ncol=2), omega.init = diag(c(0.5,0.5)),
                               transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2))
# Running saemix

saemix.options<-list(seed=123456,save=FALSE,save.graphs=FALSE, fim=FALSE, displayProgress=FALSE)
count.fit<-try(saemix(saemix.model.true,saemix.data,saemix.options))

########################################################################
# TTE - Lung cancer

## Creating the dataset for lung cancer
if(FALSE) {
  library(survival)
  data(cancer)
  cancer$cens<-as.integer(cancer$status==1) # censored=1, non-censored=0
  cancer$status<-cancer$status-1 # dead=1, alive=0
  cancer<-cbind(id=1:dim(cancer)[1],cancer)
  cancer2<-cancer
  cancer2$time<-0
  cancer2$status<-0
  cancer2$cens<-0
  lung.saemix<-rbind(cancer2, cancer)
  lung.saemix<-lung.saemix[order(lung.saemix$id, lung.saemix$time),]
  lung.saemix$sex<-lung.saemix$sex-1
  lung.saemix<-lung.saemix[,c("id","time","status","cens","inst","age", 
                              "sex", "ph.ecog", "ph.karno", "pat.karno", "wt.loss","meal.cal")]
  hasnoNA<-function(xmat) 
    apply(xmat,1,function(x) sum(is.na(x))==0)
  lung.saemix<-lung.saemix[hasnoNA(lung.saemix[,5:9]),]
  write.table(lung.saemix, file.path(datDir, "lung.saemix.tab"), quote=F, row.names=F)
}

lung.saemix

## 

weibulltte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  lambda <- psi[id,1] # Parameters of the Weibull model
  beta <- psi[id,2]
  Nj <- length(T)
  
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  hazard <- (beta/lambda)*(T/lambda)^(beta-1) # H'
  H <- (T/lambda)^beta # H
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))
  return(logpdf)
}

######################################################################## TTE - OTHER
# TTE - PBC
data("pbc2.id", package = "JM")

pbc<-pbc2.id
pbc$cens<-as.integer(pbc$status!="dead") # censored=1, non-censored=0
pbc$status<-as.integer(pbc$status=="dead")  # dead=1, alive=0
pbc2<-pbc
pbc2$years<-0
pbc2$status<-pbc2$cens<-0
pbc2<-rbind(pbc2, pbc)
colnames(pbc2)[colnames(pbc2)=="years"]<-"time"
pbc.saemix<-pbc2[order(pbc2$id, pbc2$time),]

saemix.data<-saemixData(name.data=pbc.saemix,header=TRUE,name.group=c("id"),
                        name.predictors=c("time","status","cens"),name.response=c("status"),
                        name.covariates=c("age", "sex", "drug", "ascites", "hepatomegaly", "spiders","alkaline"),
                        units=list(x="days",y="",covariates=c("yr","","-","","","","")))

exptte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  beta <- psi[id,1]
  Nj <- length(T)
  
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  hazard <- 1/beta # H'
  H <- T/beta # H
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))
  return(logpdf)
}

saemix.Weibmodel<-saemixModel(model=weibulltte.model,description="time model",modeltype="likelihood",
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","beta"))),
                          transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE))
saemix.Expmodel<-saemixModel(model=exptte.model,description="time model",modeltype="likelihood",
                              psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL, c("beta","dummy"))),
                              transform.par=c(1,0), fixed.estim = c(1,0),
                             covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE))

saemix.WeibCov<-saemixModel(model=weibulltte.model,description="time model",modeltype="likelihood",
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("lambda","beta"))),
                          covariate.model=matrix(c(1,0,1,0,1,0,rep(0,8)),ncol=2, byrow=TRUE),
                          transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE))
saemix.ExpCov<-saemixModel(model=exptte.model,description="time model",modeltype="likelihood",
                           psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL, c("beta","dummy"))),
                           transform.par=c(1,0), fixed.estim = c(1,0),
                           covariate.model=matrix(c(1,0,1,0,1,0,rep(0,8)),ncol=2, byrow=TRUE),
                           covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
tte.Weibfit<-saemix(saemix.Weibmodel,saemix.data,saemix.options)
tte.Weibcov<-saemix(saemix.WeibCov,saemix.data,saemix.options)
tte.Expfit<-saemix(saemix.Expmodel,saemix.data,saemix.options)
tte.Expcov<-saemix(saemix.ExpCov,saemix.data,saemix.options)
plot(tte.fit, plot.type="convergence")

ypred<-predict(tte.fit)

if(FALSE) 
  write.table(pbc.saemix, file.path(datDir, "pbc.saemix.tab"), quote=F, row.names=F)

## Call:
## survreg(formula = Surv(years, status2) ~ drug + sex + age, data = pbc2.id)
##                  Value Std. Error     z      p
## (Intercept)    4.13862    0.48260  8.58 <2e-16
## drugD-penicil  0.13046    0.15576  0.84  0.402
## sexfemale      0.44472    0.20147  2.21  0.027
## age           -0.03874    0.00793 -4.89  1e-06
## Log(scale)    -0.10223    0.07423 -1.38  0.168

# Exponential, JM
##                 Value Std. Error     z       p
## (Intercept)    4.3282     0.5159  8.39 < 2e-16
## drugD-penicil  0.1438     0.1721  0.84    0.40
## sexfemale      0.4816     0.2217  2.17    0.03
## age           -0.0420     0.0084 -5.00 5.7e-07

########################################################################
# RTTE - Gaucher ? (few events...)

########################################################################
# RTTE - Simulated data

# For repeated time-to-events, a simulation algorithm is to simulate TTE repeatedly in an individual until the time to censoring has been reached (Penichou et al. 2014). Here the simulation for each successive event is perfomed by simulating a random variable in $\mathcal{U}[0,1]$ and using the inverse of the cumulative density function to generate the corresponding time to event.

simul.rtte.unif<-function(psi) { # xidep, id not important, we only use psi
  censoringtime <- 3
  maxevents <- 30
  lambda <- psi[,1]
  beta <- psi[,2]
  simdat<-NULL
  N<-nrow(psi)
  for(i in 1:N) {
    eventTimes<-c(0)
    T<-0
    Vj<-runif(1)
    #    T <- (-log(Vj)*lambda[i])^(beta[i])
    T<-lambda[i]*(-log(Vj))^(1/beta[i])
    nev<-0
    while (T < censoringtime & nev<maxevents){
      eventTimes <- c(eventTimes, T)  
      nev<-nev+1
      Vj<-runif(1)
      #      T <- T+(-log(Vj)*lambda[i])^(beta[i])
      #      T<-(-log(Vj)*lambda[i] + T^(1/beta[i]))^(beta[i])
      T<-lambda[i]*(-log(Vj) + (T/lambda[i])^(beta[i]))^(1/beta[i])
    }
    if(nev==maxevents) {
      message("Reached maximum number of events\n")
    }
    eventTimes<-c(eventTimes, censoringtime)
    cens<-rep(1,length(eventTimes))
    cens[1]<-cens[length(cens)]<-0
    simdat<-rbind(simdat,
                  data.frame(id=i, T=eventTimes, status=cens))
  }
  return(simdat)
}

# Subjects
set.seed(12345)
param<-c(2, 1.5, 0.5)
# param<-c(4, 1.2, 0.3)
omega<-c(0.25,0.25)
nsuj<-200
risk<-rep(0,nsuj)
risk[(nsuj/2+1):nsuj]<-1
psiM<-data.frame(lambda=param[1]*exp(rnorm(nsuj,sd=omega[1])), beta=param[2]*exp(param[3]*risk+rnorm(nsuj,sd=omega[2])))
simdat <- simul.rtte.unif(psiM)
simdat$risk<-as.integer(simdat$id>(nsuj/2))

if(FALSE) { # Check simulated parameters
  summary(psiM)
  apply(log(psiM),2,sd)
}

if(FALSE) {
  xtim<-seq(0:20)
  rtte.data<-data.frame(id=rep(1:nsuj,each=length(xtim)),time=rep(xtim,nsuj))
  preds <- simul.rtte.rweib(psiM, rtte.data$id, rtte.data[,c("time")])
  
  par(mfrow=c(1,2))
  hist(preds)
  hist(simdat$T[simdat$T>0])
  
  rtte.data$tlat<-preds
  rtte.data$status<-as.integer(rtte.data$tlat<3)
  
  # Remove duplicated censored times
  dat1<-NULL
  for(i in 1:nsuj) {
    idat<-rtte.data[rtte.data$id==i,]
    idat<-idat[!duplicated(idat$tlat),,drop=FALSE]
    dat1<-rbind(dat1, c(i,0,0), idat[,-c(2)])
  }
  rtte.data<-dat1
  table(tapply(rtte.data$id, rtte.data$id, length))
}

if(FALSE)
  write.table(simdat,file.path(ecoDir,"simulatedRTTE.csv"), quote=F, row.names=F)

saemix.data<-saemixData(name.data=simdat, name.group=c("id"), name.predictors=c("T"), name.response="status", name.covariates="risk")

### Fit RTTE
rtte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  N <- nrow(psi) # nb of subjects
  Nj <- length(T) # nb of events (including 0 and censoring times)
  # censoringtime = 6
  censoringtime = max(T) # same censoring for everyone
  lambda <- psi[id,1]
  beta <- psi[id,2]
  tinit <- which(T==0) # indices of beginning of observation period
  tcens <- which(T==censoringtime) # indices of censored events 
  tevent <- setdiff(1:Nj, append(tinit,tcens)) # indices of non-censored event times
  hazard <- (beta/lambda)*(T/lambda)^(beta-1)
  H <- (T/lambda)^beta
  logpdf <- rep(0,Nj)
  logpdf[tcens] <- -H[tcens] + H[tcens-1]
  logpdf[tevent] <- -H[tevent] + H[tevent-1] + log(hazard[tevent])
  return(logpdf)
}

saemix.model.base<-saemixModel(model=rtte.model,description="Repeated TTE model",modeltype="likelihood",
                               psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("lambda","beta"))),
                               transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, byrow=TRUE))
saemix.model<-saemixModel(model=rtte.model,description="Repeated TTE model",modeltype="likelihood",
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("lambda","beta"))),
                          transform.par=c(1,1),covariate.model=matrix(c(0,1),ncol=2),
                          covariance.model=matrix(c(1,0,0,1),ncol=2, byrow=TRUE))
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, displayProgress=FALSE)
rtte.fit<-saemix(saemix.model,saemix.data,saemix.options)
plot(rtte.fit, plot.type="convergence")
