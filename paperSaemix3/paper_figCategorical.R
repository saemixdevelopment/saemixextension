# Loading libraries
library(xtable)
library(ggplot2)
# library(tidyr)

# Loading saemix
library(saemix)

# Folders
# workDir<-"/home/eco/work/saemix/saemixextension/paperSaemix3"
workDir <- getwd() # working directory should be paperSaemix3
# setwd(workDir)

figDir <- file.path(workDir, "figs")
saveFigs <- TRUE

# Font size correction for saving to file
rsize.text <- 2
rsize.ticks <- 1.5
theme_set(theme_bw(base_size = 15)) 

###################################################### Data exploration
# Data
data(knee.saemix)

# Data
ordknee.data<-saemixData(name.data=knee.saemix,name.group=c("id"),
                         name.predictors=c("y", "time"), name.X=c("time"),
                         name.covariates = c("Age","Sex","treatment","Age2"),
                         units=list(x="d",y="", covariates=c("yr","-","-","yr2")), verbose=FALSE)

plotDiscreteData(ordknee.data, outcome="categorical", which.cov="treatment")

plotDiscreteData(ordknee.data, outcome="categorical", which.cov="Sex")

if(saveFigs) {
  plot1 <- plotDiscreteData(ordknee.data, outcome="categorical", which.cov="treatment")
#  plot1 +  theme(text = element_text(size=rel(rsize.text)), legend.text=element_text(size=rel(rsize.text)), strip.text.x = element_text(size=rel(rsize.ticks)), strip.text.y = element_text(size=rel(rsize.ticks)))
  namfig<-"knee_rawDataPropTime.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
#  plot(plot2)
  plot(plot1)
  dev.off()
}

###################################################### Fit proportional odds model

# Model for ordinal responses
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
# simulate function
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

# Fitting
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, nb.chains=10, nbiter.saemix=c(600,100), print=FALSE)
#saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, nb.chains=10, fim=FALSE)

# Saemix model without covariates
saemix.model<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                            simulate.function=simulateOrdinal, psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5, byrow=TRUE, 
                                                                           dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), transform.par=c(0,1,1,1,1),
                            omega.init=diag(c(100, 1, 1, 1, 1)), covariance.model = diag(c(1,0,0,0,1)), verbose=FALSE)
ord.fit<-saemix(saemix.model,ordknee.data,saemix.options)
summary(ord.fit)

###################################################### Covariate models
# Model with covariates resulting from the stepwise algorithm
## IIV: all alphas, none on beta :-/
## Covariates:   alp1(Age2)alp2(treatment)beta(treatment)
covariate.model <- matrix(data=0, nrow=4, ncol=5)
covariate.model[3,2]<-covariate.model[3,5]<-covariate.model[4,1]<-1
ordmodel.cov<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                          simulate.function=simulateOrdinal, psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5, byrow=TRUE, 
                                                                         dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), transform.par=c(0,1,1,1,1),
                          omega.init=diag(c(100, 1, 1, 1, 1)), covariate.model=covariate.model, covariance.model = diag(c(1,1,1,1,0)), verbose=FALSE)

# Fitting
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, nb.chains=10, nbiter.saemix=c(600,100), print=FALSE)
#saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, nb.chains=10, fim=FALSE)

ord.fit.cov<-saemix(ordmodel.cov,ordknee.data,saemix.options)
summary(ord.fit.cov)

###################################################### VPC
# VPC, stratified by treatment and by sex
nsim<-100
yfit<-ord.fit.cov
yfit<-simulateDiscreteSaemix(yfit, nsim=nsim)
discreteVPC(yfit, outcome="categorical")
discreteVPC(yfit, outcome='categorical',covsplit=TRUE, which.cov="treatment")
discreteVPC(yfit, outcome='categorical',covsplit=TRUE, which.cov="Sex")

if(saveFigs) {
  plot1 <- discreteVPC(yfit, outcome='categorical',covsplit=TRUE, which.cov="treatment")
  plot2 <- plot1 # + theme(text = element_text(size=rel(rsize.text)), legend.text=element_text(size=rel(rsize.text)), strip.text.x = element_text(size=rel(rsize.ticks)), strip.text.y = element_text(size=rel(rsize.ticks)))
  namfig<-"knee_VPCbytreatment.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot2)
  dev.off()
  plot1 <-  discreteVPC(yfit, outcome='categorical',covsplit=TRUE, which.cov="Sex")
   plot2 <- plot1 # + theme(text = element_text(size=rel(rsize.text)), legend.text=element_text(size=rel(rsize.text)), strip.text.x = element_text(size=rel(rsize.ticks)), strip.text.y = element_text(size=rel(rsize.ticks)))
 namfig<-"knee_VPCbySex.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot2)
  dev.off()
}

# VPC on median score, stratified by treatment

knee3<-NULL
for(ivis in sort(unique(knee.saemix$time))) {
  for(itrt in 0:1) {
    yvec<-knee.saemix$y[knee.saemix$treatment==itrt & knee.saemix$time==ivis]
    knee3 <- rbind(knee3,
                   data.frame(time=ivis, treatment=itrt, mean=mean(yvec)))
  }
}

if(FALSE) {
  knee3 <- knee.saemix %>%
    group_by(time, treatment) %>%
    summarise(mean=mean(y))
  
}

simdat <-yfit@sim.data@datasim
simdat$time<-rep(yfit@data@data$time,nsim)
simdat$treatment<-rep(yfit@data@data$treatment,nsim)
ytab<-NULL
for(irep in 1:nsim) {
  xtab<-simdat[simdat$irep==irep,]
  xtab1<-NULL
  for(ivis in sort(unique(xtab$time))) {
    for(itrt in 0:1) {
      yvec<-xtab$ysim[xtab$treatment==itrt & xtab$time==ivis]
      xtab1 <- rbind(xtab1,
                     data.frame(time=ivis, treatment=itrt, mean=mean(yvec)))
    }
  }
  
  # suppressMessages(
  #   xtab1 <- xtab %>%
  #     group_by(time, treatment) %>%
  #     summarise(mean=mean(ysim))
  # )
  ytab<-rbind(ytab,xtab1[,c("time","treatment","mean")])
}
gtab<-NULL
for(ivis in sort(unique(ytab$time))) {
  for(itrt in 0:1) {
    yvec<-ytab$mean[ytab$treatment==itrt & ytab$time==ivis]
    gtab <- rbind(gtab,
                   data.frame(time=ivis, treatment=itrt, lower=quantile(yvec, c(0.05)), mean=median(yvec), upper=quantile(yvec, c(0.95))))
  }
}


# gtab <- ytab %>%
#   group_by(time, treatment) %>%
#   summarise(lower=quantile(mean, c(0.05)), mean=median(mean), upper=quantile(mean, c(0.95)))

kneeMedvpc <- ggplot(data = knee3, aes(x = time, y=mean, group=treatment)) + 
  geom_ribbon(data=gtab, aes(x=time, ymin=lower, ymax=upper), alpha=0.5, fill="lightblue") +
  geom_point(colour='blue') + theme_bw() + 
  scale_fill_brewer(palette = "Blues") + theme(legend.position = "top") +
  labs(fill = "Score") + xlab("Time (d)") + ylab("Median value of score over time") + facet_wrap(.~treatment) 
plot2 <- kneeMedvpc #+ theme(text = element_text(size=rel(rsize.text)), legend.text=element_text(size=rel(rsize.text)), strip.text.x = element_text(size=rel(rsize.ticks)), strip.text.y = element_text(size=rel(rsize.ticks)))

print(kneeMedvpc)
if(saveFigs) {
  namfig<-"knee_medianScoreVPC.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot2)
  dev.off()
}

###################################################### SE
case.ordinal <- read.table(file.path(workDir,"bootstrapRes","bootstrapCase_knee.res"), header=T)
cond.ordinal <- read.table(file.path(workDir,"bootstrapRes","bootstrapCond_knee.res"), header=T)

nboot<-dim(case.ordinal)[1]
case.ordinal <- case.ordinal[!is.na(case.ordinal[,2]),]
cond.ordinal <- cond.ordinal[!is.na(cond.ordinal[,2]),]

par.estim<-format(c(ord.fit@results@fixed.effects,diag(ord.fit@results@omega)[ord.fit@results@indx.omega]), digits=2, nsmall=1)
df2<-data.frame(parameter=colnames(case.ordinal)[-c(1)], saemix=par.estim)

for(i in 1:2) {
  if(i==1) {
    resboot1<-case.ordinal
    namboot<-"case"
  } else {
    resboot1<-cond.ordinal
    namboot <-"cNP"
  }
  mean.bootDist<-apply(resboot1, 2, mean, na.rm=T)[-c(1)]
  sd.bootDist<-apply(resboot1, 2, sd, na.rm=T)[-c(1)]
  quant.bootDist<-apply(resboot1[-c(1)], 2, quantile, c(0.025, 0.975), na.rm=T)
  l1<-paste0(format(mean.bootDist, digits=1, scientific=FALSE)," (",format(sd.bootDist,digits=1, nsmall=0, trim=T, scientific=FALSE),")")
  l2<-paste0("[",format(quant.bootDist[1,], digits=1, trim=T, scientific=FALSE),", ",format(quant.bootDist[2,],digits=1, trim=T, scientific=FALSE),"]")
  df2<-cbind(df2, l1, l2)
  i1<-3+2*(i-1)
  colnames(df2)[i1:(i1+1)]<-paste0(namboot,".",c("estimate","CI"))
}
print(df2)
df2.nocov<-df2

case.ordinal <- read.table(file.path(workDir,"bootstrapRes","bootstrapCase_kneeCov.res"), header=T)
cond.ordinal <- read.table(file.path(workDir,"bootstrapRes","bootstrapCond_kneeCov.res"), header=T)

nboot<-dim(case.ordinal)[1]
case.ordinal <- case.ordinal[!is.na(case.ordinal[,2]),]
cond.ordinal <- cond.ordinal[!is.na(cond.ordinal[,2]),]

par.estim<-format(c(ord.fit.cov@results@fixed.effects,diag(ord.fit.cov@results@omega)[ord.fit.cov@results@indx.omega]), digits=2, nsmall=1)
df2<-data.frame(parameter=colnames(case.ordinal)[-c(1)], saemix=par.estim)
for(i in 1:2) {
  if(i==1) {
    resboot1<-case.ordinal
    namboot<-"case"
  } else {
    resboot1<-cond.ordinal
    namboot <-"cNP"
  }
  mean.bootDist<-apply(resboot1, 2, mean, na.rm=T)[-c(1)]
  sd.bootDist<-apply(resboot1, 2, sd, na.rm=T)[-c(1)]
  quant.bootDist<-apply(resboot1[-c(1)], 2, quantile, c(0.025, 0.975), na.rm=T)
  l1<-paste0(format(mean.bootDist, digits=1, scientific=FALSE)," (",format(sd.bootDist,digits=1, nsmall=0, trim=T, scientific=FALSE),")")
  l2<-paste0("[",format(quant.bootDist[1,], digits=1, trim=T, scientific=FALSE),", ",format(quant.bootDist[2,],digits=1, trim=T, scientific=FALSE),"]")
  df2<-cbind(df2, l1, l2)
  i1<-3+2*(i-1)
  colnames(df2)[i1:(i1+1)]<-paste0(namboot,".",c("estimate","CI"))
}
print(df2)
df2.cov<-df2


df3<-df2.cov[,c(2,6)] # conditional bootstrap
df3<-rbind(df3, c("-", "-"))
rownames(df3)<-c("$\\alpha_1$","$\\beta_{Age^2,\\alpha_1}$","$\\alpha_2$","$\\beta_{Trt,\\alpha_2}$","$\\alpha_3$","$\\alpha_4$","$\\beta$", "$\\beta_{Trt,\\beta}$", "$\\omega_{\\alpha_1}$","$\\omega_{\\alpha_2}$", "$\\omega_{\\alpha_3}$", "$\\omega_{\\alpha_4}$", "$\\omega_{\\beta}$")
df4<-df2.nocov[,c(2,6)] # model without covariates
df3.nocov<-rbind(df4[1,], c("-", "-"), df4[2,],c("-", "-"),df4[3:5,],c("-", "-"),df4[6,],c("-", "-"),c("-", "-"),c("-", "-"),df4[7,])
df3<-cbind(df3.nocov, df3)
rownames(df3)<-c("$\\alpha_1$","$\\beta_{Age^2,\\alpha_1}$","$\\alpha_2$","$\\beta_{Trt,\\alpha_2}$","$\\alpha_3$","$\\alpha_4$","$\\beta$", "$\\beta_{Trt,\\beta}$", "$\\omega_{\\alpha_1}$","$\\omega_{\\alpha_2}$", "$\\omega_{\\alpha_3}$", "$\\omega_{\\alpha_4}$", "$\\omega_{\\beta}$")

print(xtable(df3), only.contents=TRUE, include.rownames=TRUE,  floating=F, sanitize.rownames.function = identity)

