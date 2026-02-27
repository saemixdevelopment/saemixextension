#################################### Setup
# When in my directory - set to FALSE to run from getwd(), must be paperSaemix3 (all paths relative)
EcoAtHome <- FALSE

# Libraries
library(saemix)

library(ggplot2)
library(gridExtra)
library(xtable)
# Automatically loaded when loading saemix (?)
# library(ggplot2)
# library(gridExtra)
# library(MASS)
# library(rlang)

# Folders
if(EcoAtHome) {
  # Load updated files
  # pending compilation and CRAN upload
  workDir<-"/home/eco/work/saemix/saemixextension/paperSaemix3"
} else {
  workDir<-getwd() 
}

setwd(workDir)
bootResDir <- file.path(workDir, "bootstrapRes")
figDir <- file.path(workDir,"figs")

# Not really needed but just to show I can use tidyverse too
# library(tidyverse)
# Weird font issues getting tiny when saving to pdf, correcting for that
#library(showtext)
#showtext_auto()
rsize.text <- 2
rsize.ticks <- 1.5
theme_set(theme_bw(base_size = 15)) 
# note: apparently corrected now  so no longer needed

# Whether to save the plots
saveFigs<-FALSE
figDir <- file.path(workDir,"figs")

# Number of bootstrap samples
runBootstrap <- FALSE # to read the results from disk
nboot <-10 # definitely not enough for bootstrap, please run more samples using the associated script (.R)
# nboot <- 200 # too slow in Rstudio

# Stepwise algorithm for covariates (do not run in Rstudio unless you have a lot more memory than I do on my computer)
runCovKnee<-FALSE

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

##################################################################################################
#################################### Binary response model
#### Data description
data(toenail.saemix)

saemix.data<-saemixData(name.data=toenail.saemix,name.group=c("id"),name.predictors=c("time","y"), name.response="y",
                        name.covariates=c("treatment"), units=list(x="d",y="-"), verbose=FALSE)

#################################### 
#### Exploring data
# Distribution of times
if(FALSE) hist(toenail.saemix$time, breaks=c(-1,0.25,1.25,2.25, 3.25, 7,10,15,20), freq=T)
table(cut(toenail.saemix$time, breaks=c(-1,0.25,1.25,2.25, 3.25, 7,10,15,20)))

# Proportion of 0's and 1's across time
plotDiscreteData(saemix.data, outcome='binary', which.cov="treatment", bin.number=7)

# Only showing the frequency of infections
toe1 <- regroupTable(toenail.saemix, value="y", c("visit","treatment"), fun=list(nev=function(x) sum(x), n=function(x) length(x)))
toe1$freq<-toe1$nev/toe1$n
toe1$sd <- sqrt((1-toe1$nev/toe1$n)*(toe1$nev/toe1$n**2))
toe1$lower <- toe1$freq-1.96*toe1$sd
toe1$upper <- toe1$freq+1.96*toe1$sd
toe1$lower[toe1$lower<0] <-0 # we should use a better approximation for CI
toe1$treatment <- factor(toe1$treatment, labels=c("A","B"))

plot1<-ggplot(toe1, aes(x=visit, y=freq, group=treatment)) + geom_line(aes(colour=treatment)) + 
  geom_point(aes(colour=treatment)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=treatment), alpha=0.2) +
  ylim(c(0,0.5)) + theme_bw() + theme(legend.position = "top") +
  xlab("Visit number") + ylab("Observed frequency of infection")
print(plot1)
plot2 <- plot1 # + theme(text = element_text(size=rel(rsize.text)), legend.text=element_text(size=rel(rsize.text)), strip.text.x = element_text(size=rel(rsize.ticks)), strip.text.y = element_text(size=rel(rsize.ticks)))

#showtext_opts(dpi = 300)
if(saveFigs) {
  namfig<-"toenail_infectionFreq.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot2)
  dev.off()
}
plotToe.infectionFreq <- plot2

#################################### 
#### Statistical model
# saemix model
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

# simulation function (used for diagnostics)
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

saemix.model<-saemixModel(model=binary.model,description="Binary model",
                          simulate.function=simulBinary, modeltype="likelihood",
                          psi0=matrix(c(-0.5,-.15,0,0),ncol=2,byrow=TRUE,
                                      dimnames=list(NULL,c("theta1","theta2"))),
                          transform.par=c(0,0), covariate.model=c(0,1),
                          covariance.model=matrix(c(1,0,0,0),ncol=2), 
                          omega.init=diag(c(0.5,0.3)), verbose=FALSE)

# saemix fit
saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, 
                     displayProgress=FALSE, nb.chains=10, fim=FALSE, print=FALSE)
binary.fit<-saemix(saemix.model,saemix.data,saemix.options)
summary(binary.fit)
plot(binary.fit, plot.type="convergence")

#################################### 
#### Diagnostics 
#  $1_{Y_{ij}=0} \times (1-P(Y_{ij}=1)) + 1_{Y_{ij}=1} \times P(Y_{ij}=1) $
# simulate from model (nsim=100)
nsim<-1000
binary.fit <- simulateDiscreteSaemix(binary.fit, nsim=nsim)

plot1 <- discreteVPC(binary.fit, outcome="binary", which.cov="treatment")

plot2 <- plot1 # + theme(text = element_text(size=rel(rsize.text)), legend.text=element_text(size=rel(rsize.text)), strip.text.x = element_text(size=rel(rsize.ticks)), strip.text.y = element_text(size=rel(rsize.ticks)))
if(saveFigs) {
  namfig<-"toenail_vpcByTreatment.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot2)
  dev.off()
}

plotToe.vpcByTreatment <- plot2

# Nice VPC
simdat <-binary.fit@sim.data@datasim
simdat$visit<-rep(toenail.saemix$visit,nsim)
simdat$treatment<-rep(toenail.saemix$treatment,nsim)
# VPC-type diagnostic
ytab<-NULL
for(irep in 1:nsim) {
  xtab<-simdat[simdat$irep==irep,]
  xtab1 <- regroupTable(xtab, value="ysim", c("visit","treatment"), fun=list(nev=function(x) sum(x), n=function(x) length(x)))
  xtab1<-cbind(xtab1, freq=xtab1[,3]/xtab1[,4])
  ytab<-rbind(ytab,xtab1[,c(1,2,5)])
}

gtab <- regroupTable(ytab, value="freq", c("visit","treatment"), fun=list(lower=function(x) quantile(x, c(0.05)), median=function(x) quantile(x,0.5), upper=function(x) quantile(x,0.95)))
gtab$treatment <- ifelse(gtab$treatment==1,"B","A")
gtab$freq<-1

plot2 <- ggplot(toe1, aes(x=visit, y=freq, group=treatment)) + geom_line(aes(colour=treatment)) + 
  geom_point(aes(colour=treatment)) + 
  geom_line(data=gtab, aes(x=visit, y=median), linetype=2, colour='lightblue') + 
  geom_ribbon(data=gtab,aes(ymin=lower, ymax=upper), alpha=0.5, fill='lightblue') +
  ylim(c(0,0.5)) + theme_bw() + theme(legend.position = "none") + facet_wrap(.~treatment) +
  xlab("Visit number") + ylab("Frequency of infection") 

plot3 <- plot2 # +  theme(text = element_text(size=rel(rsize.text)), legend.text=element_text(size=rel(rsize.text)), strip.text.x = element_text(size=rel(rsize.ticks)), strip.text.y = element_text(size=rel(rsize.ticks)))

print(plot2)
if(saveFigs) {
  namfig<-"toenail_ggplot2VPCTreatment.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot3)
  dev.off()
}
plot.ggplot2VPCTreatment <- plot2


# npd TODO

#################################### 
#### Standard errors of estimation

##### Case bootstrap
if(!runBootstrap)  {
  case.bin <- read.table(file.path(workDir,"bootstrapRes","bootstrapCase_toenail.res"), header=T)
  #  case.bin <- read.table(file.path(saemixDir,"bootstrap","results","toenail_caseBootstrap.res"), header=T)
  nboot<-dim(case.bin)[1]
}  else case.bin <- saemix.bootstrap(binary.fit, method="case", nboot=nboot) 
head(case.bin)

# Bootstrap distributions
if(nboot<200) cat("The number of bootstrap samples is too low to provide good estimates of the confidence intervals\n") else {
  resboot1<-case.bin
  ypd2<-NULL
  for(icol in 1:4) {
    ypd2<-rbind(ypd2,data.frame(rep=resboot1[,1],Param=colnames(resboot1)[(icol+1)],value=resboot1[,(icol+1)], Bootstrap="Case", stringsAsFactors=FALSE))
  }
  
  ypd2$Param<-factor(ypd2$Param, levels = unique(ypd2$Param))
  ypd2.fix<-ypd2[ypd2$Param %in% unique(ypd2$Param)[1:3],]
  ypd2.iiv<-ypd2[ypd2$Param %in% unique(ypd2$Param)[4],]
  ypd <- ypd2
  
  par.estim<-c(binary.fit@results@fixed.effects,diag(binary.fit@results@omega)[binary.fit@results@indx.omega])
  mean.bootDist<-apply(resboot1, 2, mean)[-c(1)]
  df<-data.frame(Param=unique(ypd2$Param), mean.boot=mean.bootDist, est.saemix=par.estim, Bootstrap="Case") 
  
  plot.density2<-ggplot(data=ypd2) + geom_density(aes(value,fill="red4"), alpha=0.5) + 
    geom_vline(data=df,aes(xintercept=est.saemix),colour="red",size=1.2) + 
    geom_vline(data=df,aes(xintercept=mean.boot),colour="blue",size=1.2) +
    theme_bw() + theme(axis.title.x = element_blank(),axis.text.x = element_text(size=11, angle=30, hjust=1), legend.position = "none") + 
    facet_wrap(~Param, ncol=2, scales = 'free')
#  print(plot.density2)
}

plot.densityCaseToenail <- plot.density2

##### Conditional bootstrap
if(!runBootstrap) {
  #  cond.bin <- read.table(file.path(saemixDir,"bootstrap","results","toenail_condBootstrap.res"), header=T)
  cond.bin <- read.table(file.path(workDir,"bootstrapRes","bootstrapCond_toenail.res"), header=T)
  nboot<-dim(cond.bin)[1]
} else 
  cond.bin <- saemix.bootstrap(binary.fit, method="conditional", nboot=nboot) 
summary(cond.bin)

# Bootstrap distributions
if(nboot<200) cat("The number of bootstrap samples is too low to provide good estimates of the confidence intervals\n") else {
  resboot1<-cond.bin
  ypd2<-NULL
  for(icol in 1:4) {
    ypd2<-rbind(ypd2,data.frame(rep=resboot1[,1],Param=colnames(resboot1)[(icol+1)],value=resboot1[,(icol+1)], Bootstrap="Conditional", stringsAsFactors=FALSE))
  }
  
  ypd2$Param<-factor(ypd2$Param, levels = unique(ypd2$Param))
  ypd2.fix<-ypd2[ypd2$Param %in% unique(ypd2$Param)[1:3],]
  ypd2.iiv<-ypd2[ypd2$Param %in% unique(ypd2$Param)[4],]
  ypd <- rbind(ypd,ypd2)
  
  par.estim<-c(binary.fit@results@fixed.effects,diag(binary.fit@results@omega)[binary.fit@results@indx.omega])
  mean.bootDist<-apply(resboot1, 2, mean)[-c(1)]
  df2<-data.frame(Param=unique(ypd2$Param), mean.boot=mean.bootDist, est.saemix=par.estim, Bootstrap="Conditional")
  df<-rbind(df,df2)
  
  plot.density2<-ggplot(data=ypd2) + geom_density(aes(value,fill="red4"), alpha=0.5) + 
    geom_vline(data=df2,aes(xintercept=est.saemix),colour="red",size=1.2) + 
    geom_vline(data=df2,aes(xintercept=mean.boot),colour="blue",size=1.2) +
    theme_bw() + theme(axis.title.x = element_blank(),axis.text.x = element_text(size=11, angle=30, hjust=1), legend.position = "none") + 
    facet_wrap(~Param, ncol=2, scales = 'free')
  print(plot.density2)
  
  plot.density3<-ggplot(data=ypd) + geom_density(aes(value,fill="red4"), alpha=0.5) + 
    geom_vline(data=df,aes(xintercept=est.saemix),colour="red",size=1.2) + 
    geom_vline(data=df,aes(xintercept=mean.boot),colour="blue",size=1.2) +
    theme_bw() + theme(axis.title.x = element_blank(),axis.text.x = element_text(size=11, angle=30, hjust=1), legend.position = "none") + 
    facet_grid(Bootstrap~Param, scales = 'free')
  #    facet_wrap(Bootstrap~Param, nrow=2, scales = 'free')
  
  print(plot.density3)
}
plot.densityBootToenail <- plot.density3

##### Bootstrap results
if(nboot<200) cat("The number of bootstrap samples is too low to provide good estimates of the confidence intervals\n") else {
  par.estim<-c(binary.fit@results@fixed.effects,diag(binary.fit@results@omega)[binary.fit@results@indx.omega], sqrt(diag(binary.fit@results@omega)[binary.fit@results@indx.omega]))
  namsd<-paste0("SD.",colnames(binary.fit@results@omega)[binary.fit@results@indx.omega])
  df2<-data.frame(parameter=c(colnames(case.bin)[-c(1)],namsd), saemix=par.estim)
  for(i in 1:2) {
    if(i==1) {
      resboot1<-case.bin
      namboot<-"case"
    } else {
      resboot1<-cond.bin
      namboot <-"cNP"
    }
    
    ncol<-dim(resboot1)[2]
    resboot1<-cbind(resboot1, sqrt(resboot1[,ncol]))
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
}

rownames(df2)<-c("$\\theta_1$ (-)","$\\theta_2$ (mo$^{(-1)}$)","$\\beta_{trt,\\theta_2}$ (-)", "$\\omega^2_{\\theta_1}$ (-)","$\\omega_{\\theta_1}$ (-)")
print(xtable(df2),only.contents=TRUE, include.rownames=TRUE)

xtab.toenailBootSE <- df2
# df2 <- df2[-c(4),-c(1)]

##################################################################################################
#################################### Categorical response model

#### Data
data(knee.saemix)

# Data
ordknee.data<-saemixData(name.data=knee.saemix,name.group=c("id"),
                         name.predictors=c("y", "time"), name.X=c("time"),
                         name.covariates = c("Age","Sex","treatment","Age2"),
                         units=list(x="d",y="", covariates=c("yr","-","-","yr2")), verbose=FALSE)

plotDiscreteData(ordknee.data, outcome="categorical", which.cov="treatment")

plotDiscreteData(ordknee.data, outcome="categorical", which.cov="Sex")
plot1 <- plotDiscreteData(ordknee.data, outcome="categorical", which.cov="treatment")
plot2 <- plot1 #+  theme(text = element_text(size=rel(rsize.text)), legend.text=element_text(size=rel(rsize.text)), strip.text.x = element_text(size=rel(rsize.ticks)), strip.text.y = element_text(size=rel(rsize.ticks)))

if(saveFigs) {
  namfig<-"knee_rawDataPropTime.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot2)
  dev.off()
}

kneepl.rawDataPropTime <-plot2

# Alternative representation as barplots
gtab <- regroupTable(knee.saemix, value="y", c("time","y"), fun=list(n=function(x) length(x)))
gtab$y <- factor(gtab$y)

plot1 <- ggplot(data = gtab, aes(x = time, y=n, group=y, fill=y)) + 
  geom_bar(stat="identity", position = "dodge") + theme_bw() + 
  scale_fill_brewer(palette = "Reds") + theme(legend.position = "top") +
  labs(fill = "Score") + xlab("Time (d)") + ylab("Counts")
print(plot1)
kneepl.barplotScoreTime <-plot1

# Gender effect ?
if(FALSE) {
  gtab2 <- regroupTable(knee.saemix, value="y", c("time","Sex"), fun=list(meanscore=function(x) mean(x)))
  
  ggplot(data = gtab2, aes(x = time, y=meanscore, group=Sex, fill=Sex)) + 
    geom_line() + theme_bw() + theme(legend.position = "top") +
    labs(fill = "Score") + xlab("Time (d)") + ylab("Mean score")
}

#################################### 
#### Model

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

# Saemix model
saemix.model<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                          simulate.function=simulateOrdinal, psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5, byrow=TRUE, 
                                                                         dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), transform.par=c(0,1,1,1,1),
                          omega.init=diag(c(100, 1, 1, 1, 1)), covariance.model = diag(c(1,0,0,0,1)), verbose=FALSE)

# Fitting
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, nb.chains=10, nbiter.saemix=c(600,100), print=FALSE)
#saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, nb.chains=10, fim=FALSE)

ord.fit<-saemix(saemix.model,ordknee.data,saemix.options)
summary(ord.fit)
plot(ord.fit, plot.type="convergence")

## Note: comparable estimates obtained with Monolix (not same, but within CI)
## quite a lot of sensitivity to distributions (when using eg normal distributions in Monolix the parameters and most importantly the SE's fluctuated quite a bit)

#################################### 
# Covariate model
# Do not run, Rstudio fails (ran in a script as "R CMD BATCH paper_kneeCovModel.R paper_kneeCovModel.out")
#if(runCovKnee) cov.ordfit <- step.saemix(ord.fit, trace=TRUE, direction='both')

# Resulting model
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

# Compare the base and covariate model 
compare.saemix(ord.fit, ord.fit.cov)

#################################### 
#### Model evaluation

### Simulations for VPC
nsim<-100
yfit<-ord.fit.cov
yfit<-simulateDiscreteSaemix(yfit, nsim=nsim)
discreteVPC(yfit, outcome="categorical")
plot1 <- discreteVPC(yfit, outcome='categorical',covsplit=TRUE, which.cov="treatment")
plot2 <-  discreteVPC(yfit, outcome='categorical',covsplit=TRUE, which.cov="Sex")

if(saveFigs) {
  #  plot1 <- plot1 +  theme(text = element_text(size=rel(rsize.text)), legend.text=element_text(size=rel(rsize.text)), strip.text.x = element_text(size=rel(rsize.ticks)), strip.text.y = element_text(size=rel(rsize.ticks)))
  namfig<-"knee_VPCbytreatment.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot2)
  dev.off()
  
  #  plot2 <- plot2 +  theme(text = element_text(size=rel(rsize.text)), legend.text=element_text(size=rel(rsize.text)), strip.text.x = element_text(size=rel(rsize.ticks)), strip.text.y = element_text(size=rel(rsize.ticks)))
  namfig<-"knee_VPCbySex.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot2)
  dev.off()
}

print(plot1)
print(plot2)
kneepl.VPCbytreatment <- plot1
kneepl.VPCbySex <- plot2


if(FALSE) {
  ### Simulations for VPC - model without covariates
  nsim<-100
  yfit<-ord.fit
  yfit<-simulateDiscreteSaemix(yfit, nsim=nsim)
  discreteVPC(yfit, outcome="categorical")
}

####################################
# VPC for mean score in each group
knee3 <- regroupTable(knee.saemix, value="y", c("time","treatment"), fun=list(mean=function(x) mean(x)))

simdat <-yfit@sim.data@datasim
simdat$time<-rep(yfit@data@data$time,nsim)
simdat$treatment<-rep(yfit@data@data$treatment,nsim)
# VPC-type diagnostic
ytab<-NULL
for(irep in 1:nsim) {
  xtab<-simdat[simdat$irep==irep,]
  xtab1 <- regroupTable(xtab, value="ysim", c("time","treatment"), fun=list(mean=function(x) mean(x)))
  ytab<-rbind(ytab,xtab1)
}

gtab <- regroupTable(ytab, value="mean", c("time","treatment"), fun=list(lower=function(x) quantile(x, c(0.05)), mean=function(x) quantile(x,0.5), upper=function(x) quantile(x,0.95)))

kneeMedvpc <- ggplot(data = knee3, aes(x = time, y=mean, group=treatment)) + 
  geom_ribbon(data=gtab, aes(x=time, ymin=lower, ymax=upper), alpha=0.5, fill="lightblue") +
  geom_point(colour='blue') + theme_bw() + scale_fill_brewer(palette = "Blues")  +
  xlab("Time (d)") + ylab("Mean value of score over time") + facet_wrap(.~treatment)
plot1 <- kneeMedvpc # +  theme(text = element_text(size=rel(rsize.text)), legend.text=element_text(size=rel(rsize.text)), strip.text.x = element_text(size=rel(rsize.ticks)), strip.text.y = element_text(size=rel(rsize.ticks)))
print(kneeMedvpc)

if(saveFigs) {
  namfig<-"knee_meanScoreVPC.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot1)
  dev.off()
}
kneepl.VPCmeanScore <- kneeMedvpc

#################################### 
#### Estimation errors

##### Model without covariates
if(runBootstrap) {
  case.ordinal <- saemix.bootstrap(ord.fit, method="case", nboot=nboot) 
  cond.ordinal <- saemix.bootstrap(ord.fit, method="conditional", nboot=nboot) 
} else {
  case.ordinal <- read.table(file.path(workDir,"bootstrapRes","bootstrapCase_knee.res"), header=T)
  cond.ordinal <- read.table(file.path(workDir,"bootstrapRes","bootstrapCond_knee.res"), header=T)
}
case.ordinal <- case.ordinal[!is.na(case.ordinal[,2]),]
cond.ordinal <- cond.ordinal[!is.na(cond.ordinal[,2]),]
nboot<-dim(case.ordinal)[1]
for(icol in 7:8) {
  case.ordinal[,icol]<-sqrt(case.ordinal[,icol]) # compute SD (omega) from omega2
  cond.ordinal[,icol]<-sqrt(cond.ordinal[,icol])
}

par.estim<-format(c(ord.fit@results@fixed.effects,sqrt(diag(ord.fit@results@omega)[ord.fit@results@indx.omega])), digits=2, nsmall=1)
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
  l1<-paste0(format(mean.bootDist, digits=2)," (",format(sd.bootDist,digits=2, trim=T),")")
  l2<-paste0("[",format(quant.bootDist[1,], digits=2),", ",format(quant.bootDist[2,],digits=2, trim=T),"]")
  df2<-cbind(df2, l1, l2)
  i1<-3+2*(i-1)
  colnames(df2)[i1:(i1+1)]<-paste0(namboot,".",c("estimate","CI"))
}
print(df2)
df3<-df2[,c(2,3,4)] # case bootstrap
rownames(df3)<-c("$\\theta_1$","$\\theta_2$","$\\theta_3$","$\\theta_4$","$\\alpha$", "$\\omega_{\\theta_1}$","$\\omega_{\\alpha}$")
# print(xtable(df3), only.contents=TRUE, include.rownames=TRUE,  floating=F, sanitize.rownames.function = identity)

df3<-df2[,c(2,5,6)] # conditional bootstrap
rownames(df3)<-c("$\\theta_1$","$\\theta_2$","$\\theta_3$","$\\theta_4$","$\\alpha$", "$\\omega_{\\theta_1}$","$\\omega_{\\alpha}$")
# print(xtable(df3), only.contents=TRUE, include.rownames=TRUE,  floating=F, sanitize.rownames.function = identity)
xtab.kneeBase <- df3

cat("Conditional bootstrap estimates computed with B=",dim(cond.ordinal)[1],"samples\n")
print(df3)

###### Knee model with covariates
if(runBootstrap) {
  case.ordinal <- saemix.bootstrap(ord.fit.cov, method="case", nboot=nboot) 
  cond.ordinal <- saemix.bootstrap(ord.fit.cov, method="conditional", nboot=nboot) 
} else {
  case.ordinal <- read.table(file.path(workDir,"bootstrapRes","bootstrapCase_kneeCov.res"), header=T)
  cond.ordinal <- read.table(file.path(workDir,"bootstrapRes","bootstrapCond_kneeCov.res"), header=T)
}
case.ordinal <- case.ordinal[!is.na(case.ordinal[,2]),]
cond.ordinal <- cond.ordinal[!is.na(cond.ordinal[,2]),]
nboot<-dim(case.ordinal)[1]
for(icol in 10:13) {
  case.ordinal[,icol]<-sqrt(case.ordinal[,icol]) # compute SD (omega) from omega2
  cond.ordinal[,icol]<-sqrt(cond.ordinal[,icol])
}


par.estim<-format(c(ord.fit.cov@results@fixed.effects,sqrt(diag(ord.fit.cov@results@omega)[ord.fit.cov@results@indx.omega])), digits=2, nsmall=1)
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
  l1<-paste0(format(mean.bootDist, digits=2)," (",format(sd.bootDist,digits=2, trim=T),")")
  l2<-paste0("[",format(quant.bootDist[1,], digits=2),", ",format(quant.bootDist[2,],digits=2, trim=T),"]")
  df2<-cbind(df2, l1, l2)
  i1<-3+2*(i-1)
  colnames(df2)[i1:(i1+1)]<-paste0(namboot,".",c("estimate","CI"))
}
print(df2)

# pretty print for paper
df3<-df2[,c(2,4)] # case bootstrap
rownames(df3)<-c("$\\theta_1$","$\\beta_{Age2,\\theta_1}$","$\\theta_2$","$\\beta_{trt,\\theta_2}$","$\\theta_3$","$\\theta_4$","$\\alpha$","$\\beta_{trt,\\alpha}$",  "$\\omega_{\\theta_1}$", "$\\omega_{\\theta_2}$",  "$\\omega_{\\theta_3}$",  "$\\omega_{\\theta_4}$")
print(xtable(df3), only.contents=TRUE, include.rownames=TRUE,  floating=F, sanitize.rownames.function = identity)

df3<-df2[,c(2,5,6)] # conditional bootstrap
rownames(df3)<-c("$\\theta_1$","$\\beta_{Age2,\\theta_1}$","$\\theta_2$","$\\beta_{trt,\\theta_2}$","$\\theta_3$","$\\theta_4$","$\\alpha$","$\\beta_{trt,\\alpha}$",  "$\\omega_{\\theta_1}$", "$\\omega_{\\theta_2}$",  "$\\omega_{\\theta_3}$",  "$\\omega_{\\theta_4}$")
print(xtable(df3), only.contents=TRUE, include.rownames=TRUE,  floating=F, sanitize.rownames.function = identity)
xtab.kneeCov <-df3

##################################### Additional analyses for knee data - interaction therapy/gender

if(FALSE) {
  knee.saemix$intSexTh <- ifelse(knee.saemix$Sex==1 & knee.saemix$treatment==1,1,0)
  
  ordknee.data<-saemixData(name.data=knee.saemix,name.group=c("id"),
                           name.predictors=c("y", "time"), name.X=c("time"),
                           name.covariates = c("Age","Sex","treatment","Age2","intSexTh"),
                           units=list(x="d",y="", covariates=c("yr","-","-","yr2","-")))
  plotDiscreteData(ordknee.data, outcome="categorical", which.cov="intSexTh")
  plotDiscreteData(ordknee.data, outcome="categorical", which.cov="Sex")
  plotDiscreteData(ordknee.data, outcome="categorical", which.cov="treatment")
  
  covmodel4<-covmodel2
  covmodel4[3,]<-1
  
  saemix.model.cov4<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",simulate.function=simulateOrdinal,
                                 psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                                 transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)),
                                 covariate.model = covmodel4)
  ord.fit.cov4<-saemix(saemix.model.cov4,ordknee.data,saemix.options)
  
  ### Simulations for VPC
  nsim<-100
  yfit<-ord.fit.cov4
  yfit<-simulateDiscreteSaemix(yfit, nsim=nsim)
  
  ### VPC
  discreteVPC(yfit, outcome="categorical")
  discreteVPC(yfit, outcome="categorical", covsplit=TRUE, which.cov="Sex")
  discreteVPC(yfit, outcome="categorical", covsplit=TRUE, which.cov="treatment")
  
  covmodel5<-rbind(covmodel4, c(1,1,1,0,0))
  
  saemix.model.cov5<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",simulate.function=simulateOrdinal,
                                 psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                                 transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)),
                                 covariate.model = covmodel5)
  ord.fit.cov5<-saemix(saemix.model.cov5,ordknee.data,saemix.options)
  ### Simulations for VPC
  nsim<-100
  yfit<-ord.fit.cov5
  yfit<-simulateDiscreteSaemix(yfit, nsim=nsim)
  
  ### VPC
  discreteVPC(yfit, outcome="categorical")
  discreteVPC(yfit, outcome="categorical", covsplit=TRUE, which.cov="Sex")
  discreteVPC(yfit, outcome="categorical", covsplit=TRUE, which.cov="treatment")
  
}
