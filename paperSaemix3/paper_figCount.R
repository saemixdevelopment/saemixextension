# Execute all models and save figures
# Tables in LaTeX format are output but not saved 
library(xtable)

source("saemix3_countModel.R")

###################################################### Main text (arXiv)
# Figure 4
# Diagnostics ZIPoisson model
namfig<-"rapi_diagnosZIP.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
#  discreteVPC(ysim.zip2, outcome="count", breaks=c(0:9,16,25,80), which.cov="gender")
plot(plotRapi.diagnosZIP)
dev.off()

###################################################### Supplementary (arXiv)
# Raw RAPI data
namfig<-"rapi_rawDataPropTime.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(plotRapi.rawDataPropTime)
dev.off()

# VPC pour Poisson model
namfig<-"rapi_poissonVPC.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(plotRapi.poissonVPC)
dev.off()

# VPC for zero counts in Poisson model
namfig<-"rapi_poissonZeroesVPC.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(plotRapi.poissonZeroesVPC)
dev.off()

# Comparison between the proportion of zeroes for the different models
namfig<-"rapi_comparePropZeroes.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
print(plotRapi.comparePropZeroes)
dev.off()

# Diagnostics for truncated Poisson in the hurdle model
namfig<-"rapi_hurdleVPC_truncatedPoiss.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(plotRapi.hurdleVPC.truncatedPoiss)
dev.off()

namfig<-"rapi_hurdleVPC_truncatedPoiss.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(plotRapi.hurdleVPC.truncatedPoiss)
dev.off()

###################################################### Tables
### Table 3 (main text)
### Parameter estimates for the Poisson, Zero-Inflated Poisson and logistic portion of the hurdle model 

print(xtable(rr.tab), only.contents=TRUE, include.rownames=T, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity)

### Table comparing parameter estimates with their SE for paper (Supplementary)
### SE computed using conditional bootstrap

# Poisson with covariates
cond.count <- read.table(file.path(bootResDir,"bootstrapCond_rapiCov_Poisson.res"), header=T)
nboot<-dim(cond.count)[1]
cond.count <- cond.count[!is.na(cond.count[,2]),]
# compute correlation instead of covariance
ncol<-dim(cond.count)[2]
cond.count[ncol]<-cond.count[ncol]/sqrt(cond.count[ncol-1]*cond.count[ncol-2])

corr<-poisson.fit.cov2@results@omega[1,2]/sqrt(poisson.fit.cov2@results@omega[1,1]*poisson.fit.cov2@results@omega[2,2])
par.estim<-c(poisson.fit.cov2@results@fixed.effects,diag(poisson.fit.cov2@results@omega), corr)
df2<-data.frame(parameter=colnames(cond.count)[-c(1)])
df2[dim(df2)[1],1]<-"corr.intSlope"
resboot1<-cond.count
sd.bootDist<-apply(resboot1, 2, sd)[-c(1)]
l1<-paste0(format(par.estim, digits=1,nsmall=0, scientific=FALSE)," (",format(sd.bootDist,digits=1, nsmall=0, trim=T, scientific=FALSE),")")
quant.bootDist<-apply(resboot1[-c(1)], 2, quantile, c(0.025, 0.975))
l2<-paste0("[",format(quant.bootDist[1,], digits=2),", ",format(quant.bootDist[2,],digits=2, trim=T),"]")
df2<-cbind(df2,l1)
colnames(df2)[2]<-"Poisson"
df2<-rbind(df2[1:4,],df2[4,],df2[5:7,])
df2[5,]<-c("p0","-")

# ZIPoisson with covariates
cond.count <- read.table(file.path(bootResDir,"bootstrapCond_rapiCov_ZIP_corr.res"), header=T)
nboot<-dim(cond.count)[1]
cond.count <- cond.count[!is.na(cond.count[,2]),]
# compute correlation instead of covariance
ncol<-dim(cond.count)[2]
cond.count[ncol]<-cond.count[ncol]/sqrt(cond.count[ncol-1]*cond.count[ncol-2])

corr<-zippoisson.fit.cov3@results@omega[1,2]/sqrt(zippoisson.fit.cov3@results@omega[1,1]*zippoisson.fit.cov3@results@omega[2,2])
par.estim<-c(zippoisson.fit.cov3@results@fixed.effects,diag(zippoisson.fit.cov3@results@omega)[zippoisson.fit.cov3@results@indx.omega], corr)
resboot1<-cond.count
sd.bootDist<-apply(resboot1, 2, sd)[-c(1)]
l1<-paste0(format(par.estim, digits=1,nsmall=0, scientific=FALSE)," (",format(sd.bootDist,digits=1, nsmall=0, trim=T, scientific=FALSE),")")
df2<-cbind(df2,l1)
colnames(df2)[3]<-"ZIP"

# Hurdle model
## binary logistic component
cond.count <- read.table(file.path(bootResDir,"bootstrapCond_rapiCov_Hurdle0_corr.res"), header=T)
nboot<-dim(cond.count)[1]
cond.count <- cond.count[!is.na(cond.count[,2]),]
# compute correlation instead of covariance
ncol<-dim(cond.count)[2]
cond.count[ncol]<-cond.count[ncol]/sqrt(cond.count[ncol-1]*cond.count[ncol-2])
corr<-hurdlefit0@results@omega[1,2]/sqrt(hurdlefit0@results@omega[1,1]*hurdlefit0@results@omega[2,2])
par.estim<-c(hurdlefit0@results@fixed.effects,diag(hurdlefit0@results@omega)[hurdlefit0@results@indx.omega], corr)

resboot1<-cond.count
sd.bootDist<-apply(resboot1, 2, sd)[-c(1)]
l1<-paste0(format(par.estim, digits=1,nsmall=0, scientific=FALSE)," (",format(sd.bootDist,digits=1, nsmall=0, trim=T, scientific=FALSE),")")
df2<-cbind(df2,c(l1[1:4],"-",l1[5:7]))
colnames(df2)[4]<-"Hurdle - binary"

## truncated Poisson component
cond.count <- read.table(file.path(bootResDir,"bootstrapCond_rapiCov_Hurdle1_corr.res"), header=T)
nboot<-dim(cond.count)[1]
cond.count <- cond.count[!is.na(cond.count[,2]),]
# compute correlation instead of covariance
ncol<-dim(cond.count)[2]
cond.count[ncol]<-cond.count[ncol]/sqrt(cond.count[ncol-1]*cond.count[ncol-2])

corr<-hurdlefit1@results@omega[1,2]/sqrt(hurdlefit1@results@omega[1,1]*hurdlefit1@results@omega[2,2])
par.estim<-c(hurdlefit1@results@fixed.effects,diag(hurdlefit1@results@omega)[hurdlefit1@results@indx.omega], corr)

resboot1<-cond.count
sd.bootDist<-apply(resboot1, 2, sd)[-c(1)]
l1<-paste0(format(par.estim, digits=1,nsmall=0, scientific=FALSE)," (",format(sd.bootDist,digits=1, nsmall=0, trim=T, scientific=FALSE),")")
df2<-cbind(df2,c(l1[1:4],"-",l1[5:7]))
colnames(df2)[5]<-"Hurdle - truncated"

df2<-rbind(c("","Value (SE)","Value (SE)","Value (SE)","Value (SE)"), df2)
df3<-df2[,-c(1)]
rownames(df3)<-c("Parameter","$\\alpha_0$ (-)","$\\beta_{male,\\alpha_0}$ (-)","$\\alpha_1$ (mo$^{(-1)}$)","$\\beta_{male,\\alpha_1}$ (-)","$p_0$ (-)", "$\\omega_{\\alpha_0}$ (-)","$\\omega_{\\alpha_1}$ (-)", "$\\rho(\\alpha_0,\\alpha_1)$ (-)")

print(xtable(df3), only.contents=TRUE, include.rownames=TRUE,  floating=F, sanitize.rownames.function = identity)

######################################################  Proportion of zeroes
if(FALSE){ # old models
  nsim<-100
  yfit1<-simulateDiscreteSaemix(poisson.fit.cov2, nsim=nsim)
  obs.prop0 <-  length(yfit1@data@data$rapi[yfit1@data@data$rapi==0])/yfit1@data@ntot.obs
  poiss.prop0 <-length(yfit1@sim.data@datasim$ysim[yfit1@sim.data@datasim$ysim==0])/length(yfit1@sim.data@datasim$ysim)
  cat("Observed proportion of 0's",obs.prop0,"\n")
  cat("      Poisson model, p=",poiss.prop0,"\n")
  
  ysim.zip<-simulateDiscreteSaemix(zippoisson.fit.cov2, nsim=nsim)
  zip.prop0 <-length(ysim.zip@sim.data@datasim$ysim[ysim.zip@sim.data@datasim$ysim==0])/length(ysim.zip@sim.data@datasim$ysim)
  cat("      ZI-Poisson model, p=",zip.prop0,"\n")
  
  yhurdle0<-simulateDiscreteSaemix(hurdlefit0, nsim=nsim)
  hurd.prop0 <-length(yhurdle0@sim.data@datasim$ysim[yhurdle0@sim.data@datasim$ysim==0])/length(yhurdle0@sim.data@datasim$ysim)
}
###################################################### Quick check
if(FALSE) { # different results from Atkins for the hurdle, trying to figure out why
  tab.men <- table(rapi.saemix$rapi[rapi.saemix$time==0 & rapi.saemix$gender==1] == 0)
  tab.women <- table(rapi.saemix$rapi[rapi.saemix$time==0 & rapi.saemix$gender==0] == 0)
  cat("Observed proportion of 0's at time 0 (in women):",tab.women[2]/sum(tab.women),"\n")
  
  xcal<-rnorm(1000, mean=hurdlefit0@results@fixed.effects[1], sd=sqrt(hurdlefit0@results@omega[1,1]))
  cat("Expected proportion of 0's at time 0 (in women):",mean(1-1/(1+exp(-xcal))),"\n")
  
  xcal1<-rnorm(1000, mean=1.56, sd=sqrt(2.30766)) # omega not reported in Atkins et al.
  cat("Expected proportion of 0's at time 0 (in women) from Atkins, assuming same variance:",mean(1-1/(1+exp(-xcal1))),"\n")
  xcal2<-rnorm(1000, mean=1.56, sd=sqrt(0.01)) # omega not reported in Atkins et al.
  mean(1-1/(1+exp(-xcal2)))
}

