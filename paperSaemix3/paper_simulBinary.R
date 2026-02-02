# Loading libraries
library(xtable)
library(ggplot2)

# Loading saemix
library(saemix)

# Folders
# workDir<-"/home/eco/work/saemix/saemixextension/paperSaemix3"
workDir <- getwd() # working directory should be paperSaemix3
# setwd(workDir)

#saemixDir<-"/home/eco/work/saemix/saemixextension"
#paperDir <- "/home/eco/xtex/saemix/saemix3"
simDir <- file.path(workDir, "discreteEval")
figDir <- file.path(workDir, "figs")
nsim<-200 # Number of simulations
saveFigure <- TRUE

###################################################### Data
# Proportion of 0's and 1's across time
data(toenail.saemix)

saemix.data<-saemixData(name.data=toenail.saemix,name.group=c("id"),name.predictors=c("time","y"), name.response="y",
                        name.covariates=c("treatment"), verbose=FALSE)

plotDiscreteData(saemix.data, outcome='binary', which.cov="treatment")

# tidyverse version (yes, I had a phase)
if(FALSE) {
  toe1 <- toenail.saemix %>%
    group_by(visit, treatment) %>%
    summarise(nev = sum(y), n=n()) %>%
    mutate(freq = nev/n, sd=sqrt((1-nev/n)*(nev/n**2))) %>%
    mutate(lower=freq-1.96*sd, upper=freq+1.96*sd)
  toe1$lower[toe1$lower<0] <-0 # we should use a better approximation for CI
  toe1$treatment <- factor(toe1$treatment, labels=c("A","B"))
}
toe1<-NULL
for(ivis in sort(unique(toenail.saemix$visit))) {
  for(itrt in 0:1) {
    yvec<-toenail.saemix$y[toenail.saemix$treatment==itrt & toenail.saemix$visit==ivis]
    toe1 <- rbind(toe1,
                  data.frame(visit=ivis, treatment=itrt, nev=sum(yvec), n=length(yvec)))
  }
}
toe1$freq <- toe1$nev/toe1$n
toe1$sd <- sqrt((1-toe1$freq)*(toe1$freq)/toe1$n)  # sqrt((1-toe1$freq)*(toe1$freq**2))
toe1$lower <- toe1$freq -1.96*toe1$sd
toe1$upper <- toe1$freq +1.96*toe1$sd
toe1$lower[toe1$lower<0] <-0 # we should use a better approximation for CI
toe1$treatment <- factor(toe1$treatment, labels=c("A","B"))

plot1<-ggplot(toe1, aes(x=visit, y=freq, group=treatment)) + geom_line(aes(colour=treatment)) + 
  geom_point(aes(colour=treatment)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=treatment), alpha=0.2) +
  ylim(c(0,0.5)) + theme_bw() + theme(legend.position = "top") +
#  theme(legend.position = "top", axis.text=element_text(size=15), legend.text=element_text(size=15),axis.title=element_text(size=17,face="bold")) +
  xlab("Visit number") + ylab("Observed frequency of infection")

#plot1 <- plot1+theme(text = element_text(size=rel(4)), legend.text=element_text(size=rel(4)))

print(plot1)

if(saveFigure) {
  namfig<-"toenail_infectionFreq.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot1)
  dev.off()
}

###################################################### REE

reepl <- reese <- NULL
for(iscenar in 1:2) {
  if(iscenar==1) {
    # Original
    scenario <- "binaryOrig"
    datDir <- file.path(simDir, scenario, "data")
    resDir <- file.path(simDir, scenario, "results")
    par.fix<-c(-1.71, -0.39, -0.15)
    par.om <- c(4.02)
  }
  if(iscenar==2) {
    # IIV
    scenario <- "binaryIIV"
    datDir <- file.path(simDir, scenario, "data")
    resDir <- file.path(simDir, scenario, "results")
    par.fix<-c(-1.71, -0.39, -0.15)
    par.om <- c(1, 0.2) # around 50% variability for both
  }
  # Parameters
  par.true <- c(par.fix, par.om)
  npar <- length(par.true)
  
  namResFiletrue <- file.path(resDir, paste0("truefit_",scenario,".tab"))
  namResFilepop <- file.path(resDir, paste0("popfit_",scenario,".tab"))
  namResFilefar <- file.path(resDir, paste0("farfit_",scenario,".tab"))
  
  tab.true <- read.table(namResFiletrue, header=T)
  tab.pop <- read.table(namResFilepop, header=T)
  tab.far <- read.table(namResFilefar, header=T)
  
  ree<-tab.true[,c(2:(1+npar))]
  for(i in 1:npar) ree[,i]<-(ree[,i]-par.true[i])/par.true[i]
  x1<-apply(ree,2,mean)
  x2<-apply(ree,2,sd)
  l1<-paste0(format(x1*100, digits=2, nsmall=2)," (",format(x2*100, digits=1, nsmall=0),")")
  if(iscenar==1) l1 <-c(l1,c(""))
  reese <- cbind(reese,l1)
  ree <- ree*100
  if(FALSE) {
    reepl <- rbind(reepl,
                   cbind(ree %>%
                           gather(key="variable",value="value"),
                         setting="True", scenario=scenario))
  }
  for(icol in 1:dim(ree)[2])
    reepl <- rbind(reepl,
                   data.frame(value=ree[,icol],variable=colnames(ree)[icol],setting="True", scenario=scenario))
  
  ree<-tab.pop[,c(2:(1+npar))]
  for(i in 1:npar) ree[,i]<-(ree[,i]-par.true[i])/par.true[i]
  x1<-apply(ree,2,mean)
  x2<-apply(ree,2,sd)
  l1<-paste0(format(x1*100, digits=2, nsmall=2)," (",format(x2*100, digits=1, nsmall=0),")")
  if(iscenar==1) l1 <-c(l1,c(""))
  reese <- cbind(reese,l1)
  ree <- ree*100
  
  if(FALSE) {
    reepl <- rbind(reepl,
                 cbind(ree %>%
                         gather(key="variable",value="value"), setting="Pop", scenario=scenario))
  }
  for(icol in 1:dim(ree)[2])
    reepl <- rbind(reepl,
                   data.frame(value=ree[,icol],variable=colnames(ree)[icol],setting="Pop", scenario=scenario))
  
  ree<-tab.far[,c(2:(1+npar))]
  for(i in 1:npar) ree[,i]<-(ree[,i]-par.true[i])/par.true[i]
  x1<-apply(ree,2,mean)
  x2<-apply(ree,2,sd)
  l1<-paste0(format(x1*100, digits=2, nsmall=2)," (",format(x2*100, digits=1, nsmall=0),")")
  if(iscenar==1) l1 <-c(l1,c(""))
  reese <- cbind(reese,l1)
  ree <- ree*100
  
  if(FALSE) {
    reepl <- rbind(reepl,
                 cbind(ree %>%
                         gather(key="variable",value="value"), setting="Far", scenario=scenario))
  }
  for(icol in 1:dim(ree)[2])
    reepl <- rbind(reepl,
                   data.frame(value=ree[,icol],variable=colnames(ree)[icol],setting="Far", scenario=scenario))
  
}
reepl$scenario<-factor(reepl$scenario, levels=c("binaryOrig","binaryIIV"), labels=c("Original model","IIV on both parameters"))
reepl$variable <- factor(reepl$variable, levels=c("theta1","theta2", "beta","omega1","omega2"))

ggplot(reepl, aes(x=variable, y=value)) + geom_violin() + 
  geom_abline(intercept=0, slope=0) + geom_abline(intercept=10, slope=0, linetype="dashed") + geom_abline(intercept=-10, slope=0, linetype="dashed") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  facet_grid(scenario~setting)

reepl1 <- reepl[reepl$setting=="True",]
plot1 <- ggplot(reepl1, aes(x=variable, y=value, fill=variable)) + geom_violin() + 
  stat_summary(fun.data=mean_sdl, mult=1,  geom="pointrange", color="DarkRed") + 
  geom_abline(intercept=0, slope=0) + geom_abline(intercept=10, slope=0, linetype="dashed") + geom_abline(intercept=-10, slope=0, linetype="dashed") +
  geom_abline(intercept=5, slope=0, linetype="dotted") + geom_abline(intercept=-5, slope=0, linetype="dotted")  +
  xlab("Parameter") + ylab("Relative estimation error (%)") + scale_fill_brewer(palette="Blues") + 
  theme_bw() + 
  #theme(legend.position = "none", base_size=14) + 
   theme(legend.position = "none") +   facet_grid(.~scenario)

# plot1 <- plot1+theme(text = element_text(size=rel(4)), legend.text=element_text(size=rel(4)), strip.text.x = element_text(size=rel(3.5)), strip.text.y = element_text(size=rel(3.5)))

# Results - Figure 1 (violin plot of the estimates for the parameters in the two scenarios)
if(saveFigure) {
  namfig<-"binarySimulation_REE.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot1)
  dev.off()
}

# Supplementary Figure 1
# Can't figure out how to add mean per group :-/
plot2 <- ggplot(reepl, aes(x=variable, y=value, fill=setting)) + geom_violin() + 
#  stat_summary(fun.data=mean_sdl, mult=1,  geom="pointrange", color="DarkRed") + 
  geom_abline(intercept=0, slope=0) + geom_abline(intercept=10, slope=0, linetype="dashed") + geom_abline(intercept=-10, slope=0, linetype="dashed") +
  geom_abline(intercept=5, slope=0, linetype="dotted") + geom_abline(intercept=-5, slope=0, linetype="dotted")  +
  xlab("Parameter") + ylab("Relative estimation error (%)") + scale_fill_brewer(palette="Blues") +  theme_bw() +    
  theme(legend.position = "top") +  facet_grid(.~scenario)

# plot2 <- plot2+theme(text = element_text(size=rel(4)), legend.text=element_text(size=rel(4)), strip.text.x = element_text(size=rel(3.5)), strip.text.y = element_text(size=rel(3.5)))

if(saveFigure) {
  namfig<-"binarySimulation_settingsREE.eps"
  cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
  plot(plot2)
  dev.off()
}

# Summary function
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
# p + stat_summary(fun.data=data_summary)

###################################################### Table with the bias and SE with 3 settings each time !
#colnames(reese)<-c("Original model","IIV on both parameters")
rownames(reese)<-c("$\\theta_1$","$\\theta_2$","$\\beta$","$\\omega_1$","$\\omega_2$")

print(xtable(reese), only.contents=TRUE, include.rownames=T, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity)

###################################################### SE
