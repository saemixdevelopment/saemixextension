# Loading libraries
library(xtable)
library(ggplot2)
library(tidyr)

# Loading saemix
#library(saemix)

# Folders
workDir<-"/home/eco/work/saemix/discreteEval"
# setwd(workDir)

saemixDir<-"/home/eco/work/saemix/saemixextension"
figDir <- file.path(workDir, "figs")
nsim<-200 # Number of simulations
saveFigures <- FALSE

###################################################### REE

reepl <- reese <- NULL
for(iscenar in 1:2) {
  if(iscenar==1) {
    # Original
    scenario <- "binaryOrig"
    datDir <- file.path(workDir, scenario, "data")
    resDir <- file.path(workDir, scenario, "results")
    par.fix<-c(-1.71, -0.39, -0.15)
    par.om <- c(4.02)
  }
  if(iscenar==2) {
    # IIV
    scenario <- "binaryIIV"
    datDir <- file.path(workDir, scenario, "data")
    resDir <- file.path(workDir, scenario, "results")
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
  reepl <- rbind(reepl,
                 cbind(ree %>%
                   gather(key="variable",value="value"),
                 setting="True", scenario=scenario))
  
  ree<-tab.pop[,c(2:(1+npar))]
  for(i in 1:npar) ree[,i]<-(ree[,i]-par.true[i])/par.true[i]
  x1<-apply(ree,2,mean)
  x2<-apply(ree,2,sd)
  l1<-paste0(format(x1*100, digits=2, nsmall=2)," (",format(x2*100, digits=1, nsmall=0),")")
  if(iscenar==1) l1 <-c(l1,c(""))
  reese <- cbind(reese,l1)
  ree <- ree*100
  
  reepl <- rbind(reepl,
                 cbind(ree %>%
                         gather(key="variable",value="value"), setting="Pop", scenario=scenario))
  ree<-tab.far[,c(2:(1+npar))]
  for(i in 1:npar) ree[,i]<-(ree[,i]-par.true[i])/par.true[i]
  x1<-apply(ree,2,mean)
  x2<-apply(ree,2,sd)
  l1<-paste0(format(x1*100, digits=2, nsmall=2)," (",format(x2*100, digits=1, nsmall=0),")")
  if(iscenar==1) l1 <-c(l1,c(""))
  reese <- cbind(reese,l1)
  ree <- ree*100
  
  reepl <- rbind(reepl,
                 cbind(ree %>%
                         gather(key="variable",value="value"), setting="Far", scenario=scenario))
}
reepl$scenario<-factor(reepl$scenario, levels=c("binaryOrig","binaryIIV"), labels=c("Original model","IIV on both parameters"))
reepl$variable <- factor(reepl$variable, levels=c("theta1","theta2", "beta","omega1","omega2"))

ggplot(reepl, aes(x=variable, y=value)) + geom_violin() + 
  geom_abline(intercept=0, slope=0) + geom_abline(intercept=10, slope=0, linetype="dashed") + geom_abline(intercept=-10, slope=0, linetype="dashed") +
  facet_grid(scenario~setting)

reepl1 <- reepl[reepl$setting=="True",]
plot1 <- ggplot(reepl1, aes(x=variable, y=value, fill=variable)) + geom_violin() + 
  stat_summary(fun.data=mean_sdl, mult=1,  geom="pointrange", color="DarkRed") + 
  geom_abline(intercept=0, slope=0) + geom_abline(intercept=10, slope=0, linetype="dashed") + geom_abline(intercept=-10, slope=0, linetype="dashed") +
  geom_abline(intercept=5, slope=0, linetype="dotted") + geom_abline(intercept=-5, slope=0, linetype="dotted")  +
  xlab("Parameter") + ylab("Relative estimation error (%)") + scale_fill_brewer(palette="Blues") + 
  theme_bw() + theme(legend.position = "none") + facet_grid(.~scenario)

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
  xlab("Parameter") + ylab("Relative estimation error (%)") + scale_fill_brewer(palette="Blues") + 
  theme_bw() + theme(legend.position = "top") + facet_grid(.~scenario)

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
