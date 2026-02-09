# Execute all models and save figures
# Tables in LaTeX format are output but not saved 
library(xtable)

source("saemix3_tteModel.R")

###################################################### Main text (arXiv)
# Figure 5 
# Lung survival data 

namfig<-"lung_exploreSurv.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(lungpl.exploreSurv)
dev.off()

# Figure 6
# 
namfig<-"lung_compareTTEfits.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(lungpl.compareTTEfits)
dev.off()

###################################################### Tables
# Table 5
# Comparing survival models on the lung cancer data

print(xtable(lungtab.comparemodels), only.contents=TRUE, include.rownames=FALSE,  floating=F, sanitize.rownames.function = identity)

# Table 6
# Parameter estimates for the model with covariates 
rownames(lungtab.parestim) <- c("T$\\_e$ (days)","$\\beta_{gender, T\\_e}$ (-)","$\\beta_{ECOG,T\\_e}$ (-)","$\\gamma$ (-)","$\\omega_{\\gamma}$ (-)")

xtab <- lungtab.parestim[,1:3]
xtab[,2]<-format(xtab[,2],digits=2, trim=T)
xtab[,3]<-paste0("[",format(lungtab.parestim[,4], digits=1,nsmall=1, trim=T),"-",format(lungtab.parestim[,5],digits=2, trim=T),"]")
colnames(xtab) <- c("Parameter","Estimate", "Bootstrap CI")
xtab<-xtab[,-c(1)]

print(xtable(xtab), only.contents=TRUE, include.rownames=TRUE,  floating=F, sanitize.rownames.function = identity)
