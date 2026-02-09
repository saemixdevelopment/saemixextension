# Execute all models and save figures
# Tables in LaTeX format are output but not saved 
library(xtable)

source("saemix3_categoricalModel.R")

###################################################### Main text (arXiv)
# Figure 3
# Diagnostics knee model with covariates

namfig<-"knee_VPCbytreatment.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(kneepl.VPCbytreatment)
dev.off()

###################################################### Supplementary (arXiv)

namfig<-"knee_rawDataPropTime.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot( kneepl.rawDataPropTime)
dev.off()

namfig<-"knee_VPCbySex.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(kneepl.VPCbySex)
dev.off()

namfig<-"knee_meanScoreVPC.eps"
cairo_ps(file = file.path(figDir, namfig), onefile = TRUE, fallback_resolution = 600, height=8.27, width=11.69)
plot(kneepl.VPCmeanScore)
dev.off()

###################################################### Tables
### Table 2 (main text)
#### combine base and covariate model estimates

nline<-dim(xtab.kneeCov)[1]+1
df <- as.data.frame(matrix(nrow=nline, ncol=4))
rownames(df)<-c(rownames(xtab.kneeCov),rownames(xtab.kneeBase)[7])
colnames(df)<-rep(c("Value","CI"),2)
df[c(1,3,5:7,9,13),1:2]<-xtab.kneeBase[,c(1,3)]
df[1:12,3:4]<-xtab.kneeCov[,c(1,3)]
for(icol in 1:4) {
  yvec<-df[,icol]
  df[is.na(yvec),icol]<-"-"
}

print(xtable(df), only.contents=TRUE, include.rownames=TRUE,  floating=F, sanitize.rownames.function = identity)

