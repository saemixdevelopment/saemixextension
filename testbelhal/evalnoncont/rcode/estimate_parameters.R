##########################################################################
### Librairies and functions - saemix
rootDir<-"/home/eco/work/monolix/rversion/newLib"
simDir<-file.path(rootDir,"zesims")
setwd(rootDir)
saemixDir<-"/home/eco/work/monolix/rversion/newLib/saemix"
source(file.path(rootDir,"saemix","testeco","helper-source.R"))

# Directories

# library(saemix)
library(MASS)
library(ggplot2)

### Librairies and functions - mlx
if(FALSE) {
  source(file.path(condbootDir,"mlxConnectors","newConnectors.R"))
  source(file.path(condbootDir,"rfiles","functions_mlx.R"))
  source(file.path(condbootDir,"rfiles","functions.R"))
  
  library(MlxConnectors)
  initializeMlxConnectors(software = "monolix")
}

##########################################################################
estim.parOrig<-FALSE # estimate parameters for original simulated datafiles

# PD model, rich design
iscenar<-1
nsim<-200
source(file.path(simDir,"rcode","select_model.R"))

cat("Estimating parameters for example",iscenar,"\n")
source(file.path(simDir,"rcode","create_saemixModel.R"))

source(file.path(simDir,"rcode","run_saemixModel.R"))


##########################################################################
