##########################################################################
### Librairies and functions - saemix
rootDir<-"/Users/karimimohammedbelhal/Desktop/R_package/oldsaemixextension/testbelhal"
simDir<-file.path(rootDir,"zesims")
setwd(rootDir)
saemixDir<-"/Users/karimimohammedbelhal/Desktop/R_package/oldsaemixextension"
source(file.path(rootDir,"saemix","testbelhal","helper-source.R"))

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

# Binary model
iscenar<-1
nsim<-20
source(file.path(simDir,"rcode","select_model.R"))

cat("Estimating parameters for example",iscenar,"\n")
source(file.path(simDir,"rcode","create_saemixModel.R"))

source(file.path(simDir,"rcode","run_saemixModel.R"))


##########################################################################
