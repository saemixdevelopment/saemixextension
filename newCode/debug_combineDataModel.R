#####################################################  Directories
saemixDir<-"/home/eco/work/saemix/saemixextension"
setwd(saemixDir)

# Defining data directory
progDir<-file.path(saemixDir,"R")
progDirExt<-file.path(saemixDir,"Rext")
datDir<-file.path(saemixDir,"data")
datDir40<-file.path(saemixDir,"data40")

##################################################### Data
# Loading library files - Data class
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDirExt,"SaemixOutcome.R"))
source(file.path(progDirExt,"SaemixData.R"))

##################################################### Model
# Loading library files - Covariate class
source(file.path(progDirExt,"SaemixCovariateModel.R"))
source(file.path(progDirExt,"SaemixCovariate.R"))

# Loading library files Variability levels
source(file.path(progDirExt,"SaemixVarLevel.R"))

# Loading library files Parameter
source(file.path(progDirExt,"SaemixParameter.R"))



##################################################### Combining Data and Model
# Loading library files - Model class
source(file.path(progDirExt,"SaemixModel.R"))
source(file.path(progDirExt,"SaemixModel-methods.R"))

# Create data object
# Create model object
# Check outcomes
# Check covariates
