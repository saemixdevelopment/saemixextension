##########################################################################
# Directories
rootDir<-"/Users/karimimohammedbelhal/Desktop/R_package/oldsaemixextension/testbelhal"
setwd(rootDir)
simDir<-file.path(rootDir,"zesims")

### Librairies and functions - saemix
saemixDir<-"/Users/karimimohammedbelhal/Desktop/R_package/oldsaemixextension"
setwd(saemixDir)
source(file.path(saemixDir,"testbelhal","helper-source.R"))

dir.create(file.path(simDir,"simData"))
dir.create(file.path(simDir,"simEstim"))

##########################################################################
# Simulation examples
#####################
simulate.data<-TRUE # Simulate original data

# Logistic model
nsim<-20
for(iscenar in 1:3) {
  cat("Simulating data for example",iscenar,"\n")
  source(file.path(simDir,"rcode","select_model.R"))
  dir.create(datDir)
  dir.create(resDir)
  # Répertoire de sauvegarde
  #dir.create(saveDir)
  source(file.path(simDir,"rcode","simulate_data.R"))
}

##########################################################################