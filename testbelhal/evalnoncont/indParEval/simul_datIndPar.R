##########################################################################
# Directories
rootDir<-"/Users/karimimohammedbelhal/Desktop/R_package/oldsaemixextension/testbelhal"
setwd(rootDir)
simDir<-file.path(rootDir,"zesims")

### Librairies and functions - saemix
saemixDir<-"/Users/karimimohammedbelhal/Desktop/R_package/oldsaemixextension"
setwd(saemixDir)
source(file.path(saemixDir,"testbelhal","helper-source.R"))

##########################################################################
iscenar<-1

simulate.data<-TRUE # Simulate original data
source(file.path(rcode,"simulate_data.R"))


##########################################################################
