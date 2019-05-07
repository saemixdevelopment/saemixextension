##########################################################################
# Directories
rootDir<-"/home/eco/work/monolix/rversion/newLib"
setwd(rootDir)
simDir<-file.path(rootDir,"zesims")

### Librairies and functions - saemix
saemixDir<-"/home/eco/work/monolix/rversion/newLib/saemix"
setwd(saemixDir)
source(file.path(saemixDir,"testeco","helper-source.R"))

##########################################################################
iscenar<-1

simulate.data<-TRUE # Simulate original data
source(file.path(rcode,"simulate_data.R"))


##########################################################################
