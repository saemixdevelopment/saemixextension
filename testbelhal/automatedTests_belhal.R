# Wiping the slate clean and positioning in the right directory
rm(list = ls())
library(testthat)

# Testing where we are :-)
saemixDir<-"/home/eco/work/monolix/rversion/newLib/saemix"
if(is.na(file.info(saemixDir)[1])) {
  saemixDir<-"/Users/karimimohammedbelhal/Desktop/Phd/R_Package/contributions/FinalsaemixExtension/ecomets/saemix"
  rootDir<-"/Users/karimimohammedbelhal/Desktop/Phd/R_Package/contributions/FinalsaemixExtension/ecomets"
  datDir<-"/Users/karimimohammedbelhal/Desktop/Phd/R_Package/contributions/FinalsaemixExtension/ecomets/data"
}

setwd(saemixDir)


#####################################################
# Loading library files and defining data directory
datDir<-file.path(saemixDir,"saemix","data")
source(file.path(saemixDir,"testbelhal","helper-source.R"))

#####################################################
###### Testing classes for TTE and ordinal model

# Data - expect many warnings (NA due to conversion)
# Model

# => Belhal

###########################################################
################### Running models with non-continuous data
###########################################################
# Running TTE model
source(file.path(saemixDir,"testbelhal","test_setup_tte.R"))

# => Belhal
test_file(file.path(saemixDir,"testbelhal","testthat_predict_tte.R"))

#####################################################
# Running ordinal data model
source(file.path(saemixDir,"testbelhal","test_setup_ord.R"))

# => Belhal
test_file(file.path(saemixDir,"testbelhal","testthat_predict_ord.R"))

#####################################################
# Running count data model

## COUNT DATA => TODO

#####################################################
