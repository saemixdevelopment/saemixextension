# Wiping the slate clean and positioning in the right directory
rm(list = ls())
library(testthat)

# Testing where we are :-)
saemixDir<-"/home/eco/work/monolix/rversion/newLib/saemix"
if(is.na(file.info(saemixDir)[1])) {
  saemixDir<-"/Users/karimimohammedbelhal/Desktop/Phd/R_Package/contributions/FinalsaemixExtension/ecomets/saemix"
}

setwd(saemixDir)


#####################################################
# Loading library files and defining data directory
source(file.path(saemixDir,"testeco","helper-source.R"))
datDir<-file.path(saemixDir,"data")

#####################################################
###### Testing classes

# Data - expect many warnings (NA due to conversion)
test_file(file.path(saemixDir,"testeco","testthat_saemixData-class.R"))

test_file(file.path(saemixDir,"testeco","testthat_saemixData-read.R"))

# Problème d'environnement - obligée de définir les data frame passés à saemixData dans l'environnement global et de nettoyer après
rm(theo.saemix)
rm(PD1.saemix)
rm(PD2.saemix)
test_file(file.path(saemixDir,"testeco","testthat_saemixData-covariates.R"))
rm(theo.saemix)
rm(PD1.saemix)
rm(PD2.saemix)
rm(tab1)
rm(tab2)
rm(tab3)
rm(pkpddat)
rm(cow.saemix)

# Model and object
test_file(file.path(saemixDir,"testeco","testthat_saemixModel.R"))

test_file(file.path(saemixDir,"testeco","testthat_saemixObject.R"))

#####################################################
###### Testing auxiliary functions

test_file(file.path(saemixDir,"testeco","testthat_functions.R"))

###########################################################
################### Running models with continuous data
###########################################################
source(file.path(saemixDir,"testeco","test_setup_cont.R"))

################
# Testing methods
test_file(file.path(saemixDir,"testeco","testthat_summary.R"))

test_file(file.path(saemixDir,"testeco","testthat_replaceData_cont.R"))

# Testing predict functions
test_file(file.path(saemixDir,"testeco","testthat_predict.R"))

################
# Testing plots - no 'validation', just need user eye

# Testing LL, 3 methods (linearisation, IS, AGQ) 

# Testing FIM

#####################################################
