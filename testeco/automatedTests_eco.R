# Wiping the slate clean and positioning in the right directory
rm(list = ls())
library(testthat)
library(rlang) # is_missing

# Testing where we are :-)
saemixDir<-"/home/eco/work/saemix/saemixextension"
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

test_file(file.path(saemixDir,"testeco","testthat_saemixData.R"))

test_file(file.path(saemixDir,"testeco","testthat_saemixData-plot.R"))

# Problème d'environnement - obligée de définir les data frame passés à saemixData dans l'environnement global et de nettoyer après, pas très clean !!!
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

# Model
test_file(file.path(saemixDir,"testeco","testthat_saemixModel-class.R"))

test_file(file.path(saemixDir,"testeco","testthat_saemixModel-function.R"))

#####################################################
###### Testing auxiliary functions

test_file(file.path(saemixDir,"testeco","testthat_ssq_combined2.R"))

###########################################################
################### Running models with continuous data
###########################################################
source(file.path(saemixDir,"testeco","test_setup_cont.R"))

# Testing saemix object with data and model
test_file(file.path(saemixDir,"testeco","testthat_saemixObject.R"))

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
# Estimating individual parameters after initialising to pop parameters from a fit




#####################################################
