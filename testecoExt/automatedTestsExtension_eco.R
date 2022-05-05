##################################################### Setup
# Wiping the slate clean and positioning in the right directory
if(FALSE) 
  rm(list = ls())
library(testthat)

# Testing where we are :-)
saemixDir<-"/home/eco/work/saemix/saemixextension"
if(is.na(file.info(saemixDir)[1])) {
  saemixDir<-"/Users/karimimohammedbelhal/Desktop/Phd/R_Package/contributions/FinalsaemixExtension/ecomets/saemix"
}

setwd(saemixDir)

# Defining data directory
progDir<-file.path(saemixDir,"R")
progDirExt<-file.path(saemixDir,"Rext")
datDir<-file.path(saemixDir,"data")
datDir40<-file.path(saemixDir,"data40")

##################################################### Test Data class

# Loading library files - Data class
source(file.path(progDir,"aaa_generics.R"))
#source(file.path(progDir,"global.R"))
source(file.path(progDirExt,"SaemixOutcome.R"))
source(file.path(progDirExt,"SaemixOutcome-methods.R"))
source(file.path(progDirExt,"SaemixData.R"))

# Current status: passes
test_file(file.path(saemixDir,"testecoExt","testthat_saemixOutcome-class.R"))
# Current status: passes
test_file(file.path(saemixDir,"testecoExt","testthat_saemixOutcome-methods.R"))

# One error which seems appriopriately caught
# Error in initialize(value, ...) : 
# 'initialize' method returned an object of class “character” instead of the required class “SaemixDiscreteOutcome”

# Current status: fails (work in progress, need to define constructor functions and the lot)
test_file(file.path(saemixDir,"testecoExt","testthat_saemixData-class.R"))

##################################################### Test Model class
# Loading library files - Covariate class
source(file.path(progDirExt,"SaemixCovariateModel.R"))
source(file.path(progDirExt,"SaemixCovariate.R"))

# Current status: passes
test_file(file.path(saemixDir,"testecoExt","testthat_saemixCovariateModel-class.R"))
test_file(file.path(saemixDir,"testecoExt","testthat_saemixCovariate-class.R"))

################################# OK
# Loading library files Variability levels
source(file.path(progDirExt,"SaemixVarLevel.R"))

# Current status: passes
test_file(file.path(saemixDir,"testecoExt","testthat_saemixVarLevel-class.R"))

################################# OK
# Loading library files Parameter
source(file.path(progDirExt,"SaemixParameter.R"))
source(file.path(progDirExt,"SaemixParameter-methods.R"))

# Current status: passes
test_file(file.path(saemixDir,"testecoExt","testthat_saemixParameter-class.R"))
# currently not created, all tests in the first testthat file
# test_file(file.path(saemixDir,"testecoExt","testthat_saemixParameter-methods.R"))

################################# In progress
# Loading library files - Model class
source(file.path(progDirExt,"SaemixModel.R"))
source(file.path(progDirExt,"SaemixModel-methods.R"))

# Current status: passes
test_file(file.path(saemixDir,"testecoExt","testthat_saemixModel-class.R"))

# Current status: fails (work in progress, legacy versus new model specs)
test_file(file.path(saemixDir,"testecoExt","testthat_saemixModel-methods.R"))

##################################################### Test Model + Data
# Loading library files - Covariate transformation methods
source(file.path(progDirExt,"SaemixCovariateTransform.R"))
# Current status: fails (work in progress)
test_file(file.path(saemixDir,"testecoExt","testthat_saemixCovariateTransform-class.R"))

##################################################### Test Res class
# Loading library files - Model class
source(file.path(progDir,"SaemixRes.R"))
source(file.path(progDir,"SaemixObject.R"))



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
