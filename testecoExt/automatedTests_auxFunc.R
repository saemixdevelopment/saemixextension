# Context: check Data and Model classes first
if(FALSE) {
  source("automatedTests_classData.R")
  source("automatedTests_classModel.R")
}

#####################################################
rm(list = ls())
library(testthat)
library(rlang) # is_missing

# Testing where we are :-)
saemixDir<-"/home/eco/work/saemix/saemixextension"
#if(is.na(file.info(saemixDir)[1])) {
#  saemixDir<-"/Users/karimimohammedbelhal/Desktop/Phd/R_Package/contributions/FinalsaemixExtension/ecomets/saemix"
#}

setwd(saemixDir)

#####################################################
# Defining folders with code and data
progDir<-file.path(saemixDir,"R")
progDirExt<-file.path(saemixDir,"Rext")
datDir<-file.path(saemixDir,"data")
datDir40<-file.path(saemixDir,"data40")
testDirExt <- file.path(saemixDir,"testecoExt")
testDir <- file.path(saemixDir,"testeco")
testDirLegacy <- file.path(testDirExt,"legacyTests")

# Loading generic definitions
source(file.path(progDir,"aaa_generics.R"))

# Reporter (location=each individual expect_xx, default=progress report)
myreporter <- default_reporter()
# myreporter <- "location"

###############################################################################################################
###### Sourcing SaemixData class and methods
##############################################
source(file.path(progDirExt,"SaemixData.R"))
source(file.path(progDirExt,"SaemixData-printmethods.R"))
source(file.path(progDirExt,"SaemixData-methods.R"))
source(file.path(progDirExt,"SaemixData-covariatemethods.R"))
source(file.path(progDirExt,"SaemixData-plotmethods.R"))

###############################################################################################################
###### Sourcing Classes for SaemixModel
##############################################
# Outcome
source(file.path(progDirExt,"SaemixOutcome.R"))
source(file.path(progDirExt,"SaemixErrorModel.R"))

# Parameter
source(file.path(progDirExt,"SaemixParameter.R"))
source(file.path(progDirExt,"SaemixParameter-methods.R"))

# VarLevel
source(file.path(progDirExt,"SaemixVarModel.R"))
source(file.path(progDirExt,"SaemixVarModel-methods.R"))

# Model for fixed effects ("population individual model", ie mu+sum_cov beta.cov without eta)
source(file.path(progDirExt,"SaemixPopModel.R"))

# Individual Model Class
source(file.path(progDirExt,"SaemixIndivModel.R"))

# Model Class
source(file.path(progDirExt,"SaemixModel.R"))
source(file.path(progDirExt,"SaemixModel-methods.R"))

# Individual parameter Class
source(file.path(progDirExt,"SaemixVarPhi.R"))
source(file.path(progDirExt,"SaemixIndivPar.R"))

###############################################################################################################
###### Testing computational functions
##############################################

### WiP !!!

###### Testing computational functions
source(file.path(progDirExt,"func_initialise.R"))
source(file.path(progDirExt,"func_aux.R"))

test_file(file.path(testDirExt,"testthat_auxFunc.R"))

###### Testing initialisation

