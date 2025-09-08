# Wiping the slate clean and positioning in the right directory
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

#####################################################
# Reporter (location=each individual expect_xx, default=progress report)
myreporter <- default_reporter()
# myreporter <- "location"


#####################################################
###### Testing Data Class
# New functions (Rext)
source(file.path(progDirExt,"SaemixData.R"))
source(file.path(progDirExt,"SaemixData-methods.R"))
source(file.path(progDirExt,"SaemixData-methods_covariates.R"))

#############################
# Legacy tests for Data
# expect many warnings (NA due to conversion)
#### saemixData class - legacy code
## caution: need to copy tests to a different directory otherwise old functions are reloaded :-/

ifile <- "testthat_saemixData-class.R"
file.copy(from=file.path(testDir,ifile), to=testDirLegacy, overwrite=TRUE)
test_file(file.path(testDirLegacy,ifile), reporter=myreporter)

#### creating an SaemixData object from files or dataframes
#### automatic recognition of column names

## ToDo: check FAIL
ifile <- "testthat_saemixData-read.R"
file.copy(from=file.path(testDir,ifile), to=testDirLegacy, overwrite=TRUE)
test_file(file.path(testDirLegacy,ifile), reporter=myreporter)

#### datasets with more complex covariates
## ToDo: check FAIL
ifile <- "testthat_saemixData-covariates.R"
file.copy(from=file.path(testDir,ifile), to=testDirLegacy, overwrite=TRUE)
test_file(file.path(testDirLegacy,ifile), reporter=myreporter)

#### transforming continuous and categorical covariates
ifile <- "testthat_saemixData-transform.R"
file.copy(from=file.path(testDir,ifile), to=testDirLegacy, overwrite=TRUE)
test_file(file.path(testDirLegacy,ifile), reporter=myreporter)

#### plots
ifile <- "testthat_saemixData-plot.R"
file.copy(from=file.path(testDir,ifile), to=testDirLegacy, overwrite=TRUE)
test_file(file.path(testDirLegacy,ifile), reporter=myreporter)

#############################
# New tests for Data

#### class created with outcome (character vector in the data class)
test_file(file.path(testDirExt,"testthat_saemixData.R"), reporter=myreporter)

# ToDo:
## add more tests for the class with IOV
## add tests with covariates read and associated to different varlevels
## debug the last 3 tests

#####################################################
###### Testing new Classes for model

#############################
# Outcome
source(file.path(progDirExt,"SaemixErrorModel.R"))
source(file.path(progDirExt,"SaemixOutcome.R"))

# 
test_file(file.path(testDirExt,"testthat_saemixErrorModel.R"), reporter=myreporter)
test_file(file.path(testDirExt,"testthat_saemixOutcome-class.R"), reporter=myreporter)

#############################
# Parameter
source(file.path(progDirExt,"SaemixParameter.R"))
source(file.path(progDirExt,"SaemixParameter-methods.R"))

# ToDo
## Add parameters with covariates on different variability levels
# ToDo: add tests for list matching

test_file(file.path(testDirExt,"testthat_saemixParameter-class.R"), reporter=myreporter)

#############################
# VarLevel
source(file.path(progDirExt,"SaemixVarModel.R"))
source(file.path(progDirExt,"SaemixVarModel-methods.R"))

# ToDo
## add more test in particular with invalid correlaton structures
## add tests for indices
test_file(file.path(testDirExt,"testthat_saemixVarModel-class.R"), reporter=myreporter)

#############################
# Model for fixed effects ("population individual model", ie mu+sum_cov beta.cov without eta)
source(file.path(progDirExt,"SaemixPopModel.R"))

test_file(file.path(testDirExt,"testthat_saemixPopModel-class.R"), reporter=myreporter)

#####################################################
###### Testing Individual Model Class
source(file.path(progDirExt,"SaemixIndivModel.R"))

test_file(file.path(testDirExt,"testthat_saemixIndivModel-class.R"), reporter=myreporter)

#####################################################
###### Testing Model Class
source(file.path(progDirExt,"SaemixModel.R"))
source(file.path(progDirExt,"SaemixModel-methods.R"))

#############################
test_file(file.path(testDirExt,"testthat_saemixModel-class.R"))

test_file(file.path(testDirExt,"testthat_auxFunc.R"))

#####################################################
