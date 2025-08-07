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
testDir <- file.path(saemixDir,"testeco")
testDirExt <- file.path(saemixDir,"testecoExt")

# Loading generic definitions
source(file.path(progDir,"aaa_generics.R"))

#####################################################
###### Testing Data Class
# New functions (Rext)
source(file.path(progDirExt,"SaemixData.R"))
source(file.path(progDirExt,"SaemixData-methods.R"))
# Covariate transformations (R)
source(file.path(progDirExt,"SaemixData-methods_covariates.R"))

#############################
# Legacy tests for Data
# expect many warnings (NA due to conversion)
#### saemixData class - legacy code
test_file(file.path(testDir,"testthat_saemixData-class.R"), reporter = "location")

#### creating an SaemixData object from files or dataframes
#### automatic recognition of column names
test_file(file.path(testDir,"testthat_saemixData-read.R"), reporter = "location")

#### datasets with more complex covariates
test_file(file.path(testDir,"testthat_saemixData-covariates.R"), reporter = "location")

#### transforming continuous and categorical covariates
test_file(file.path(testDir,"testthat_saemixData-transform.R"), reporter = "location")

#### plots
test_file(file.path(testDir,"testthat_saemixData-plot.R"), reporter = "location")

#############################
# New tests for Data

#### class created with outcome
test_file(file.path(testDirExt,"testthat_saemixData.R"), reporter="location")

# ToDo:
## add more tests for the class with IOV
## add tests with covariates read and associated to different varlevels
## debug the last 3 tests

#####################################################
###### Testing new Classes

#############################
# Outcome
source(file.path(progDirExt,"SaemixOutcome.R"))

# 
test_file(file.path(testDirExt,"testthat_saemixOutcome-class.R"), reporter="location")

#############################
# Parameter
source(file.path(progDirExt,"SaemixParameter.R"))
source(file.path(progDirExt,"SaemixParameter-methods.R"))

# ToDo
## Add parameters with covariates on different variability levels
# ToDo: add tests for list matching

test_file(file.path(testDirExt,"testthat_saemixParameter-class.R"), reporter="location")

#############################
# VarLevel
source(file.path(progDirExt,"SaemixVarModel.R"))
source(file.path(progDirExt,"SaemixVarModel-methods.R"))

# ToDo
## add more test in particular with invalid correlaton structures
## add tests for indices
test_file(file.path(testDirExt,"testthat_saemixVarModel-class.R"), reporter="location")

#############################
# Model for fixed effects ("population individual model", ie mu+sum_cov beta.cov without eta)
source(file.path(progDirExt,"SaemixPopModel.R"))

test_file(file.path(testDirExt,"testthat_saemixPopModel-class.R"), reporter="location")

#####################################################
###### Testing Individual Model Class
source(file.path(progDirExt,"SaemixIndivModel.R"))

test_file(file.path(testDirExt,"testthat_saemixIndivModel-class.R"), reporter="location")

#####################################################
###### Testing Model Class
source(file.path(progDirExt,"SaemixModel.R"))
source(file.path(progDirExt,"SaemixModel-methods.R"))

#############################
# Legacy tests for Model
test_file(file.path(testDir,"testthat_saemixModel-class.R"))

test_file(file.path(testDir,"testthat_saemixModel-function.R"))

#####################################################
