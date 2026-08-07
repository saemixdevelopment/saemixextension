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


###############################################################################################################
###### Testing Data Class
# New functions (Rext)
source(file.path(progDirExt,"SaemixData.R"))
source(file.path(progDirExt,"SaemixData-printmethods.R"))
source(file.path(progDirExt,"SaemixData-methods.R"))
source(file.path(progDirExt,"SaemixData-covariatemethods.R"))
source(file.path(progDirExt,"SaemixData-plotmethods.R"))

#### creating an SaemixData object from files or dataframes
#### automatic recognition of column names
ifile <- "testthat_saemixData-read.R"
test_file(file.path(testDirExt,ifile), reporter=myreporter)

#### datasets with more complex covariates
ifile <- "testthat_saemixData-covariates.R"
test_file(file.path(testDirExt,ifile), reporter=myreporter)

#############################
# Legacy tests for Data (copied from previous classes, should still work/pass)
# expect many warnings (NA due to conversion)
#### saemixData class - legacy code
## caution: need to copy tests to a different directory otherwise old functions are reloaded :-/

ifile <- "testthat_saemixData-class.R"
file.copy(from=file.path(testDir,ifile), to=testDirLegacy, overwrite=TRUE)
test_file(file.path(testDirLegacy,ifile), reporter=myreporter)

#### transforming continuous and categorical covariates
ifile <- "testthat_saemixData-transform.R"
file.copy(from=file.path(testDir,ifile), to=testDirLegacy, overwrite=TRUE)
test_file(file.path(testDirLegacy,ifile), reporter=myreporter)

#### plots - ToDo
ifile <- "testthat_saemixData-plot.R"
file.copy(from=file.path(testDir,ifile), to=testDirLegacy, overwrite=TRUE)
test_file(file.path(testDirLegacy,ifile), reporter=myreporter)

##########################################################
# New tests for Data

#### class created with outcome (character vector in the data class)
test_file(file.path(testDirExt,"testthat_saemixData.R"), reporter=myreporter)

# ToDo:
## add more tests for the class with IOV
## add tests with covariates read and associated to different varlevels
## debug the last 3 tests
