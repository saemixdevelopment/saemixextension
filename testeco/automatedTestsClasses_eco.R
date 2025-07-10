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
###### Testing Data Class

# Data - expect many warnings (NA due to conversion)
test_file(file.path(saemixDir,"testeco","testthat_saemixData-class.R"))

test_file(file.path(saemixDir,"testeco","testthat_saemixData-read.R"))

test_file(file.path(saemixDir,"testeco","testthat_saemixData-covariates.R"))

test_file(file.path(saemixDir,"testeco","testthat_saemixData-transform.R"))


test_file(file.path(saemixDir,"testeco","testthat_saemixData-plot.R"))


# Problème d'environnement - obligée de définir les data frame passés à saemixData dans l'environnement global et de nettoyer après, pas très clean !!!
## solved ?

#####################################################
###### Testing Model Class

test_file(file.path(saemixDir,"testeco","testthat_saemixModel-class.R"))

test_file(file.path(saemixDir,"testeco","testthat_saemixModel-function.R"))

#####################################################
