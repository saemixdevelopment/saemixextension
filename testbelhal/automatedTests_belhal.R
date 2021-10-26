# Wiping the slate clean and positioning in the right directory
rm(list = ls())
library(testthat)

# Testing where we are :-)
saemixDir<-"/home/eco/work/saemix/saemixextension/"
if(is.na(file.info(saemixDir)[1])) {
  saemixDir<-"/Users/karimimohammedbelhal/Desktop/ML_Research/Saemix - R_Packages/saemix/saemixextension/"
  rootDir<-"/Users/karimimohammedbelhal/Desktop/ML_Research/Saemix - R_Packages/saemix/saemixextension/testbelhal"
  datDir<-"/Users/karimimohammedbelhal/Desktop/ML_Research/Saemix - R_Packages/saemix/saemixextension/data"
}

setwd(saemixDir)


#####################################################
# Loading library files and defining data directory
datDir<-file.path(saemixDir,"testbelhal")
source(file.path(saemixDir,"testbelhal","helper-source.R"))

#####################################################
###### Testing classes for TTE and ordinal model

# Data - expect many warnings (NA due to conversion)
test_file(file.path(saemixDir,"testbelhal","testthat_saemixData-class.R"))
test_file(file.path(saemixDir,"testbelhal","testthat_saemixData-read.R"))

# Model and object
test_file(file.path(saemixDir,"testbelhal","testthat_saemixModel.R"))
test_file(file.path(saemixDir,"testbelhal","testthat_saemixObject.R"))

#####################################################
###### Testing auxiliary functions
test_file(file.path(saemixDir,"testbelhal","testthat_functions.R"))



###########################################################
################### Running models with continuous data
###########################################################
# Running Model1cpt model
source(file.path(saemixDir,"testbelhal","test_setup_cont.R"))
test_file(file.path(saemixDir,"testbelhal","testthat_predict_cont.R"))


###########################################################
################### Running models with continuous data and new kernel
###########################################################
# Running Model1cpt model with 4th kernel
source(file.path(saemixDir,"testbelhal","test_setup_cont_kernel.R"))
test_file(file.path(saemixDir,"testbelhal","testthat_predict_cont_kernel.R"))


###########################################################
################### Running models with non-continuous data
###########################################################
# Running TTE model
source(file.path(saemixDir,"testbelhal","test_setup_tte.R"))
test_file(file.path(saemixDir,"testbelhal","testthat_predict_tte.R"))

#####################################################
# Running ordinal data model
source(file.path(saemixDir,"testbelhal","test_setup_ord.R"))
test_file(file.path(saemixDir,"testbelhal","testthat_predict_ord.R"))

#####################################################
# Running ordinal data model
source(file.path(saemixDir,"testbelhal","test_setup_binomial.R"))
# test_file(file.path(saemixDir,"testbelhal","testthat_predict_binomial.R"))

#####################################################
# Running count data model
source(file.path(saemixDir,"testbelhal","test_setup_count.R"))
### WHEN ONLY ONE PARAM TO ESTIMATE (fixed.estim=c(1,0)) OBTAIN: 
# Error in cbind(blocA, t(blocC)) : 
#   le nombre de lignes des matrices doit correspondre (voir argument 2)
test_file(file.path(saemixDir,"testbelhal","testthat_predict_count.R"))

#####################################################

