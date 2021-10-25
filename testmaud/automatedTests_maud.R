
rm(list = ls())
library(testthat)

# Loading library files and defining data directory

saemixDir<-"/Users/mdelattre/Documents/Boulot/Recherche/Travaux packages R/NewSaemix/saemixextension"

progDir<-file.path(saemixDir,"R")
datDir<-file.path(saemixDir,"data")
datDir2<-file.path(saemixDir,"testbelhal")

setwd(progDir)
repo.files <- list.files()

for (j in 1:length(repo.files)){
  file <- repo.files[j]
  source(file)
}


####################################################
################### Testing bic and related methods 
####################################################

source(file.path(saemixDir,"testmaud","test_setup_cont.R"))

# Check for correct computation of bic.covariate quantities

test_file(file.path(saemixDir,"testmaud","testthat_bic-covariate.R"))

# Check compare.saemix method

test_file(file.path(saemixDir,"testmaud","testthat_compare.R"))

# Check step.saemix method

test_file(file.path(saemixDir,"testmaud","testthat_stepwise.R"))
