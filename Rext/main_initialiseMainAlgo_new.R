# Filling the elements of SaemixIndividualModel as a list
# Use extractDesignMatrix ?

############################### Initialising main algorithm #############################
initialiseMainAlgo.new<-function(saemix.data,saemix.model,saemix.options) {
  # Function to reformat covariance structure and initialise lists used in the main algorithm
  # Input: data, model and options
  # Output:
  ### saemix.model: added elements betaest, covariate model, indices (); modified/formatted psi0 adjusting to the nb of covariates
  ### Dargs: data elements - passed on to functions, unchanged
  ### Uargs: list of indices and variables (fixed) - passed on to functions, unchanged
  ### varList: variability-related elements - passed to functions and optimised
  ### opt: list of options and settings (fixed) - passed on to functions, unchanged
  ### DYF: used for the acceptance/rejection algorithm
  ### parameters optimised during the fit: phiM, mean.phi, betas, fixedpsi.ini
  ### allpar0: array holding the successive values of population parameters

  ivarlev <- 1
  # Individual variability model(s) 
  ## TODO extend for more than 1 level
  modelStructure <- extractModelStructure(saemix.data, saemix.model, level=ivarlev, verbose=saemix.options$warnings)
  Uargs <- modelStructure$Uargs
  saemix.model@ind.model <- list(modelStructure$indiv.model)

  # Starting values for individual parameters
  mean.phi<- saemix.model@ind.model[[ivarlev]]@COV %*% saemix.model@ind.model[[ivarlev]]@MCOV
  initialPhi <- initialisePhi(saemix.data, saemix.model, mean.phi=mean.phi, nb.chains=saemix.options$nb.chains, verbose=saemix.options$warnings)
  saemix.model <- initialPhi$saemix.model # adjust omega (if non conform or some parameters have no IIV)
  phiM <- initialPhi$phiM
  Dargs <- initialPhi$Dargs
  DYF <- initialPhi$DYF
  
}