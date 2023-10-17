


#---Sourcing saemix functions
saemixDir<- file.path(rootDir, "R") 
{
  source(file.path(saemixDir,"aaa_generics.R"))
  #source(file.path(saemixDir,"global.R"))
  source(file.path(saemixDir,"SaemixData.R"))
  source(file.path(saemixDir,"SaemixRes.R"))
  source(file.path(saemixDir,"SaemixModel.R"))
  source(file.path(saemixDir,"SaemixObject.R"))
  source(file.path(saemixDir,"main.R"))
  source(file.path(saemixDir,"func_aux.R"))
  source(file.path(saemixDir,"main_initialiseMainAlgo.R"))
  source(file.path(saemixDir,"main_estep.R"))
  source(file.path(saemixDir,"main_mstep.R"))
  source(file.path(saemixDir,"func_FIM.R"))
  source(file.path(saemixDir,"func_npde.R"))
  source(file.path(saemixDir,"func_plots.R"))
  source(file.path(saemixDir,"func_distcond.R"))
  source(file.path(saemixDir,"func_simulations.R"))
  source(file.path(saemixDir,"compute_LL.R"))
  source(file.path(saemixDir,"func_estimParam.R"))
  source(file.path(saemixDir,"func_bootstrap.R"))
  
  source(file.path(saemixDir,"backward.R"))
  source(file.path(saemixDir,"forward.R"))
  source(file.path(saemixDir,"stepwise.R"))
  source(file.path(saemixDir,"func_stepwise.R"))
  source(file.path(saemixDir,"func_compare.R"))
}


# #---Sourcing joint functions
# jointDir <- file.path(rootDir, "joint")
# {
#   source(file.path(jointDir,"multi_aux.R"))
#   source(file.path(jointDir,"multi_initializeMainAlgo.R"))
#   source(file.path(jointDir,"multi_estep.R"))
#   source(file.path(jointDir,"multi_mstep.R"))
#   source(file.path(jointDir,"multi_main.R"))
#   source(file.path(jointDir,"multi_map.R"))
#   source(file.path(jointDir,"compute_LL_multi.R"))
# }
# 



#---R files 
sourceRDir <- file.path(rootDir, "FIM_stochastique") 
{
  source(file.path(sourceRDir,"stocha_aux.R"))
  source(file.path(sourceRDir,"stocha_mstep.R"))
  source(file.path(sourceRDir,"stocha_main.R"))
  source(file.path(sourceRDir,"stocha_SaemixRes.R"))
}

library(LaplacesDemon)  # as.symmetric.matrix as.positive.definite

varCovMat <- function(mat){
  # Check if mat is symetric and make it if not 
  # Check if mat is positive definite definite and make it if not
  
  if(!(is.symmetric.matrix(mat))){
    mat = as.symmetric.matrix(mat)
  }
  if(!(is.positive.definite(mat))){
    mat = as.positive.definite(mat)
  }
  
  return(mat)
}