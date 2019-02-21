progDir<-file.path(saemixDir,"R")
datDir<-file.path(saemixDir,"data")

{
  source(file.path(progDir,"aaa_generics.R"))
  #source(file.path(progDir,"global.R"))
  source(file.path(progDir,"SaemixData.R"))
  source(file.path(progDir,"SaemixRes.R"))
  source(file.path(progDir,"SaemixModel.R"))
  source(file.path(progDir,"SaemixObject.R"))
  source(file.path(progDir,"main.R"))
  source(file.path(progDir,"func_aux.R"))
  source(file.path(progDir,"main_initialiseMainAlgo.R"))
  source(file.path(progDir,"main_estep.R"))
  source(file.path(progDir,"main_mstep.R"))
  source(file.path(progDir,"func_FIM.R"))
  source(file.path(progDir,"func_plots.R"))
  source(file.path(progDir,"func_distcond.R"))
  source(file.path(progDir,"func_simulations.R"))
  source(file.path(progDir,"compute_LL.R"))
  source(file.path(progDir,"func_estimParam.R"))
}

# Files not yet transferred and validated
if(FALSE) {
  newCodeDir<-file.path(saemixDir,"newCode")
  source(file.path(newCodeDir,"SaemixObject-predict2.R"))
  source(file.path(newCodeDir,"estimateParameters.R"))
}


