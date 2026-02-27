cat("This code needs to be run in the root folder (folder paperSaemix3 on github), all paths are relative\n")

# workDir<-"/home/eco/work/saemix/saemixextension/paperSaemix3"
workDir <- getwd()

source(file.path(workDir, "saemix3_categoricalModel.R"))

source(file.path(workDir, "saemix3_countModel.R"))

source(file.path(workDir, "saemix3_tteModel.R"))
