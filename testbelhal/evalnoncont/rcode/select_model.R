if(iscenar %in% c(1:9)) {
  source(file.path(simDir,"rcode","model_ordinal.R"))
}

# Scenario Binary logistic 1 - Emax model
if(iscenar %in% c(1:3)) {
  namsimdat<-"ordinal"
  parpop<-parpop[-c(4)]
  nampar<-nampar[-c(4)]
  modfun<-binary.model
}


namsimDir<-paste(namsimdat,design,sep=".")
datDir<-file.path(simDir,"simData",namsimDir)
resDir<-file.path(simDir,"simEstim",namsimDir)



