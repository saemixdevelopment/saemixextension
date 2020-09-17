#####################################################
# Testing where we are :-)
saemixDir<-"/home/eco/work/monolix/rversion/newLib/saemix"
setwd(saemixDir)

# Loading library files and defining data directory
source(file.path(saemixDir,"testeco","helper-source.R"))
datDir<-file.path(saemixDir,"data")

# Code bootstrap
source(file.path(saemixDir,"bootstrap","caseBootstrap.R"))

#####################################################
# Continuous model

source(file.path(saemixDir,"testeco","test_setup_cont.R"))
saemix.fit<-theo.fit1

# nboot<-100 # Nb of samples from the conditional distribution
nboot<-20 # Nb of samples from the conditional distribution
source(file.path(saemixDir,"bootstrap","bootstrapDistribution.R"))

# Bootstrap distribution of the parameters
par(mfrow=c(2,2))
for(i in 1:4) hist(res.boot[,(i+1)],main="",xlab=colnames(res.boot)[i+1])

par(mfrow=c(2,2))
for(i in 1:4) hist(res.boot[,(i+5)],main="",xlab=colnames(res.boot)[i+5])

#####################################################
# Longitudinal binary model
source(file.path(saemixDir,"testbelhal","test_setup_binomial.R"))
saemix.fit<-binary.fit

# nboot<-100 # Nb of samples from the conditional distribution
nboot<-50 # Nb of samples from the conditional distribution
source(file.path(saemixDir,"bootstrap","bootstrapDistribution.R"))

par(mfrow=c(2,2))
for(i in 1:4) hist(res.boot[,(i+1)],main="",xlab=colnames(res.boot)[i+1])

# TODO - simulation study with binary model

#####################################################
