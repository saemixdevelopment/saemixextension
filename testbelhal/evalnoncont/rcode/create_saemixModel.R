
# Monolix model
if(iscenar %in% c(1:3)) {
saemix.model<-saemixModel(model=binary.model,description="Binary model",modeltype="likelihood",
        psi0=matrix(parpop,ncol=2,byrow=TRUE,dimnames=list(NULL,nampar)),
        transform.par=c(0,0),covariance.model=matrix(c(1,0,0,1),ncol=2, byrow=TRUE),omega.init = omega)

saemix.options<-list(fix.seed=FALSE,directory=file.path(simDir,"simEstim","current"))
}
