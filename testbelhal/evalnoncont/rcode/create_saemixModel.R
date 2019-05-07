
# Monolix model
if(iscenar %in% c(1,4,7)) {
  saemix.model<-saemixModel(model=modfun, description="PD Emax/Hill model",
                            psi0=matrix(parpop, ncol=3, byrow = TRUE, dimnames=list(NULL, nampar)),transform.par=c(1,1,1),
                            covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),
                            omega.init = omega, error.init=c(0,respar), error.model="proportional")
  saemix.options<-list(fix.seed=FALSE,directory=file.path(simDir,"simEstim","current"))
}
if(iscenar %in% c(2:3,5:6,8:9)) {
  omega2<-diag(rep(1,4))
  omega2[1:3,1:3]<-omega
  saemix.model<-saemixModel(model=modfun, description="PD Emax/Hill model",
                            psi0=matrix(parpop, ncol=4, byrow = TRUE, dimnames=list(NULL, nampar)),transform.par=c(1,1,1,1),
                            covariance.model=matrix(c(1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0),ncol=4,byrow=TRUE),
                            omega.init = omega2, error.init=c(0,respar), error.model="proportional")
  saemix.options<-list(fix.seed=FALSE,directory=file.path(simDir,"simEstim","current"))
}
