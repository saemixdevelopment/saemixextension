# PK model
model1cpt<-function(psi,id,xidep) { 
  tim<-xidep[,1]
  dose<-xidep[,2]
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypk<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypk)
}
# PK/PD model
model1cptdirect<-function(psi,id,xidep) { 
  tim<-xidep[,1]
  dose<-xidep[,2]
  ytype<-xidep$ytype
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  ic50<-psi[id, 4]
  k<-CL/V
  ypk<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  ypd<-100*(1-ypk/(ypk+ic50))
  ypk[ytype==2]<-ypd[ytype==2]
  return(ypk)
}

chooseDebugModel <- function(model=1) {
  if(model %in% c(1:11)) {
    # Model outcomes
    out1<-list(conc=saemixOutcome(unit="mg/L", model="combined2", start=c(1, 0.2)))
    out2<-list(conc=saemixOutcome(unit="mg/L", model="proportional", start=c(0.5)))
    pkpd.outcome<-list(conc=saemixOutcome(unit="mg/L", model="combined2", start=c(1, 0.2)), 
                       effect=saemixOutcome(unit="%", start=1))
    
    # Mean value for population parameters
    pk.psi0<-c(ka=1, vd=5, cl=0.1)
    pkpd.psi0<-c(ka=1, vd=5, cl=0.1, ic50=5)
    # Model parameters
    
    # PK models
    if(model==1) {
      lpar1 <- list(ka=saemixPar(mu.start=pk.psi0[1], omega.start=0.3),
                    cl=saemixPar(mu.start=pk.psi0[2], omega.start=0.5),
                    vd=saemixPar(mu.start=pk.psi0[3], omega.start=0.7))
      saemix.model<-saemixModel(model=model1cpt, outcome=out1, parameter=lpar1)
    }
    if(model==2) {
      lpar2 <- list(ka=saemixPar(mu.start=pk.psi0[1]*1.5, omega=0),
                    cl=saemixPar(mu.start=pk.psi0[2], omega.start=0.5),
                    vd=saemixPar(mu.start=pk.psi0[3], omega.start=0.7,rho.param=c("cl"), rho=0.5))
      saemix.model<-saemixModel(model=model1cpt, outcome=out2, parameter=lpar2)
    }
    if(model==3) {
      lpar3 <- list(ka=saemixPar(mu.start=pk.psi0[1], omega=0, fixed=TRUE),
                    cl=saemixPar(mu.start=pk.psi0[2], omega.start=0.5),
                    vd=saemixPar(mu.start=pk.psi0[3], omega.start=0.7,rho.param=c("cl"), rho=0.5))
      saemix.model<-saemixModel(model=model1cpt, outcome=out2, parameter=lpar3)
    }
    if(model==4) {
      lpar4 <- list(ka=saemixPar(mu.start=pk.psi0[1], omega.start=1, covariate=c(sex=binCov(beta=0.2))),
                    cl=saemixPar(mu.start=pk.psi0[2], omega.start=0.5, covariate=c(wt=contCov(beta=0.75, beta.fix=1), age=contCov(beta=-0.5))),
                    vd=saemixPar(mu.start=pk.psi0[3], omega.start=0.7,rho.param=c("cl"), rho=0.5, covariate=c(wt=contCov(beta=1, beta.fix=1), sex=binCov(beta=0.4))))
      saemix.model<-saemixModel(model=model1cpt, outcome=out1, parameter=lpar4)
    }
    if(model==5) {
      lpar5 <- list(ka=saemixPar(mu.start=pk.psi0[1]*1.5, fixed=TRUE,  omega.start=1),
                    cl=saemixPar(mu.start=pk.psi0[2], omega.start=0.5, covariate=c(sex=binCov(beta=0.2), wt=contCov(beta=0.75, beta.fix=1), age=contCov(beta=-0.5))),
                    vd=saemixPar(mu.start=pk.psi0[3], omega.start=0.7,rho.param=c("cl"), rho=0.5, covariate=c(wt=contCov(beta=1, beta.fix=1), age=contCov(beta=-0.15))))
      saemix.model<-saemixModel(model=model1cpt, outcome=out1, parameter=lpar5)
    }
    if(model==6) {
      lpar6 <- list(ka=saemixPar(mu.start=pkpd.psi0[1], fixed=TRUE, omega=0, covariate=c(sex=binCov(beta=0.2))),
                    cl=saemixPar(mu.start=pkpd.psi0[2], omega.start=0.5, covariate=c(wt=contCov(beta=0.75, beta.fix=1), age=contCov(beta=-0.5))),
                    vd=saemixPar(mu.start=pkpd.psi0[3], omega.start=0.7,rho.param=c("cl"), rho=0.5, covariate=c(wt=contCov(beta=1, beta.fix=1), sex=binCov(beta=0.15))))
      saemix.model<-saemixModel(model=model1cpt, outcome=out1, parameter=lpar6)
    }
    if(model==7) {
      lpar7 <- list(ka=saemixPar(mu.start=pkpd.psi0[1], fixed=TRUE, omega=0),
                    cl=saemixPar(mu.start=pkpd.psi0[2], omega.start=0.5, covariate=c(sex=binCov(beta=0.2), wt=contCov(beta=0.75), age=contCov(beta=-0.5))),
                    vd=saemixPar(mu.start=pkpd.psi0[3], omega.start=0.7,rho.param=c("cl"), rho=0.5, covariate=c(wt=contCov(beta=1), sex=binCov(beta=0.15))))
      saemix.model<-saemixModel(model=model1cpt, outcome=out1, parameter=lpar7)
      
    }
    if(model==8) {
      lpar8 <- list(ka=saemixPar(mu.start=pkpd.psi0[1]*1.5, omega=0, covariate=c(sex=binCov(beta=0.2), wt=contCov(beta=0.5, beta.fix=1))),
                    cl=saemixPar(mu.start=pkpd.psi0[2], omega.start=0.5, covariate=c(sex=binCov(beta=0.2), wt=contCov(beta=0.75), age=contCov(beta=-0.5))),
                    vd=saemixPar(mu.start=pkpd.psi0[3], omega.start=0.7,rho.param=c("cl"), rho=0.5, covariate=c(wt=contCov(beta=1), sex=binCov(beta=0.15))))
      saemix.model<-saemixModel(model=model1cpt, outcome=out1, parameter=lpar8)
    }
    if(model==9) {
      lpar9 <- list(ka=saemixPar(mu.start=pkpd.psi0[1], omega=0, covariate=c(sex=binCov(beta=0.2), age=contCov(beta=-0.5))),
                    cl=saemixPar(mu.start=pkpd.psi0[2], omega.start=0.5, covariate=c(wt=contCov(beta=0.75), age=contCov(beta=-0.5))),
                    vd=saemixPar(mu.start=pkpd.psi0[3], omega.start=0.7,rho.param=c("cl"), rho=0.5, covariate=c(wt=contCov(beta=1), sex=binCov(beta=0.15))))
      saemix.model<-saemixModel(model=model1cpt, outcome=out1, parameter=lpar9)
    }
    
    # PK/PD models
    if(model==10) {
      lpar <- list(ka=saemixPar(mu.start=pkpd.psi0[1], omega.start=0.3),
                   cl=saemixPar(mu.start=pkpd.psi0[2], omega.start=0.5, covariate=c(sex=binCov(beta=0.5))),
                   vd=saemixPar(mu.start=pkpd.psi0[3], omega.start=0.7, covariate=c(wt=contCov(beta=1, beta.fix=1)), rho.param=c("cl"), rho=0.5),
                   ic50=saemixPar(mu.start=pkpd.psi0[4]))
      saemix.model<-saemixModel(model=model1cptdirect, outcome=pkpd.outcome, parameter=lpar)
    }
    if(model==11) {
      lpar <- list(ka=saemixPar(mu.start=pkpd.psi0[1], omega.start=0.2, covariate=c(sex=binCov(beta=0.3))),
                   cl=saemixPar(mu.start=pkpd.psi0[2], omega.start=c(0.3,0.4), covariate=c(sex=binCov(beta=0.2), lwt=contCov(name="wt",transform=log, beta=0.75, beta.fix=1), lage=contCov(name="age", transform=log, beta=-0.5))),
                   vd=saemixPar(mu.start=pkpd.psi0[3], omega.start=c(0.3,0.4), covariate=c(lwt=contCov(name="wt",transform=log, beta=1, beta.fix=1), cyp=binCov(beta=0.5)), rho.param=c("cl"), rho=0.5),
                   ic50=saemixPar(mu.start=pkpd.psi0[4], omega.start=0.4, covariate=c(lage=contCov(name="age", transform=log, beta=0.5), trt=catCov(beta=c(0.1,0.3)))))
      saemix.model<-saemixModel(model=model1cptdirect, outcome=pkpd.outcome, parameter=lpar)
    }
    return(saemix.model)
  }
}
