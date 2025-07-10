############### Theophylline example
library(saemix)
data(theo.saemix)

saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
                        name.group=c("Id"),name.predictors=c("Dose","Time"),
                        name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
                        units=list(x="hr",y="mg/L", covariates=c("kg","-")), name.X="Time", verbose=FALSE)
print(saemix.data)

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

saemix.model<-saemixModel(model=model1cpt,
                          description="One-compartment model with first-order absorption", 
                          psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
                                      dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
                          covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
                          covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
                          omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant", verbose=FALSE)
theo.fit<-saemix(model=saemix.model,
                   data=saemix.data,
                   control=list(seed=632545,directory="newtheo", save=FALSE,save.graphs=FALSE, warnings=FALSE))

############### 
saemixObject <- theo.fit


############### Initialising
#nnodes.gq<-saemixObject["options"]$nnodes.gq
nnodes.gq<-3 # sufficient in Ueckert et al. # number of nodes on each 1-D grid
nsd.gq<-saemixObject["options"]$nsd.gq   # the integral is computed on the interval [E(eta|y) +- nsd_gq*SD(eta|y)]

saemix.model<-saemixObject["model"]
saemix.data<-saemixObject["data"]
saemix.res<-saemixObject["results"]
xind<-saemix.data["data"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
yobs<-saemix.data["data"][,saemix.data["name.response"]]

#  covariance.model<-0*saemix.model["covariance.model"]
covariance.model<-saemix.model["covariance.model"]
omega<-saemix.res["omega"]
omega.null<-0*omega
#  diag(covariance.model)<-mydiag(saemix.model["covariance.model"])
#  omega<-0*saemix.res["omega"] # Why use only diag(omega) ???
#  diag(omega)<-mydiag(saemix.res["omega"])
hat.phi<-saemix.res["cond.mean.phi"]
nphi<-dim(hat.phi)[2]
nomega<-sum(covariance.model[lower.tri(covariance.model,diag=TRUE)])
if (saemixObject["model"]["modeltype"]=="structural"){
  nres<-length(saemix.res["indx.res"])
} else{
  nres <- 0
}
nytype<-length(unique(saemix.data["data"]["ytype"]))

############### Nodes

i1.omega2<-saemixObject["model"]["indx.omega"]
nphi1<-length(i1.omega2)
IOmega.phi1<-solve(omega[i1.omega2,i1.omega2,drop=FALSE])
mean.phi1<-saemix.res["mean.phi"][,i1.omega2,drop=FALSE]

y<-gqg.mlx(nphi1,nnodes.gq)
x<-(y$nodes-0.5)*2
w<-(y$weights)*(2**nphi1)

pres<-saemix.res["respar"]
cond.var.phi<-saemix.res["cond.var.phi"]
cond.mean.phi<-saemix.res["cond.mean.phi"]

############### Computing gradient of F
# 
nx<-dim(x)[1]
condsd.eta<-sqrt(cond.var.phi[,i1.omega2,drop=FALSE])
xmin<-cond.mean.phi[,i1.omega2,drop=FALSE]-nsd.gq*condsd.eta
xmax<-cond.mean.phi[,i1.omega2,drop=FALSE]+nsd.gq*condsd.eta
a<-(xmin+xmax)/2
b<-(xmax-xmin)/2



log.const<-0
idx.exp<-which(saemixObject["model"]["error.model"]=="exponential")
if(length(idx.exp)>0)
  #	if(saemixObject["model"]["error.model"]=="exponential")
  log.const<-(-sum(yobs[saemix.data["data"][,"ytype"] %in% idx.exp]))

Q<-0
for (j in 1:nx) {
  phi[,i1.omega2] <- a+b*matrix(rep(x[j,],saemix.data["N"]),ncol=nphi1,byrow=TRUE)
  psi<-transphi(phi,saemixObject["model"]["transform.par"])
  if(saemixObject["model"]["modeltype"]=="structural"){
    f<-saemixObject["model"]["model"](psi, saemix.data["data"][,"index"], xind)
    for(i in idx.exp) f[saemix.data["data"][,"ytype"]==i]<-log(cutoff(f[saemix.data["data"][,"ytype"]==i]))
    g<-error(f,pres,saemix.data["data"][,"ytype"])
    DYF[ind.io] <- -0.5*((yobs-f)/g)**2 - log(g)
    ly<-colSums(DYF)
  } else {
    f<-saemixObject["model"]["model"](psi, saemix.data["data"][,"index"], xind)
    DYF[ind.io] <- f
    ly<-colSums(DYF)
  }
  dphi1<-phi[,i1.omega2,drop=FALSE]-saemix.res["mean.phi"][,i1.omega2,drop=FALSE]
  lphi1<-(-0.5)*rowSums((dphi1%*%IOmega.phi1)*dphi1)
  ltot<-ly+lphi1
  ltot[is.na(ltot)]<-(-Inf)
  Q<-Q+w[j]*exp(ltot)
}



dphi<-cutoff(abs(colMeans(hat.phi))*1e-4,1e-10)
coefphi<-c(0,-1,1)
F<-array(data=0,dim=c(saemix.data["ntot.obs"],nphi,length(coefphi)))
gs<-matrix(0,saemix.data["ntot.obs"],4)
etype.exp<-which(saemix.model["error.model"]=='exponential')

for (l in 1:length(coefphi)) {
  for (j in 1:nphi) {
    phi<-hat.phi
    phi[,j]<-phi[,j]+coefphi[l]*dphi[j]
    psi<-transphi(phi,saemix.model["transform.par"])
    f <- saemix.model["model"](psi, saemix.data["data"][,"index"],xind)
    for(ityp in etype.exp) f[saemix.data["data"][,saemix.data["name.ytype"]]==ityp]<-log(cutoff(f[saemix.data["data"][,saemix.data["name.ytype"]]==ityp]))    
    F[,j,l]<-f
  }
}

ind.covariates<-which(saemix.model["betaest.model"]>0)
f0<-F[,1,1]
#  DF<-(F[,,3]-F[,,2])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi))/2 
DF<-(F[,,3]-F[,,1])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi)) #gradient of f (changed from F[,,2] to F[,,1])
z<-matrix(0,saemix.data["ntot.obs"],1)

# g0<-cutoff(saemix.res["respar"][1]+saemix.res["respar"][2]*abs(f0))
# if (saemixObject["model"]["modeltype"]=="structural")
#   g0<-error(f0,saemix.res@respar,saemix.data["data"][["ytype"]]) 

############### Parameters in the model, including betas and covariances
covariate.estim<-saemix.model["betaest.model"]
covariate.estim[1,]<-saemix.model["fixed.estim"]

j<-which(saemix.model["betaest.model"]>0)
ind.fixed.est<-(covariate.estim[j]>0)
npar<-sum(ind.fixed.est)
# Tracking indices for covariances
myidx.omega<-c()
myidx.cor<-c()
name.rand1<-name.rand2<-c()
myidx.track<-NULL
ipar<-npar
for(iom in 1:dim(covariance.model)[1]) {
  for(jom in iom:dim(covariance.model)[1]) {
    if(covariance.model[iom,jom]==1) {
      ipar<-ipar+1
      myidx.track<-rbind(myidx.track,c(ipar,iom,jom))
      if(iom==jom) {
        myidx.omega<-c(myidx.omega,ipar)
        name.rand1<-c(name.rand1,paste("Var",saemixObject@model@name.modpar[iom],sep="."))
        name.rand2<-c(name.rand2,paste("SD",saemixObject@model@name.modpar[iom],sep="."))
      }
      else {
        myidx.cor<-c(myidx.cor,ipar)
        name.rand1<-c(name.rand1,paste("Cov",saemixObject@model@name.modpar[iom],saemixObject@model@name.modpar[jom],sep="."))
        name.rand2<-c(name.rand2,paste("Corr",saemixObject@model@name.modpar[iom],saemixObject@model@name.modpar[jom],sep="."))
      }
    }
  }
}
if(length(myidx.cor)>0) {
  track.var<-myidx.track[myidx.track[,1] %in% myidx.omega,]
  for(i in myidx.cor) {
    ij<-which(myidx.track[,1]==i)
    myidx.track[ij,2]<-track.var[track.var[,2]==myidx.track[ij,2],1]
    myidx.track[ij,3]<-track.var[track.var[,2]==myidx.track[ij,3],1]
  }
}
if(saemixObject@model@modeltype=="structural")
  namallpar<-c(saemixObject@results@name.fixed,name.rand1, saemixObject@results@name.sigma[saemixObject@results@indx.res], name.rand2) else
    namallpar<-c(saemixObject@results@name.fixed,name.rand1, name.rand2)
