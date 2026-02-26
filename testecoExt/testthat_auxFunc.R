

context("Computing error and ssq")
test_that("Testing function error", {
  lerr <- list(new(Class="SaemixErrorModel"))
  fpred<-1:10
  gpred <- error(fpred,lerr,rep(1,length(fpred)))
  expect_equal(rep(1,10),gpred)
  lerr <- list(new(Class="SaemixErrorModel", name="proportional", param=0.2))
  gpred <- error(fpred,lerr,rep(1,length(fpred)))
  expect_equal(gpred,0.2*fpred)
})

test_that("Testing function ssq", {
  xerr <- new(Class="SaemixErrorModel")
  fpred<-1:10
  yobs<-fpred+0.2
  xcal <- ssq(xerr@param, yobs, fpred, xerr)
  expect_equal(xcal, 10*(.2**2))
  xerr <- new(Class="SaemixErrorModel", name="proportional", param=0.2)
  xcal <- ssq(xerr@param, yobs, fpred, xerr)
  gpred <- error(fpred,list(xerr),rep(1,length(fpred)))
  expect_equal(xcal, sum((yobs-fpred)**2/(gpred**2) + 2*log(gpred)))
})

context("Testing computational functions")

############################################ Recoding func_aux.R
# Generate design matrices at each varlevel
## dans data ou dans une fonction qui associe data et model ?

# Old version of lists
if(FALSE) {
  Dargs.old<-list(IdM=IdM, XM=XM, yM=yM, NM=NM, N=N, nobs=saemix.data["ntot.obs"],
              yobs=saemix.data["data"][,saemix.data["name.response"]],transform.par=saemix.model["transform.par"],
              error.model=saemix.model["error.model"],structural.model=structural.model , 
              etype.exp=which(saemix.model["error.model"] == "exponential"),modeltype=saemix.model["modeltype"])
  
  Uargs.old<-list(nchains=saemix.options$nb.chains,nb.parameters=nb.parameters, nb.betas=nb.betas, nb.etas=nb.etas, 
              nb.parest=nb.parest,indx.betaC=indx.betaC, indx.betaI=indx.betaI, ind.res=ind.res,indest.omega=indest.omega,
              i0.omega2=i0.omega2, i1.omega2=i1.omega2,j.covariate=j.covariate, j0.covariate=j0.covariate,
              ind.fix10=ind.fix10, ind.fix11=ind.fix11, ind.fix1=ind.fix1, ind.fix0=ind.fix0,
              MCOV0=MCOV0, COV=COV, COV0=COV0, COV1=COV1, LCOV=LCOV, COV2=COV2, dstatCOV=dstatCOV, 
              Mcovariates=Mcovariates, ind.ioM=ind.ioM)
  
}
# New version of lists
if(FALSE) {
  # exponential and modeltype are now contained in the list of error.model
  Dargs <- list(IdM=IdM, XM=XM, yM=yM, NM=NM, N=N, nobs=saemix.data["ntot.obs"],
                yobs=saemix.data["data"][,saemix.data["name.response"]],
                transform.par = saemix.model@transform,
                structural.model=structural.model, 
                error.model = list.error.model)
}

# WiP
## define elements needed from Dargs and args
## change to Dargs (data), Margs (model), and potentially other args (Vargs for variances ?)
context("Computing LL_y")

datDir <- "/home/eco/work/saemix/saemixextension/data"

test_that("Testing compute.LLy, IIV only - computing likelihood of predicted (not simulated) data, expecting roughly 0", {
  theonocov<-read.csv(file.path(datDir,"../","data40","theoSimul_nocov.csv"),header=T,na=".")
  x<-saemixData(name.data=theonocov, name.predictors=c("tim","dose"), name.response="ypred",varlevel=c("id","occ")) 
  psinocov<-read.csv(file.path(datDir,"../","data40","theoSimul_psinocov.csv"),header=T,na=".")
  # Data related list
  dlist <- createStructData(x, nb.chains=1)
  Dargs <- dlist$Dargs
  DYF <- dlist$DYF
  # Defining matching model
  parameters<-list(ka=saemixParam(mu.init=1.5),
                   vd=saemixParam(mu.init=35, covariate="wt", covariate.init=c(1), covariate.estim=c("fixed")),  
                   cl=saemixParam(name="cl",mu.init=1.5, corr = list(iiv=c("vd")), covariate=c("wt","sex"), covariate.init=c(0.75,0), covariate.estim=c("fixed","estimated"), corr.init=list(iiv=c(0.5))))
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
  xmodel <- saemixModel(model=theo1cpt, parameters=parameters)
  mlist <- createStructModel(xmodel)
  Margs <- mlist$Margs
  psiM<-as.matrix(psinocov[,2:4])
  phiM<-transformPar(psiM, Margs$invtransform)
  LLy <- compute.LLy(phiM, Margs, Dargs, DYF) 
  expect_equal(round(sum(LLy),digits=2),0)
})

test_that("Testing compute.LLy, IIV only", {
  theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T,na=".")
  x<-saemixData(name.data=theo.saemix, name.group="wrongid",name.predictors=c("Dose","Time"), name.response="Concentration")
  # Data related list
  dlist <- createStructData(x, nb.chains=1)
  Dargs <- dlist$Dargs
  DYF <- dlist$DYF
  expect_true(max(summary(Dargs$IdM-x@var.index[[1]]))<10-6)
  
  parameters<-list(ka=saemixParam(),
               vd=saemixParam(mu.init=20, covariate="wt", covariate.init=c(1), covariate.estim=c("fixed")),  
               cl=saemixParam(name="cl",mu.init=2, corr = list(iiv=c("vd")), covariate=c("wt","sex"), covariate.init=c(0.75,0), covariate.estim=c("fixed","estimated"), corr.init=list(iiv=c(0.5))))
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
  theomodel.iiv <- saemixModel(model=theo1cpt, parameters=parameters)
  mlist <- createStructModel(theomodel.iiv)
  Margs <- mlist$Margs
  phiM<-as.matrix(data.frame(ka=log(seq(1,2,length.out=x@N)), V=log(seq(25,35,length.out=x@N)), CL=log(seq(1.7,1.3,length.out=x@N))))
  LLy <- compute.LLy(phiM, Margs, Dargs, xstr$DYF) 
  expect_equal(round(sum(LLy),digits=2),2162.2)
})
