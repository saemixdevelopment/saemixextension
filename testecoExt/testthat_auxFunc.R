
context("Testing computational functions")

test_that("Converting from phi to psi and back using transphi", {
  param2<-list(ka=saemixParam(mu.init=c(2), sd.init=c(0.5), covariate=c("lage"),covariate.init=c(0.2)),
               vd=saemixParam(mu.init=c(20), sd.init=0.7, covariate="lwt", covariate.init=c(1), covariate.estim=c("fixed")))
  indivmodel <- new(Class="SaemixIndivModel", param2)
  nsuj<-10
  cdesign <- matrix(c(rep(1,nsuj), log(seq(from=50,to=(50+2*(nsuj-1)), by=2)/60),
                      log(seq(from=90, to=(90-4*(nsuj-1)), by=-4)/70)), ncol=3)
  colnames(cdesign)<-c("pop","lage","lwt")
  phipop <- cdesign %*% indivmodel@popmodel[[1]]@phi
  psipop<- matrix(c(indivmodel@transform[[1]](phipop[,1]),indivmodel@transform[[2]](phipop[,2])), ncol=2)
  
  expect_equivalent(transphi(phipop, indivmodel@transform), psipop)
  expect_equivalent(transphi(psipop, indivmodel@invtransform), phipop)
})


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

context("Computing error and ssq")
test_that("Testing function error", {
  lerr <- list(new(Class="SaemixErrorModel"))
  fpred<-1:10
  gpred <- error(fpred,lerr,rep(1,length(fpred)))
  expect_equal(fpred,gpred)
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


# WiP
## define elements needed from Dargs and args
context("Computing LL_y")

test_that("Testing compute.LLy", {
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
  saemix.model <- new(Class="SaemixModel", parameters=parameters, model=model1cpt)
  
  IdM
  XM
  yM
  ind.ioM
  
  list.error.model<-vector(mode="list",length=saemix.model@noutcome)
  for(i in 1:saemix.model@noutcome)
    if(saemix.model@outcome[[i]]@type.outcome=="continuous") list.error.model[[i]]<-saemix.model@outcome[[i]]@error.model else list.error.model[[i]]<-"none"
    
    Dargs <- list(IdM=IdM, XM=XM, yM=yM,
                  transform.par = saemix.model@transform,
                  structural.model=saemix.model@model, 
                  error.model = list.error.model)
    args<-list(ind.ioM=ind.ioM)
    
})
