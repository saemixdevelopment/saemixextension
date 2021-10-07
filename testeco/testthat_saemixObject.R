context("Default control options\n")

test_that("Default option list", {
  #  expect_error(x<-saemixModel());
  xopt<-saemixControl();
  expect_is(xopt,"list")
  expect_equal(xopt$map,TRUE)
  expect_equal(xopt$fim,TRUE)
  expect_equal(xopt$ll.is,TRUE)
  expect_equal(xopt$ll.gq,FALSE)
  expect_equal(xopt$nb.simpred,100)
})

test_that("Option list - changes", {
  xopt<-saemixControl(nb.chains=5,nbiter.saemix = c(500,300), ipar.lmcmc = 100)
  expect_is(xopt,"list")
  expect_equal(xopt$nb.chains,5)
  expect_equal(xopt$nbiter.saemix,c(500,300))
  expect_equal(xopt$ipar.lmcmc,100)
  expect_equal(xopt$nbiter.sa,250)
})

context("Creating an empty SaemixObject object\n")

test_that("Create an SaemixObject with an empty results slot", {
  smx.data<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.group=c("Id"),name.predictors=c("Dose","Time"),name.covariates=c("Weight","Sex"), name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time",verbose=F)
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
  smx.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))), transform.par=c(1,1,1), covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="combined",error.init=c(1,0.5))
  smx.opt<-saemixControl(nb.chains=5,nbiter.saemix = c(500,300), ipar.lmcmc = 100)
  x<-createSaemixObject.empty(smx.model,smx.data,smx.opt)
  expect_equal(x@results@status,"empty")
  expect_length(x@results@fixed.effects,0)
  expect_length(x@results@name.fixed,0)
})

context("Initialising an SaemixObject object with initial population estimates\n")

test_that("Fill in the results of an SaemixObject with initial estimates", {
  smx.data<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.group=c("Id"),name.predictors=c("Dose","Time"),name.covariates=c("Weight","Sex"), name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time",verbose=F)
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
  smx.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))), transform.par=c(1,1,1), covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="combined",error.init=c(1,0.5))
  smx.opt<-saemixControl(nb.chains=5,nbiter.saemix = c(500,300), ipar.lmcmc = 100)
  x<-createSaemixObject.initial(smx.model,smx.data,smx.opt)
  expect_equal(x@results@status,"initial")
  expect_equal(x@results@name.fixed,x@model@name.fixed)
  expect_equivalent(x@results@fixed.effects,c(1.,20,0.5,-0.01))
  expect_equivalent(x@results@omega,diag(c(1,1,1)))
})

