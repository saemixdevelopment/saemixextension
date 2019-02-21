context("Testing creation of SaemixModel for structural models\n")

test_that("Creating an empty SaemixModel object with class", {
  x<-new(Class="SaemixModel")
  expect_is(x,"SaemixModel")
  expect_null(body(x@model))
  expect_is(x@psi0,"matrix")
  expect_equal(dim(x@psi0),c(0,0))
  #  expect_error(new(Class="SaemixModel"));    # did not test for particular error message (not sure if good practice), but could be changed if desired
})

test_that("Errors in using saemixModel - no model", {
#  expect_error(x<-saemixModel());
  x<-saemixModel();
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixModel failed")
})

test_that("Errors in using saemixModel - no initial values for parameters", {
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
  x<-saemixModel(model=model1cpt);    
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixModel failed")
})

test_that("Minimal SaemixModel object", {
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
  x<-saemixModel(model=model1cpt, psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))), modeltype="structural")
  expect_is(x, "SaemixModel") # tests for particular class
  expect_equal(x@nb.parameters,3)
})

test_that("Successful creation of a SaemixModel object of type structural", {
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
  x<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  expect_is(x, "SaemixModel") # tests for particular class
  expect_equal(x@nb.parameters,3)
})

test_that("Successful creation of a SaemixModel object with 2 responses", {
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
  x<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))), name.response=c("PK","PD"),name.sigma=c("add","prop"), error.model=c("combined","proportional"))
  expect_is(x, "SaemixModel") # tests for particular class
  expect_equal(x@nb.parameters,3)
  expect_equal(length(x@indx.res),3)
})

test_that("Successful creation of a SaemixModel object with all arguments", {
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
  x<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))), transform.par=c(1,1,1), covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="combined",error.init=c(1,0.5))
  expect_is(x, "SaemixModel") # tests for particular class
  expect_equal(x@nb.parameters,3)
  expect_identical(colnames(x@omega.init),colnames(x@psi0))
  expect_equal(dim(x@omega.init),c(3,3))
  expect_identical(x@error.init,c(1,0.5))
  expect_identical(sum(x@covariate.model),1)
  expect_identical(sum(x@betaest.model),4)
  expect_identical(sum(x@transform.par),3)
  expect_identical(sum(x@fixed.estim),3)
})
