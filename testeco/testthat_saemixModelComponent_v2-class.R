context("Structural model - object SaemixStructuralModel\n")

test_that("Creating an empty SaemixModel object with class (just creates the slots)", {
  x<-new(Class="SaemixStructuralModel")
  expect_is(x,"SaemixStructuralModel")
  expect_null(body(x@model))
  #  expect_error(new(Class="SaemixModel"));    # did not test for particular error message (not sure if good practice), but could be changed if desired
})

test_that("Testing class with a model", {
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
  x<-new(Class="SaemixStructuralModel", model=model1cpt)
  expect_equal(x@nb.outcome, 1)
  expect_equal(x@name.outcome, "y1")
  expect_equal(x@outcome[[1]]@error.model, "constant")
  expect_equal(x@outcome[[1]]@error.parameters, 1)
  expect_equal(x@outcome[[1]]@type, "continuous")
})

test_that("Testing class with a wrong model (wrong arguments)", {
  model1cpt.pb<-function(tim, dose, ka, CL, V) { 
    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  expect_error(x<-new(Class="SaemixStructuralModel", model=model1cpt.pb))
})

test_that("Testing class with a model, explicit outcome", {
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
  x<-new(Class="SaemixStructuralModel", model=model1cpt, outcome=list(concentration=continuousOutcome(model='combined1')))
  expect_equal(x@nb.outcome, 1)
  expect_equal(x@name.outcome, "concentration")
  expect_equal(x@outcome[[1]]@error.model, "combined1")
  expect_equal(x@outcome[[1]]@error.npar, 2)
  expect_equal(x@outcome[[1]]@type, "continuous")
})


test_that("Testing class with a model, two continuous outcome", {
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
  x<-new(Class="SaemixStructuralModel", model=model1cptdirect, 
         outcome=list(concentration=continuousOutcome(model='combined1'), pca=continuousOutcome(model="proportional")))
  expect_equal(x@nb.outcome, 2)
  expect_equal(x@name.outcome, c("concentration","pca"))
  expect_equal(x@outcome[[1]]@error.model, "combined1")
  expect_equal(x@outcome[[1]]@error.npar, 2)
  expect_equal(x@outcome[[1]]@type, "continuous")
})

test_that("Testing class with a model (non functional), one continuous outcome + two discrete outcomes", {
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
  x<-new(Class="SaemixStructuralModel", model=model1cptdirect, 
         outcome=list(concentration=continuousOutcome(model='combined1'), bleeding=discreteOutcome(type="event"),
                      pain=discreteOutcome()))
  expect_equal(x@nb.outcome, 3)
  expect_equal(x@name.outcome, c("concentration","bleeding","pain"))
  expect_equal(x@outcome[[1]]@error.model, "combined1")
  expect_equal(x@outcome[[1]]@error.npar, 2)
  expect_equal(x@outcome[[1]]@type, "continuous")
  expect_equal(x@outcome[[2]]@type, "event")
  expect_equal(x@outcome[[3]]@type, "binary")
})

test_that("Testing class with a model (non functional), one continuous outcome + two discrete outcomes", {
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
  psi0 <- matrix(c(1.,20,0.5, 5), ncol=4,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL","IC50")))
  
  x<-new(Class="SaemixModel", model=model1cptdirect, 
         outcome=list(concentration=continuousOutcome(model='combined1'), bleeding=discreteOutcome(type="event"),
                      pain=discreteOutcome()),
         psi0=psi0)
  expect_equal(x@name.outcome, c("concentration","bleeding","pain"))
})

context("Class SaemixModel \n")

test_that("Testing class with a model, explicit outcome", {
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
  x<-new(Class="SaemixModel", model=model1cpt, outcome=list(concentration=continuousOutcome(model='combined1')), parameter=c("ka","V","CL"))
  expect_equal(x@nb.outcome, 1)
  expect_equal(x@name.outcome, "concentration")
  expect_equal(x@outcome[[1]]@error.model, "combined1")
  expect_equal(x@outcome[[1]]@error.npar, 2)
  expect_equal(x@outcome[[1]]@type, "continuous")
  expect_equal(x@npar, 3)
  expect_equal(sum(x@transform.par), 0)
  expect_equal(x@var.level[[1]]@omega.model, diag(3))
  expect_equal(x@var.level[[1]]@name.level, "iiv")
  expect_equal(x@var.level[[1]]@variable, "id")
})


context("Auxiliary functions for SaemixModel objects - covariance model\n")
