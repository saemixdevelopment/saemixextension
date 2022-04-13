context("Structural model - object SaemixStructuralModel\n")

test_that("Creating an empty SaemixModel object with class (just creates the slots)", {
  x<-new(Class="SaemixStructuralModel")
  expect_error(validObject(x))
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
  expect_equal(x@name.outcome, "y")
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
  out1<-createSaemixOutcome(continuousOutcome(model='combined1'))
  x<-new(Class="SaemixStructuralModel", model=model1cpt, outcome=list(concentration=out1))
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
         outcome=list(concentration=createSaemixOutcome(continuousOutcome(model='combined1')), pca=createSaemixOutcome(continuousOutcome(model="proportional"))))
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
  out1<-createSaemixOutcome(continuousOutcome(model='combined1'))
  out2<-createSaemixOutcome(discreteOutcome(type="event"))
  out3<-createSaemixOutcome(discreteOutcome())
  x<-new(Class="SaemixStructuralModel", model=model1cptdirect, 
         outcome=list(concentration=out1, bleeding=out2,pain=out3))
  expect_equal(x@nb.outcome, 3)
  expect_equal(x@name.outcome, c("concentration","bleeding","pain"))
  expect_equal(x@outcome[[1]]@error.model, "combined1")
  expect_equal(x@outcome[[1]]@error.npar, 2)
  expect_equal(x@outcome[[1]]@type, "continuous")
  expect_equal(x@outcome[[2]]@type, "event")
  expect_equal(x@outcome[[3]]@type, "binary")
})


context("Class SaemixParameterModel \n")

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
  out1<-createSaemixOutcome(continuousOutcome(model='combined1'))
  out2<-createSaemixOutcome(discreteOutcome(type="event"))
  out3<-createSaemixOutcome(discreteOutcome())
  
  x<-new(Class="SaemixParameterModel", model=model1cptdirect, 
         outcome=list(concentration=out1, bleeding=out2,pain=out3), parameter=colnames(psi0))
  expect_equal(x@name.outcome, c("concentration","bleeding","pain"))
  expect_equal(x@npar, 4)
  expect_equal(sum(x@transform.par), 0)
  expect_equal(x@name.modpar, colnames(psi0))
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
  out1<-createSaemixOutcome(continuousOutcome(model='combined1'))
  out2<-createSaemixOutcome(discreteOutcome(type="event"))
  out3<-createSaemixOutcome(discreteOutcome())
  
  x<-new(Class="SaemixParameterModel", model=model1cptdirect, 
         outcome=list(concentration=out1, bleeding=out2,pain=out3), parameter=colnames(psi0), mu.start=psi0[1,])
  expect_equal(x@name.outcome, c("concentration","bleeding","pain"))
  expect_equal(x@npar, 4)
  expect_equal(sum(x@transform.par), 0)
  expect_equal(x@name.modpar, colnames(psi0))
  expect_equal(x@mu.start, psi0[1,])
})


context("Class SaemixModel \n")

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
  out1<-createSaemixOutcome(continuousOutcome(model='combined1'))
  out2<-createSaemixOutcome(discreteOutcome(type="event"))
  out3<-createSaemixOutcome(discreteOutcome())
  
  x<-new(Class="SaemixModel", model=model1cptdirect, 
         outcome=list(concentration=out1, bleeding=out2,pain=out3), parameter=colnames(psi0))
  expect_equal(x@name.outcome, c("concentration","bleeding","pain"))
  expect_equal(x@npar, 4)
  expect_equal(sum(x@transform.par), 0)
  expect_equal(sum(x@var.model[[1]]@omega), 4)
  expect_equal(x@name.modpar, colnames(psi0))
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
  out1<-createSaemixOutcome(continuousOutcome(model='combined1'))

  x<-new(Class="SaemixModel", model=model1cpt, outcome=list(concentration=out1), parameter=c("ka","V","CL"))
  expect_equal(x@nb.outcome, 1)
  expect_equal(x@name.outcome, "concentration")
  expect_equal(x@outcome[[1]]@error.model, "combined1")
  expect_equal(x@outcome[[1]]@error.npar, 2)
  expect_equal(x@outcome[[1]]@type, "continuous")
  expect_equal(x@npar, 3)
  expect_equal(sum(x@transform.par), 0)
  expect_equal(x@var.model[[1]]@omega.model, diag(3))
  expect_equal(x@var.model[[1]]@name.level, "iiv")
  expect_equal(x@var.model[[1]]@variable, "id")
})

test_that("Testing class with one variability level", {
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
  out1<-createSaemixOutcome(continuousOutcome(model='combined1'))
  var1<-saemixVarModel(size=3, omega=diag(c(0,0.5,0.5)), omega.model=diag(c(0,1,1)))
  
  x<-new(Class="SaemixModel", model=model1cpt, outcome=list(concentration=out1), parameter=c("ka","V","CL"), var.model=var1)
  expect_equal(x@nb.outcome, 1)
  expect_equal(x@name.outcome, "concentration")
  expect_equal(x@outcome[[1]]@error.model, "combined1")
  expect_equal(x@outcome[[1]]@error.npar, 2)
  expect_equal(x@outcome[[1]]@type, "continuous")
  expect_equal(x@npar, 3)
  expect_equal(sum(x@transform.par), 0)
  expect_equal(x@var.model[[1]]@omega.model, diag(c(0,1,1)))
  expect_equal(sum(x@var.model[[1]]@omega.model), 2)
  expect_equal(x@var.model[[1]]@name.level, "iiv")
  expect_equal(x@var.model[[1]]@variable, "id")
  expect_equal(names(x@var.model)[1], "iiv")
})

test_that("Testing class specifying most of the components", {
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
  out1<-createSaemixOutcome(continuousOutcome(model='combined1', start=c(1,0.5)))
  out2<-createSaemixOutcome(continuousOutcome())
  mat1<-diag(c(0,1,1,1))
  mat1[2,3]<-mat1[3,2]<-1
  var1<-saemixVarModel(size=3, omega=diag(c(0,0.5,0.5,1)), omega.model=mat1, omega.model.fix=diag(c(0,0,0,1)))
  covmod<-matrix(c(0,1,1,0,0,0,0,1), ncol=4, byrow=T, dimnames=list(c("wt","sex"), NULL))
  
  x<-new(Class="SaemixModel", model=model1cptdirect, parameter=c("ka","V","CL","IC50"),
         outcome=list(concentration=out1, effect=out2),  
         mu.start = c(1.,20,0.5, 5), transform.par=c(1,1,1,1), mu.fix=c(1,0,0,0),
         covariate.model=covmod, covariate.model.fix=matrix(c(rep(0,7),1), ncol=4, byrow=T), beta.start=c(1,0.75),
         var.model=var1, verbose=FALSE)

  expect_equal(x@nb.outcome, 2)
  expect_equal(x@name.outcome, c("concentration","effect"))
  expect_equal(x@outcome[[1]]@error.model, "combined1")
  expect_equal(x@outcome[[1]]@error.npar, 2)
  expect_equal(x@outcome[[1]]@type, "continuous")
  expect_equal(x@npar, 4)
  expect_equal(sum(x@transform.par), 4)
  expect_equal(x@var.model[[1]]@omega.model, mat1)
  expect_equal(sum(x@var.model[[1]]@omega.model), 5)
  expect_equal(x@var.model[[1]]@name.level, "iiv")
  expect_equal(x@var.model[[1]]@variable, "id")
  expect_equal(names(x@var.model)[1], "iiv")
  expect_equal(sum(x@var.model[[1]]@omega.model.fix), 1)
  expect_equal(sum(x@covariate.model), 3)
  expect_equal(sum(x@covariate.model.fix), 1)
})


test_that("Testing class with two variability levels", {
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
  out1<-createSaemixOutcome(continuousOutcome(model='combined1'))
  mat1<-diag(c(0,1,1))
  mat1[2,3]<-mat1[3,2]<-1
  var1<-saemixVarModel(size=3, omega=diag(c(0,0.5,0.5)), omega.model=mat1)
  var2<-saemixVarModel(size=3, name.level="iov", variable="occ", omega=diag(c(0,0,0.3)), omega.model=diag(c(0,0,1)))
  
  x<-new(Class="SaemixModel", model=model1cpt, outcome=list(concentration=out1), parameter=c("ka","V","CL"), var.model=list(var1, var2))
  expect_equal(x@nb.outcome, 1)
  expect_equal(x@name.outcome, "concentration")
  expect_equal(x@outcome[[1]]@error.model, "combined1")
  expect_equal(x@outcome[[1]]@error.npar, 2)
  expect_equal(x@outcome[[1]]@type, "continuous")
  expect_equal(x@npar, 3)
  expect_equal(sum(x@mu.start), 3)
  expect_equal(sum(x@transform.par), 0)

  expect_equal(x@nvarlevel, 2)
  expect_equal(sum(x@var.model[[1]]@omega.model), 4)
  expect_equal(sum(x@var.model[[2]]@omega.model), 1)
  expect_equal(x@var.model[[1]]@name.level, "iiv")
  expect_equal(x@var.model[[1]]@variable, "id")
  expect_equal(names(x@var.model)[1], "iiv")
  expect_equal(x@var.model[[2]]@name.level, "iov")
  expect_equal(x@var.model[[2]]@variable, "occ")
})


context("Auxiliary functions for SaemixModel objects - covariance model\n")

context("Using saemixModel() to create model\n")

test_that("Missing elements, model or parameters", {
  expect_message(x<-saemixModel(verbose=TRUE))
  expect_equal(x,"Creation of SaemixModel failed")
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
  expect_message(x<-saemixModel(model=model1cpt, verbose=TRUE))
  expect_equal(x,"Creation of SaemixModel failed")
})


test_that("Minimal creation, only model and parameter, check defaults", {
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
  x1<-saemixModel(model=model1cpt, parameter=c("ka","V","CL"))
  x2<-saemixModel(model=model1cpt, parameter=c(ka=1, V=20, CL=0.5))
  expect_equal(x1@npar, 3)
  expect_equal(x1@name.outcome, "y")
  expect_equal(x1@mu.start,rep(1,3))
  expect_equal(x1@mu.fix,rep(0,3))
  expect_equal(x1@transform.par,rep(0,3))
  expect_equal(x1@var.model[[1]]@omega, diag(3))
  expect_equal(x1@var.model[[1]]@omega.model, diag(3))
  expect_equal(x1@covariate.model, matrix(data=0,nrow=0, ncol=0))
  expect_equal(x1@covariate.model.fix, matrix(data=0,nrow=0, ncol=0))
  expect_equal(x1@name.modpar,c("ka","V","CL"))
  expect_equal(length(x1@name.thetas),0)
  expect_equal(length(x1@name.X),0)
  expect_equal(length(x1@name.cov),0)
  expect_equal(length(x1@name.predictors),0)
  expect_equal(x2@npar, 3)
  expect_equal(x2@mu.start,c(1,20,0.5))
  expect_equal(x2@transform.par,rep(0,3))
})

test_that("Full model", {
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
  out1<-createSaemixOutcome(continuousOutcome(model='combined1', start=c(1,0.5)))
  out2<-createSaemixOutcome(continuousOutcome())
  mat1<-diag(c(0,1,1,1))
  mat1[2,3]<-mat1[3,2]<-1
  var1<-saemixVarModel(size=3, omega=diag(c(0,0.5,0.5,1)), omega.model=mat1, omega.model.fix=diag(c(0,0,0,1)))
  covmod<-matrix(c(0,1,1,0,0,0,0,1), ncol=4, byrow=T, dimnames=list(c("wt","sex"), NULL))
  
  x<-saemixModel(model=model1cptdirect, description="One compartment model with direct Imax effect",
                 outcome=list(concentration=out1, effect=out2), 
                 parameter=c(ka=1, V=20, CL=0.5, IC50=5), transform.par=c(1,1,1,1), mu.fix=c(1,0,0,0),
                 covariate.model=covmod, covariate.model.fix=matrix(c(rep(0,7),1), ncol=4, byrow=T), beta.start=c(1,0.75),
                 var.model=var1, verbose=FALSE)
  expect_equal(x@nb.outcome, 2)
  expect_equal(x@name.outcome, c("concentration","effect"))
  expect_equal(x@outcome[[1]]@error.model, "combined1")
  expect_equal(x@outcome[[1]]@error.npar, 2)
  expect_equal(x@outcome[[1]]@type, "continuous")
  expect_equal(x@npar, 4)
  expect_equal(sum(x@transform.par), 4)
  expect_equal(x@var.model[[1]]@omega.model, mat1)
  expect_equal(sum(x@var.model[[1]]@omega.model), 5)
  expect_equal(x@var.model[[1]]@name.level, "iiv")
  expect_equal(x@var.model[[1]]@variable, "id")
  expect_equal(names(x@var.model)[1], "iiv")
  expect_equal(sum(x@var.model[[1]]@omega.model.fix), 1)
  expect_equal(sum(x@covariate.model), 3)
  expect_equal(sum(x@covariate.model.fix), 1)
})

test_that("Full model, mismatch between covariate model size (3 columns) and nb of parameters (4)", {
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
  out1<-createSaemixOutcome(continuousOutcome(model='combined1', start=c(1,0.5)))
  out2<-createSaemixOutcome(continuousOutcome())
  mat1<-diag(c(0,1,1,1))
  mat1[2,3]<-mat1[3,2]<-1
  var1<-saemixVarModel(size=3, omega=diag(c(0,0.5,0.5,1)), omega.model=mat1, omega.model.fix=diag(c(0,0,0,1)))
  covmod<-matrix(c(0,1,1,0,0,1), ncol=3, byrow=T, dimnames=list(c("wt","sex"), NULL))
  
  expect_message(x<-saemixModel(model=model1cptdirect, description="One compartment model with direct Imax effect",
                 outcome=list(concentration=out1, effect=out2), 
                 parameter=c(ka=1, V=20, CL=0.5, IC50=5), transform.par=c(1,1,1,1), mu.fix=c(1,0,0,0),
                 covariate.model=covmod, covariate.model.fix=matrix(c(rep(0,5),1), ncol=3, byrow=T), beta.start=c(1,0.75),
                 var.model=var1, verbose=TRUE))
  expect_equal(x@nb.outcome, 2)
  expect_equal(x@name.outcome, c("concentration","effect"))
  expect_equal(x@outcome[[1]]@error.model, "combined1")
  expect_equal(x@outcome[[1]]@error.npar, 2)
  expect_equal(x@outcome[[1]]@type, "continuous")
  expect_equal(x@npar, 4)
  expect_equal(sum(x@transform.par), 4)
  expect_equal(x@var.model[[1]]@omega.model, mat1)
  expect_equal(sum(x@var.model[[1]]@omega.model), 5)
  expect_equal(x@var.model[[1]]@name.level, "iiv")
  expect_equal(x@var.model[[1]]@variable, "id")
  expect_equal(names(x@var.model)[1], "iiv")
  expect_equal(sum(x@var.model[[1]]@omega.model.fix), 1)
  expect_equal(length(x@covariate.model), 0)
  expect_equal(length(x@covariate.model.fix), 0)
})

test_that("Legacy saemix from 3.0 - check same as new saemixModel", {
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
  x1<-saemixModel(model=model1cpt, c(ka=1, V=20, CL=0.5, IC50=5))
  x2<-saemixModel.legacy(model=model1cpt,psi0=matrix(c(1, 20, 0.5), ncol=3, dimnames=list(NULL,c("ka","V","CL"))))
  expect_identical(x1, x2)
})


