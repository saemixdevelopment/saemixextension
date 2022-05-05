
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
  lout<-list(out1=saemixOutcome(model='combined1', start=c(1,0.5)), 
             out2=saemixOutcome()) 
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
  x1<-saemixModel(model=model1cpt, parameter=c(ka=1, V=20, CL=0.5, IC50=5))
  x2<-saemixModel.legacy(model=model1cpt,psi0=matrix(c(1, 20, 0.5), ncol=3, dimnames=list(NULL,c("ka","V","CL"))))
  expect_equal(x1, x2)
  #  expect_identical(x1, x2)
})

