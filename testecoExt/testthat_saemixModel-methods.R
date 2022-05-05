
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
  expect_equal(x1@name.outcome, "y1")
  expect_equal(x1@mu.start,rep(1,3))
  expect_equal(x1@mu.fix,rep(0,3))
  expect_equal(x1@transform.par,rep(1,3))
  expect_equal(unname(x1@var.model[[1]]@omega), diag(3))
  expect_equal(unname(x1@var.model[[1]]@omega.model), diag(3))
  expect_equal(x1@covariate.model, matrix(data=0,nrow=0, ncol=0))
  expect_equal(x1@covariate.model.fix, matrix(data=0,nrow=0, ncol=0))
  expect_equal(x1@name.modpar,c("ka","V","CL"))
  expect_equal(length(x1@name.thetas),0)
  expect_equal(length(x1@name.X),0)
  expect_equal(length(x1@name.cov),0)
  expect_equal(length(x1@name.predictors),0)
  expect_equal(x2@npar, 3)
  expect_equal(unname(x2@mu.start),c(1,20,0.5))
  expect_equal(x2@transform.par,rep(1,3))
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
  lout<-list(out1=saemixOutcome(model='combined1', start=c(1,0.5)), 
             out2=saemixOutcome())
  lpar<-list(ka=saemixPar(omega.level=c()), 
             V=saemixPar(mu.start=20, covariate=c(wt=contCov(beta=1, beta.fix=1))), 
             CL=saemixPar(mu.start=0.5, rho.param=c("V"), covariate=c(wt=contCov(beta=0.75, beta.fix=1))), 
             IC50=saemixPar(mu.start=5,  covariate=c(sex=binCov())))
#  getCovariateModel(lpar)
  x<-saemixModel(model=model1cptdirect, description="One compartment model with direct Imax effect",
                 outcome=list(concentration=saemixOutcome(model='combined1', start=c(1,0.5)), effect=saemixOutcome()), 
                 parameter=lpar, verbose=FALSE)
  
  mat1<-diag(c(0,1,1,1))
  mat1[2,3]<-mat1[3,2]<-1
  covmod<-matrix(c(0,1,1,0,0,0,0,1), ncol=4, byrow=T, dimnames=list(c("wt","sex"), NULL))
  
  expect_equal(x@nb.outcome, 2)
  expect_equal(x@name.outcome, c("concentration","effect"))
  expect_equal(x@outcome[[1]]@error.model, "combined1")
  expect_equal(x@outcome[[1]]@error.npar, 2)
  expect_equal(x@outcome[[1]]@type, "continuous")
  expect_equal(x@npar, 4)
  expect_equal(sum(x@transform.par), 4)
  expect_equal(unname(x@var.model[[1]]@omega.model), mat1)
  expect_equal(sum(x@var.model[[1]]@omega.model), 5)
  expect_equal(x@var.model[[1]]@name.level, "iiv")
  expect_equal(x@var.model[[1]]@variable, "id")
  expect_equal(names(x@var.model)[1], "iiv")
  expect_equal(sum(x@var.model[[1]]@omega.model.fix), 0)
  expect_equal(sum(x@covariate.model), 3)
  expect_equal(unname(x@covariate.model), unname(covmod))
  expect_equal(sum(x@covariate.model.fix), 2)
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
  lout<-list(out1=saemixOutcome(model='combined1', start=c(1,0.5)), 
             out2=saemixOutcome())
  lpar<-list(ka=saemixPar(omega.level=c()), 
             V=saemixPar(mu.start=20, covariate=c(wt=contCov(beta=1, beta.fix=1))), 
             CL=saemixPar(mu.start=0.5, rho.param=c("V"), covariate=c(wt=contCov(beta=0.75, beta.fix=1))), 
             IC50=saemixPar(mu.start=5,  covariate=c(sex=binCov())))
  
  
  x1<-saemixModel(model=model1cpt, parameter=c(ka=1, V=20, CL=0.5), outcome=c("out1"))
  x2<-saemixModel.legacy(model=model1cpt,psi0=matrix(c(1, 20, 0.5), ncol=3, dimnames=list(NULL,c("ka","V","CL"))))
  expect_equal(x1, x2)
  #  expect_identical(x1, x2)
})

